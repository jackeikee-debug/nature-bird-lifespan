"""Run NCBI Gene symbol/taxid pilot mapping for Phase 2 targeted queue."""

from __future__ import annotations

import argparse
import json
import pathlib
import time
import urllib.parse
import urllib.request
from collections import Counter

import pandas as pd


ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def cache_key(symbol: str, taxid: str) -> str:
    safe_symbol = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in symbol)
    return f"taxid_{taxid}__symbol_{safe_symbol}.json"


def fetch_json(url: str, params: dict[str, str], retries: int = 4) -> dict:
    full_url = f"{url}?{urllib.parse.urlencode(params)}"
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(full_url, timeout=30) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            time.sleep(2 + attempt * 2)
    raise RuntimeError(f"NCBI request failed: {last_error}")


def lookup(symbol: str, taxid: str, cache_dir: pathlib.Path, email: str, sleep_s: float) -> dict:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / cache_key(symbol, taxid)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))

    term = f'("{symbol}"[Gene Name] OR "{symbol}"[Preferred Symbol]) AND txid{taxid}[Organism]'
    params = {
        "db": "gene",
        "term": term,
        "retmode": "json",
        "retmax": "10",
        "tool": "bird_lifespan_phase2",
    }
    if email:
        params["email"] = email
    search = fetch_json(ESEARCH_URL, params)
    ids = search.get("esearchresult", {}).get("idlist", [])
    summary = {}
    if ids:
        sparams = {
            "db": "gene",
            "id": ",".join(ids),
            "retmode": "json",
            "tool": "bird_lifespan_phase2",
        }
        if email:
            sparams["email"] = email
        summary = fetch_json(ESUMMARY_URL, sparams)
    payload = {"query": {"symbol": symbol, "taxid": taxid, "term": term}, "search": search, "summary": summary}
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    if sleep_s:
        time.sleep(sleep_s)
    return payload


def best_candidate(payload: dict, symbol: str) -> dict:
    ids = payload.get("search", {}).get("esearchresult", {}).get("idlist", [])
    result = payload.get("summary", {}).get("result", {})
    candidates = []
    for gene_id in ids:
        item = result.get(gene_id)
        if not isinstance(item, dict):
            continue
        name = str(item.get("name", ""))
        desc = str(item.get("description", ""))
        exact = name.upper() == symbol.upper()
        contains = symbol.upper() in name.upper() or symbol.upper() in desc.upper()
        candidates.append(
            {
                "gene_id": gene_id,
                "symbol": name,
                "description": desc,
                "exact_symbol": exact,
                "contains_symbol": contains,
            }
        )
    if not candidates:
        return {}
    candidates.sort(key=lambda x: (not x["exact_symbol"], not x["contains_symbol"], x["gene_id"]))
    return candidates[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--cache-dir", type=pathlib.Path, default=pathlib.Path("data/interim/phase2/ncbi_gene_pilot_cache"))
    parser.add_argument("--email", default="")
    parser.add_argument("--sleep", type=float, default=0.34)
    parser.add_argument("--max-rows", type=int, default=None)
    args = parser.parse_args()

    rows = pd.read_csv(args.input, sep="\t")
    if args.max_rows is not None:
        rows = rows.head(args.max_rows).copy()

    mapped = []
    counts: Counter[str] = Counter()
    for _, row in rows.iterrows():
        rec = row.to_dict()
        symbol = str(row["human_gene_symbol"])
        taxid = str(row["species_taxid"])
        try:
            payload = lookup(symbol, taxid, args.cache_dir, args.email, args.sleep)
            candidate = best_candidate(payload, symbol)
        except Exception as exc:  # noqa: BLE001
            rec.update(
                {
                    "ncbi_pilot_status": "query_error",
                    "ncbi_pilot_reason": str(exc),
                    "ncbi_gene_id": "",
                    "ncbi_gene_symbol": "",
                    "ncbi_gene_description": "",
                    "ncbi_gene_url": "",
                    "ncbi_symbol_confidence": "none",
                }
            )
            counts["query_error"] += 1
            mapped.append(rec)
            continue

        if not candidate:
            rec.update(
                {
                    "ncbi_pilot_status": "no_candidate",
                    "ncbi_pilot_reason": "no_ncbi_gene_symbol_taxid_hit",
                    "ncbi_gene_id": "",
                    "ncbi_gene_symbol": "",
                    "ncbi_gene_description": "",
                    "ncbi_gene_url": "",
                    "ncbi_symbol_confidence": "none",
                }
            )
            counts["no_candidate"] += 1
        else:
            confidence = "medium" if candidate["exact_symbol"] else "low"
            rec.update(
                {
                    "ncbi_pilot_status": "candidate_found",
                    "ncbi_pilot_reason": "exact_symbol" if candidate["exact_symbol"] else "non_exact_symbol",
                    "ncbi_gene_id": candidate["gene_id"],
                    "ncbi_gene_symbol": candidate["symbol"],
                    "ncbi_gene_description": candidate["description"],
                    "ncbi_gene_url": f"https://www.ncbi.nlm.nih.gov/gene/{candidate['gene_id']}",
                    "ncbi_symbol_confidence": confidence,
                }
            )
            counts["candidate_found"] += 1
        mapped.append(rec)

    out = pd.DataFrame(mapped)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    summary = (
        out.groupby(["human_gene_symbol", "maintenance_module_v2", "submodule_v2"], as_index=False)
        .agg(
            rows=("scientific_name", "size"),
            candidate_found=("ncbi_pilot_status", lambda x: int((x == "candidate_found").sum())),
            no_candidate=("ncbi_pilot_status", lambda x: int((x == "no_candidate").sum())),
            query_error=("ncbi_pilot_status", lambda x: int((x == "query_error").sum())),
            medium_confidence=("ncbi_symbol_confidence", lambda x: int((x == "medium").sum())),
        )
        .assign(candidate_fraction=lambda d: d["candidate_found"] / d["rows"])
        .sort_values(["candidate_fraction", "human_gene_symbol"])
    )
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary_output, sep="\t", index=False)

    lines = [
        "# Phase 2 P2.2 NCBI Gene Pilot Report",
        "",
        f"Input rows queried: {len(out)}",
        "",
        "## Status Counts",
        "",
    ]
    for status, count in sorted(counts.items()):
        lines.append(f"- {status}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "NCBI Gene symbol/taxid hits are first-pass candidates only. Genes with low candidate fractions should be prioritized for external orthology databases or domain-level sequence checks before inclusion in strict v2 scoring.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} with {len(out)} rows")


if __name__ == "__main__":
    main()
