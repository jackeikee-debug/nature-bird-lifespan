"""Run lightweight NCBI Protein searches for P2.2 external review rows."""

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


def cache_key(symbol: str, species: str) -> str:
    raw = f"{symbol}__{species}"
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in raw) + ".json"


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
    raise RuntimeError(f"NCBI protein request failed: {last_error}")


def query_protein(symbol: str, species: str, cache_dir: pathlib.Path, sleep_s: float) -> dict:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / cache_key(symbol, species)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    term = f'("{symbol}"[Protein Name] OR "{symbol}"[All Fields]) AND "{species}"[Organism]'
    params = {
        "db": "protein",
        "term": term,
        "retmode": "json",
        "retmax": "10",
        "tool": "bird_lifespan_phase2",
    }
    search = fetch_json(ESEARCH_URL, params)
    ids = search.get("esearchresult", {}).get("idlist", [])
    summary = {}
    if ids:
        summary = fetch_json(
            ESUMMARY_URL,
            {
                "db": "protein",
                "id": ",".join(ids),
                "retmode": "json",
                "tool": "bird_lifespan_phase2",
            },
        )
    payload = {"query": {"symbol": symbol, "species": species, "term": term}, "search": search, "summary": summary}
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    if sleep_s:
        time.sleep(sleep_s)
    return payload


def summarize(payload: dict, symbol: str) -> tuple[str, str, str, str, str]:
    ids = payload.get("search", {}).get("esearchresult", {}).get("idlist", [])
    result = payload.get("summary", {}).get("result", {})
    if not ids:
        return "no_protein_hit", "", "", "", ""
    titles = []
    accessions = []
    exactish = False
    for pid in ids[:5]:
        item = result.get(pid, {})
        title = str(item.get("title", ""))
        accession = str(item.get("accessionversion", ""))
        titles.append(title)
        accessions.append(accession)
        if symbol.upper() in title.upper():
            exactish = True
    status = "protein_symbol_like_hit" if exactish else "protein_broad_hit"
    return status, ";".join(accessions), ";".join(titles), str(len(ids)), ids[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--gene-summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--cache-dir", type=pathlib.Path, default=pathlib.Path("data/interim/phase2/ncbi_protein_external_review_cache"))
    parser.add_argument("--sleep", type=float, default=0.34)
    args = parser.parse_args()

    rows = pd.read_csv(args.input, sep="\t")
    out_rows = []
    counts: Counter[str] = Counter()
    for _, row in rows.iterrows():
        rec = row.to_dict()
        symbol = str(row["human_gene_symbol"])
        species = str(row["scientific_name"])
        try:
            payload = query_protein(symbol, species, args.cache_dir, args.sleep)
            status, accessions, titles, hit_count, first_id = summarize(payload, symbol)
        except Exception as exc:  # noqa: BLE001
            status, accessions, titles, hit_count, first_id = "query_error", "", str(exc), "0", ""
        rec.update(
            {
                "ncbi_protein_review_status": status,
                "ncbi_protein_hit_count": hit_count,
                "ncbi_protein_accessions_top": accessions,
                "ncbi_protein_titles_top": titles,
                "ncbi_protein_first_id": first_id,
            }
        )
        counts[status] += 1
        out_rows.append(rec)

    out = pd.DataFrame(out_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    summary = (
        out.groupby("human_gene_symbol", as_index=False)
        .agg(
            rows=("scientific_name", "size"),
            symbol_like_hits=("ncbi_protein_review_status", lambda x: int((x == "protein_symbol_like_hit").sum())),
            broad_hits=("ncbi_protein_review_status", lambda x: int((x == "protein_broad_hit").sum())),
            no_hits=("ncbi_protein_review_status", lambda x: int((x == "no_protein_hit").sum())),
            query_errors=("ncbi_protein_review_status", lambda x: int((x == "query_error").sum())),
        )
        .assign(symbol_like_fraction=lambda d: d["symbol_like_hits"] / d["rows"])
        .sort_values(["symbol_like_fraction", "human_gene_symbol"])
    )
    args.gene_summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.gene_summary_output, sep="\t", index=False)

    lines = [
        "# Phase 2 P2.2 NCBI Protein External Review Report",
        "",
        f"Rows queried: {len(out)}",
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
            "Protein search hits are supportive clues only. They do not establish orthology. Zero protein hits strengthen the case for lineage-specific absence, naming divergence, or annotation gaps, but still require external/database or domain-aware validation before strict scoring.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
