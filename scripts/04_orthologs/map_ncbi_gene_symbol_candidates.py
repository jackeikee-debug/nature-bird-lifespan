"""Fill an ortholog scaffold with NCBI Gene symbol/taxid candidates.

This is a feasibility mapper, not a final orthology caller. It searches NCBI
Gene for each human gene symbol within each species taxid and records candidates
as symbol-and-taxid matches. Final manuscript claims should cross-check these
calls against an orthology source such as Ensembl Compara, OMA, or OrthoDB.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import time
import urllib.parse
import urllib.request
from collections import Counter


DEFAULT_INPUT = pathlib.Path("data/processed/ortholog_matrix_primary.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/ortholog_matrix_primary_ncbi_gene_candidates.tsv")
DEFAULT_CACHE = pathlib.Path("data/raw/ncbi_gene_cache")
DEFAULT_REPORT = pathlib.Path("results/reports/ncbi_gene_candidate_mapping_report.md")

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def cache_key(gene_symbol: str, taxid: str) -> str:
    safe_symbol = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in gene_symbol)
    return f"taxid_{taxid}__symbol_{safe_symbol}.json"


def fetch_json(url: str, params: dict[str, str], retries: int = 4) -> dict:
    encoded = urllib.parse.urlencode(params)
    full_url = f"{url}?{encoded}"
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(full_url, timeout=30) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as exc:  # noqa: BLE001 - preserve retry context for report.
            last_error = exc
            time.sleep(2 + attempt * 3)
    raise RuntimeError(f"NCBI request failed after {retries} attempts: {last_error}")


def ncbi_gene_lookup(gene_symbol: str, taxid: str, cache_dir: pathlib.Path, email: str) -> dict:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / cache_key(gene_symbol, taxid)
    if cache_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))

    term = f'("{gene_symbol}"[Gene Name] OR "{gene_symbol}"[Preferred Symbol]) AND txid{taxid}[Organism]'
    params = {
        "db": "gene",
        "term": term,
        "retmode": "json",
        "retmax": "5",
        "tool": "bird_lifespan_feasibility",
    }
    if email:
        params["email"] = email
    search = fetch_json(ESEARCH_URL, params)
    ids = search.get("esearchresult", {}).get("idlist", [])
    summary = {}
    if ids:
        summary_params = {
            "db": "gene",
            "id": ",".join(ids),
            "retmode": "json",
            "tool": "bird_lifespan_feasibility",
        }
        if email:
            summary_params["email"] = email
        summary = fetch_json(ESUMMARY_URL, summary_params)

    payload = {"query": {"gene_symbol": gene_symbol, "taxid": taxid, "term": term}, "search": search, "summary": summary}
    cache_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    time.sleep(0.35)
    return payload


def best_candidate(payload: dict, requested_symbol: str) -> dict[str, str]:
    ids = payload.get("search", {}).get("esearchresult", {}).get("idlist", [])
    if not ids:
        return {}
    result = payload.get("summary", {}).get("result", {})
    candidates = []
    for gene_id in ids:
        item = result.get(gene_id)
        if not isinstance(item, dict):
            continue
        name = str(item.get("name", ""))
        description = str(item.get("description", ""))
        exact_symbol = name.upper() == requested_symbol.upper()
        candidates.append(
            {
                "gene_id": gene_id,
                "symbol": name,
                "description": description,
                "exact_symbol": exact_symbol,
            }
        )
    if not candidates:
        return {}
    candidates.sort(key=lambda item: (not item["exact_symbol"], item["gene_id"]))
    return candidates[0]


def map_rows(
    rows: list[dict[str, str]],
    cache_dir: pathlib.Path,
    email: str,
    max_queries: int | None,
) -> tuple[list[dict[str, str]], Counter[str]]:
    mapped = []
    status_counts: Counter[str] = Counter()
    query_count = 0
    for row in rows:
        updated = dict(row)
        if max_queries is not None and query_count >= max_queries:
            updated["ortholog_query_status"] = "not_queried_max_queries_reached"
            status_counts[updated["ortholog_query_status"]] += 1
            mapped.append(updated)
            continue
        gene_symbol = row["human_gene_symbol"]
        taxid = row["species_taxid"]
        query_count += 1
        try:
            payload = ncbi_gene_lookup(gene_symbol, taxid, cache_dir, email)
            candidate = best_candidate(payload, gene_symbol)
        except Exception as exc:  # noqa: BLE001 - keep row-level failure visible.
            updated["ortholog_query_status"] = "query_error"
            updated["notes"] = str(exc)
            status_counts[updated["ortholog_query_status"]] += 1
            mapped.append(updated)
            continue
        if not candidate:
            updated["ortholog_query_status"] = "no_ncbi_gene_candidate"
            updated["ortholog_status"] = "not_found"
            updated["ortholog_source_database"] = "NCBI Gene"
        else:
            updated["ortholog_query_status"] = "candidate_found"
            updated["ortholog_gene_id"] = candidate["gene_id"]
            updated["ortholog_gene_symbol"] = candidate["symbol"]
            updated["ortholog_status"] = "candidate_by_symbol_taxid"
            updated["ortholog_source_database"] = "NCBI Gene"
            updated["ortholog_source_url"] = f"https://www.ncbi.nlm.nih.gov/gene/{candidate['gene_id']}"
            updated["ortholog_confidence"] = "medium" if candidate["exact_symbol"] else "low"
            updated["copy_number_estimate"] = "not_estimated"
            updated["notes"] = candidate["description"]
        status_counts[updated["ortholog_query_status"]] += 1
        mapped.append(updated)
    return mapped, status_counts


def write_report(path: pathlib.Path, rows: list[dict[str, str]], status_counts: Counter[str]) -> None:
    clade_found = Counter(row["clade"] for row in rows if row["ortholog_query_status"] == "candidate_found")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# NCBI Gene Candidate Mapping Report",
                "",
                f"Input rows: {len(rows)}",
                "",
                "## Query Status",
                *[f"- {status}: {count}" for status, count in sorted(status_counts.items())],
                "",
                "## Candidate Found by Clade",
                *([f"- {clade}: {count}" for clade, count in sorted(clade_found.items())] or ["- none: 0"]),
                "",
                "## Interpretation",
                "These are NCBI Gene symbol/taxid candidates, not final ortholog assertions. Use them as a practical first-pass fill for feasibility and cross-check with Ensembl Compara, OMA, or OrthoDB before pathway-level claims.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cache-dir", type=pathlib.Path, default=DEFAULT_CACHE)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    parser.add_argument("--email", default="")
    parser.add_argument("--max-queries", type=int, default=None)
    args = parser.parse_args()

    rows = read_tsv(args.input)
    mapped, status_counts = map_rows(rows, args.cache_dir, args.email, args.max_queries)
    fields = list(rows[0].keys())
    write_tsv(args.output, mapped, fields)
    write_report(args.report, mapped, status_counts)
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
