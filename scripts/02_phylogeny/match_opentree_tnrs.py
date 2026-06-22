"""Match prepared species names to OpenTree TNRS.

This creates a reproducible bridge between local species labels and OpenTree
taxonomy identifiers. Network access is required.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import urllib.error
import urllib.request


DEFAULT_AUDIT = pathlib.Path("data/processed/tree_label_audit.tsv")
DEFAULT_MATCHES = pathlib.Path("data/processed/opentree_tnrs_matches.tsv")
DEFAULT_UPDATED_AUDIT = pathlib.Path("data/processed/tree_label_audit_opentree.tsv")
DEFAULT_ENDPOINT = "https://api.opentreeoflife.org/v3/tnrs/match_names"


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def post_tnrs(endpoint: str, names: list[str], timeout: int) -> dict:
    payload = {
        "names": names,
        "do_approximate_matching": True,
        "include_suppressed": False,
    }
    request = urllib.request.Request(
        endpoint,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "User-Agent": "bird-lifespan-pgls-prep/0.1",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"OpenTree TNRS HTTP {exc.code}: {body}") from exc


def best_match(result: dict) -> dict | None:
    matches = result.get("matches") or []
    if not matches:
        return None
    return sorted(
        matches,
        key=lambda match: (
            0 if match.get("is_synonym") is False else 1,
            0 if match.get("is_approximate_match") is False else 1,
            -float(match.get("score") or 0),
        ),
    )[0]


def flatten_match(local_by_search: dict[str, dict[str, str]], result: dict) -> dict[str, str]:
    search_string = result.get("name") or result.get("id") or result.get("search_string") or ""
    local = local_by_search.get(search_string, {})
    match = best_match(result)
    if not match:
        return {
            "scientific_name": local.get("scientific_name", ""),
            "tree_search_name": search_string,
            "match_status": "unmatched",
            "ott_id": "",
            "matched_name": "",
            "unique_name": "",
            "rank": "",
            "score": "",
            "is_synonym": "",
            "is_approximate_match": "",
            "tax_sources": "",
        }
    taxon = match.get("taxon", {})
    return {
        "scientific_name": local.get("scientific_name", ""),
        "tree_search_name": search_string,
        "match_status": "matched",
        "ott_id": str(taxon.get("ott_id") or match.get("ott_id") or ""),
        "matched_name": taxon.get("name", "") or match.get("matched_name", ""),
        "unique_name": taxon.get("unique_name", "") or match.get("unique_name", ""),
        "rank": taxon.get("rank", ""),
        "score": str(match.get("score", "")),
        "is_synonym": str(match.get("is_synonym", "")),
        "is_approximate_match": str(match.get("is_approximate_match", "")),
        "tax_sources": ";".join(taxon.get("tax_sources", []) or []),
    }


def match_names(audit_rows: list[dict[str, str]], endpoint: str, chunk_size: int, timeout: int) -> list[dict[str, str]]:
    names = [row["tree_search_name"] for row in audit_rows]
    local_by_search = {row["tree_search_name"]: row for row in audit_rows}
    flattened = []
    for start in range(0, len(names), chunk_size):
        chunk = names[start : start + chunk_size]
        response = post_tnrs(endpoint, chunk, timeout)
        results = response.get("results", [])
        returned = set()
        for result in results:
            row = flatten_match(local_by_search, result)
            returned.add(row["tree_search_name"])
            flattened.append(row)
        for missing in sorted(set(chunk) - returned):
            flattened.append(
                {
                    "scientific_name": local_by_search[missing].get("scientific_name", ""),
                    "tree_search_name": missing,
                    "match_status": "unmatched",
                    "ott_id": "",
                    "matched_name": "",
                    "unique_name": "",
                    "rank": "",
                    "score": "",
                    "is_synonym": "",
                    "is_approximate_match": "",
                    "tax_sources": "",
                }
            )
    return flattened


def update_audit(audit_rows: list[dict[str, str]], matches: list[dict[str, str]]) -> list[dict[str, str]]:
    by_search = {row["tree_search_name"]: row for row in matches}
    updated = []
    for row in audit_rows:
        match = by_search.get(row["tree_search_name"], {})
        enriched = dict(row)
        enriched["tree_match_status"] = match.get("match_status", "unmatched")
        enriched["ott_id"] = match.get("ott_id", "")
        enriched["matched_tree_tip"] = match.get("matched_name", "")
        if match.get("is_approximate_match") == "True":
            notes = enriched.get("notes", "")
            enriched["notes"] = (notes + "; " if notes else "") + "OpenTree approximate match"
        if match.get("is_synonym") == "True":
            notes = enriched.get("notes", "")
            enriched["notes"] = (notes + "; " if notes else "") + "OpenTree synonym match"
        updated.append(enriched)
    return updated


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=pathlib.Path, default=DEFAULT_AUDIT)
    parser.add_argument("--matches-output", type=pathlib.Path, default=DEFAULT_MATCHES)
    parser.add_argument("--updated-audit-output", type=pathlib.Path, default=DEFAULT_UPDATED_AUDIT)
    parser.add_argument("--endpoint", default=DEFAULT_ENDPOINT)
    parser.add_argument("--chunk-size", type=int, default=50)
    parser.add_argument("--timeout", type=int, default=60)
    args = parser.parse_args()

    audit_rows = read_tsv(args.audit)
    matches = match_names(audit_rows, args.endpoint, args.chunk_size, args.timeout)
    match_fields = [
        "scientific_name",
        "tree_search_name",
        "match_status",
        "ott_id",
        "matched_name",
        "unique_name",
        "rank",
        "score",
        "is_synonym",
        "is_approximate_match",
        "tax_sources",
    ]
    write_tsv(args.matches_output, matches, match_fields)
    updated = update_audit(audit_rows, matches)
    write_tsv(args.updated_audit_output, updated, list(updated[0].keys()))
    matched = sum(row["match_status"] == "matched" for row in matches)
    print(f"Wrote {args.matches_output} and {args.updated_audit_output}; matched {matched}/{len(matches)}")


if __name__ == "__main__":
    main()

