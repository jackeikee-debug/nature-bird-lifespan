"""Merge GFF annotation rescue hits into the primary candidate ortholog matrix."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import Counter


DEFAULT_MATRIX = pathlib.Path("data/processed/ortholog_matrix_primary_ncbi_gene_candidates.tsv")
DEFAULT_RESCUE = pathlib.Path("data/processed/annotation_rescue_gene_hits.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/ortholog_matrix_primary_combined_candidates.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/ortholog_matrix_primary_combined_candidates_report.md")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def rescue_key(row: dict[str, str]) -> tuple[str, str]:
    return row["scientific_name"], row["human_gene_symbol"]


def merge_rows(matrix_rows: list[dict[str, str]], rescue_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rescue_hits = {
        rescue_key(row): row
        for row in rescue_rows
        if row["rescue_status"] == "gff_symbol_hit"
    }
    merged = []
    for row in matrix_rows:
        updated = dict(row)
        hit = rescue_hits.get(rescue_key(row))
        if row["ortholog_query_status"] == "candidate_found":
            updated["combined_candidate_status"] = "ncbi_gene_candidate"
            updated["combined_candidate_source"] = "NCBI Gene"
            updated["combined_candidate_confidence"] = row["ortholog_confidence"] or "medium"
        elif hit:
            updated["ortholog_query_status"] = "rescued_from_gff"
            updated["ortholog_gene_symbol"] = hit["feature_name"] or row["human_gene_symbol"]
            updated["ortholog_status"] = "candidate_by_gff_symbol"
            updated["ortholog_source_database"] = "NCBI GenBank GFF"
            updated["ortholog_source_url"] = ""
            updated["ortholog_confidence"] = "low"
            updated["copy_number_estimate"] = "not_estimated"
            updated["notes"] = hit["matched_attribute"]
            updated["combined_candidate_status"] = "gff_rescue_candidate"
            updated["combined_candidate_source"] = "NCBI GenBank GFF"
            updated["combined_candidate_confidence"] = "low"
        else:
            updated["combined_candidate_status"] = "not_found"
            updated["combined_candidate_source"] = ""
            updated["combined_candidate_confidence"] = ""
        merged.append(updated)
    return merged


def write_report(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    status_counts = Counter(row["combined_candidate_status"] for row in rows)
    clade_counts = Counter(
        row["clade"] for row in rows if row["combined_candidate_status"] != "not_found"
    )
    total_found = sum(count for status, count in status_counts.items() if status != "not_found")
    coverage = total_found / len(rows) if rows else 0.0
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Primary Combined Candidate Matrix Report",
                "",
                f"Input rows: {len(rows)}",
                f"Combined candidate hits: {total_found}",
                f"Combined candidate coverage: {coverage:.3f}",
                "",
                "## Candidate Status",
                *[f"- {status}: {count}" for status, count in sorted(status_counts.items())],
                "",
                "## Combined Hits by Clade",
                *[f"- {clade}: {count}" for clade, count in sorted(clade_counts.items())],
                "",
                "## Interpretation",
                "This matrix combines NCBI Gene symbol/taxid candidates with low-confidence GFF symbol rescue hits. GFF rescue prevents database missingness from being misread as gene loss, but these rows still need orthology validation before mechanism claims.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, default=DEFAULT_MATRIX)
    parser.add_argument("--rescue", type=pathlib.Path, default=DEFAULT_RESCUE)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    matrix_rows = read_tsv(args.matrix)
    rescue_rows = read_tsv(args.rescue)
    merged = merge_rows(matrix_rows, rescue_rows)
    fields = list(matrix_rows[0].keys()) + [
        "combined_candidate_status",
        "combined_candidate_source",
        "combined_candidate_confidence",
    ]
    write_tsv(args.output, merged, fields)
    write_report(args.report, merged)
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
