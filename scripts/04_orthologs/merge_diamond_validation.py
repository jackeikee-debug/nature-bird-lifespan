"""Merge DIAMOND validation into the primary final candidate matrix."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import Counter


DEFAULT_MATRIX = pathlib.Path("data/processed/ortholog_matrix_primary_final_candidates.tsv")
DEFAULT_VALIDATION = pathlib.Path("data/processed/erythrura_diamond_validation.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/ortholog_matrix_primary_diamond_validated_candidates.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/ortholog_matrix_primary_diamond_validated_candidates_report.md")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def key(row: dict[str, str]) -> tuple[str, str]:
    return row["scientific_name"], row["human_gene_symbol"]


def merge_rows(matrix_rows: list[dict[str, str]], validation_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    validation = {key(row): row for row in validation_rows}
    merged = []
    for row in matrix_rows:
        updated = dict(row)
        val = validation.get(key(row), {})
        updated["diamond_validation_status"] = val.get("diamond_validation_status", "")
        updated["diamond_validation_reason"] = val.get("diamond_validation_reason", "")
        updated["diamond_pident"] = val.get("diamond_pident", "")
        updated["diamond_query_coverage"] = val.get("diamond_query_coverage", "")
        updated["diamond_subject_coverage"] = val.get("diamond_subject_coverage", "")
        updated["diamond_evalue"] = val.get("diamond_evalue", "")
        updated["diamond_bitscore"] = val.get("diamond_bitscore", "")
        updated["diamond_expected_target_rank"] = val.get("expected_target_rank", "")

        if row["final_candidate_status"] == "protein_similarity_candidate":
            if val.get("diamond_validation_status", "").startswith("diamond_validated"):
                confidence = "high" if val["diamond_validation_status"].endswith("high") else "medium"
                updated["final_candidate_status"] = "diamond_validated_protein_candidate"
                updated["final_candidate_source"] = "DIAMOND Taeniopygia reference protein vs Erythrura protein FASTA"
                updated["final_candidate_confidence"] = confidence
                updated["ortholog_confidence"] = confidence
                updated["notes"] = (
                    f"{row['notes']}; diamond_status={val['diamond_validation_status']}; "
                    f"pident={val['diamond_pident']}; qcov={val['diamond_query_coverage']}; "
                    f"evalue={val['diamond_evalue']}; bitscore={val['diamond_bitscore']}"
                )
            else:
                updated["final_candidate_status"] = "not_found_after_diamond_validation"
                updated["final_candidate_source"] = ""
                updated["final_candidate_confidence"] = ""
                updated["ortholog_status"] = "not_validated"
                updated["ortholog_confidence"] = ""
                updated["notes"] = (
                    f"{row['notes']}; diamond_status={val.get('diamond_validation_status', 'missing')}; "
                    f"reason={val.get('diamond_validation_reason', 'missing_validation')}"
                )
        merged.append(updated)
    return merged


def write_report(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    status_counts = Counter(row["final_candidate_status"] for row in rows)
    found_statuses = {
        "ncbi_gene_candidate",
        "gff_rescue_candidate",
        "diamond_validated_protein_candidate",
    }
    total_found = sum(count for status, count in status_counts.items() if status in found_statuses)
    coverage = total_found / len(rows) if rows else 0.0
    ery = [row for row in rows if row["scientific_name"] == "Erythrura gouldiae"]
    ery_found = [row for row in ery if row["final_candidate_status"] in found_statuses]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Primary DIAMOND-Validated Candidate Matrix Report",
                "",
                f"Input rows: {len(rows)}",
                f"Validated candidate hits: {total_found}",
                f"Validated candidate coverage: {coverage:.3f}",
                f"Erythrura gouldiae validated candidate hits: {len(ery_found)}/{len(ery)}",
                "",
                "## Final Candidate Status",
                *[f"- {status}: {count}" for status, count in sorted(status_counts.items())],
                "",
                "## Interpretation",
                "Protein-similarity rescue rows have been filtered through DIAMOND validation. Rows that fail DIAMOND thresholds are not counted as found in the validated feasibility matrix.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, default=DEFAULT_MATRIX)
    parser.add_argument("--validation", type=pathlib.Path, default=DEFAULT_VALIDATION)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    matrix_rows = read_tsv(args.matrix)
    validation_rows = read_tsv(args.validation)
    merged = merge_rows(matrix_rows, validation_rows)
    fields = list(matrix_rows[0].keys()) + [
        "diamond_validation_status",
        "diamond_validation_reason",
        "diamond_pident",
        "diamond_query_coverage",
        "diamond_subject_coverage",
        "diamond_evalue",
        "diamond_bitscore",
        "diamond_expected_target_rank",
    ]
    write_tsv(args.output, merged, fields)
    write_report(args.report, merged)
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
