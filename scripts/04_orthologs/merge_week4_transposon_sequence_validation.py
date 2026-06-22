"""Merge Week 4 transposon sequence validation back into the primary matrix."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DEFAULT_MATRIX = pathlib.Path("data/processed/ortholog_matrix_primary_diamond_validated_candidates.tsv")
DEFAULT_VALIDATION = pathlib.Path("data/processed/week4_transposon_reciprocal_validation.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/ortholog_matrix_primary_week4_sequence_validated.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/ortholog_matrix_primary_week4_sequence_validated_report.md")

SUPPORTED = {"reciprocal_supported", "reciprocal_weak"}


def sequence_status(row: pd.Series) -> str:
    if row["maintenance_module"] != "transposon_suppression":
        return "not_tested_non_transposon"
    if pd.isna(row.get("week4_diamond_status")):
        return "not_in_first_transposon_batch"
    if row["week4_diamond_status"] == "reciprocal_supported":
        return "sequence_supported"
    if row["week4_diamond_status"] == "reciprocal_weak":
        return "sequence_supported_weak"
    return "sequence_not_supported"


def week4_status(row: pd.Series) -> str:
    original = row["final_candidate_status"]
    seq = row["week4_sequence_status"]
    if seq == "sequence_supported":
        return "week4_sequence_supported_candidate"
    if seq == "sequence_supported_weak":
        return "week4_sequence_weak_candidate"
    if seq == "sequence_not_supported" and original in {"gff_rescue_candidate", "diamond_validated_protein_candidate"}:
        return "week4_sequence_not_supported"
    return original


def week4_confidence(row: pd.Series) -> str:
    seq = row["week4_sequence_status"]
    if seq == "sequence_supported":
        return "high"
    if seq == "sequence_supported_weak":
        return "medium"
    if row["week4_candidate_status"] == "week4_sequence_not_supported":
        return "low"
    return row.get("final_candidate_confidence", "")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, default=DEFAULT_MATRIX)
    parser.add_argument("--validation", type=pathlib.Path, default=DEFAULT_VALIDATION)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t")
    validation = pd.read_csv(args.validation, sep="\t")
    keep = [
        "scientific_name",
        "human_gene_symbol",
        "week4_diamond_status",
        "week4_diamond_reason",
        "forward_target_protein_id",
        "forward_pident",
        "forward_qcovhsp",
        "forward_scovhsp",
        "forward_evalue",
        "forward_bitscore",
        "reciprocal_top_reference",
        "reciprocal_gene",
        "reciprocal_pident",
        "reciprocal_qcovhsp",
        "reciprocal_scovhsp",
        "reciprocal_evalue",
        "reciprocal_bitscore",
    ]
    merged = matrix.merge(validation[keep], on=["scientific_name", "human_gene_symbol"], how="left")
    merged["week4_sequence_status"] = merged.apply(sequence_status, axis=1)
    merged["week4_candidate_status"] = merged.apply(week4_status, axis=1)
    merged["week4_candidate_confidence"] = merged.apply(week4_confidence, axis=1)
    merged["week4_candidate_source"] = merged["final_candidate_source"]
    supported_mask = merged["week4_sequence_status"].isin(["sequence_supported", "sequence_supported_weak"])
    merged.loc[supported_mask, "week4_candidate_source"] = "local reciprocal DIAMOND"

    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)

    trans = merged[merged["maintenance_module"] == "transposon_suppression"].copy()
    status_counts = trans["week4_sequence_status"].value_counts().sort_index()
    candidate_counts = trans["week4_candidate_status"].value_counts().sort_index()
    lines = [
        "# Ortholog Matrix Primary Week 4 Sequence-Validated Report",
        "",
        f"Matrix rows: {len(merged)}",
        f"Transposon rows: {len(trans)}",
        f"First-batch sequence-supported rows: {(trans['week4_sequence_status'] == 'sequence_supported').sum()}",
        f"First-batch weak sequence-supported rows: {(trans['week4_sequence_status'] == 'sequence_supported_weak').sum()}",
        "",
        "## Transposon Sequence Status",
        "",
    ]
    for status, count in status_counts.items():
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Week 4 Candidate Status", ""])
    for status, count in candidate_counts.items():
        lines.append(f"- {status}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This matrix keeps the original Week 3 candidate calls and adds Week 4 sequence-validation columns. Sequence-supported rows can be used for stricter transposon-suppression sensitivity scores.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
