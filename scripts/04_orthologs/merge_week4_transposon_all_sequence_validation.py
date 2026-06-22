"""Merge rescue and NCBI-direct transposon sequence validation into one matrix."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DEFAULT_MATRIX = pathlib.Path("data/processed/ortholog_matrix_primary_diamond_validated_candidates.tsv")
DEFAULT_RESCUE_VALIDATION = pathlib.Path("data/processed/week4_transposon_reciprocal_validation.tsv")
DEFAULT_NCBI_VALIDATION = pathlib.Path("data/processed/week4_transposon_ncbi_crosscheck_validation.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/ortholog_matrix_primary_week4_full_sequence_validated.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/ortholog_matrix_primary_week4_full_sequence_validated_report.md")


VALIDATION_KEEP = [
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


def source_status(row: pd.Series) -> str:
    if row["maintenance_module"] != "transposon_suppression":
        return "not_tested_non_transposon"
    if pd.notna(row.get("rescue_week4_diamond_status")):
        return "rescue_or_unresolved_batch"
    if pd.notna(row.get("ncbi_week4_diamond_status")):
        return "ncbi_direct_batch"
    return "not_tested_transposon"


def sequence_status(row: pd.Series) -> str:
    if row["maintenance_module"] != "transposon_suppression":
        return "not_tested_non_transposon"
    status = row.get("rescue_week4_diamond_status")
    if pd.isna(status):
        status = row.get("ncbi_week4_diamond_status")
    if pd.isna(status):
        return "not_tested_transposon"
    if status == "reciprocal_supported":
        return "sequence_supported"
    if status == "reciprocal_weak":
        return "sequence_supported_weak"
    return "sequence_not_supported"


def candidate_status(row: pd.Series) -> str:
    original = row["final_candidate_status"]
    seq = row["week4_sequence_status"]
    if seq == "sequence_supported":
        return "week4_sequence_supported_candidate"
    if seq == "sequence_supported_weak":
        return "week4_sequence_weak_candidate"
    if seq == "sequence_not_supported":
        return "week4_sequence_not_supported"
    return original


def confidence(row: pd.Series) -> str:
    status = row["week4_candidate_status"]
    if status == "week4_sequence_supported_candidate":
        return "high"
    if status == "week4_sequence_weak_candidate":
        return "medium"
    if status == "week4_sequence_not_supported":
        return "low"
    return row.get("final_candidate_confidence", "")


def source(row: pd.Series) -> str:
    if row["week4_candidate_status"].startswith("week4_sequence"):
        if row["week4_validation_batch_source"] == "ncbi_direct_batch":
            return "NCBI Gene plus local reciprocal DIAMOND"
        return "local reciprocal DIAMOND"
    return row.get("final_candidate_source", "")


def prefix_validation(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    keep = [col for col in VALIDATION_KEEP if col in df.columns]
    out = df[keep].copy()
    rename = {
        col: f"{prefix}_{col}"
        for col in keep
        if col not in {"scientific_name", "human_gene_symbol"}
    }
    return out.rename(columns=rename)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, default=DEFAULT_MATRIX)
    parser.add_argument("--rescue-validation", type=pathlib.Path, default=DEFAULT_RESCUE_VALIDATION)
    parser.add_argument("--ncbi-validation", type=pathlib.Path, default=DEFAULT_NCBI_VALIDATION)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t")
    rescue = prefix_validation(pd.read_csv(args.rescue_validation, sep="\t"), "rescue")
    ncbi = prefix_validation(pd.read_csv(args.ncbi_validation, sep="\t"), "ncbi")
    merged = matrix.merge(rescue, on=["scientific_name", "human_gene_symbol"], how="left")
    merged = merged.merge(ncbi, on=["scientific_name", "human_gene_symbol"], how="left")
    merged["week4_validation_batch_source"] = merged.apply(source_status, axis=1)
    merged["week4_sequence_status"] = merged.apply(sequence_status, axis=1)
    merged["week4_candidate_status"] = merged.apply(candidate_status, axis=1)
    merged["week4_candidate_confidence"] = merged.apply(confidence, axis=1)
    merged["week4_candidate_source"] = merged.apply(source, axis=1)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)

    trans = merged[merged["maintenance_module"] == "transposon_suppression"].copy()
    batch_counts = trans["week4_validation_batch_source"].value_counts().sort_index()
    seq_counts = trans["week4_sequence_status"].value_counts().sort_index()
    cand_counts = trans["week4_candidate_status"].value_counts().sort_index()
    gene_seq = trans.groupby(["human_gene_symbol", "week4_sequence_status"]).size().reset_index(name="rows")

    lines = [
        "# Ortholog Matrix Primary Week 4 Full Sequence-Validated Report",
        "",
        f"Matrix rows: {len(merged)}",
        f"Transposon rows: {len(trans)}",
        "",
        "## Validation Batch Source",
        "",
    ]
    for status, count in batch_counts.items():
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Sequence Status", ""])
    for status, count in seq_counts.items():
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Candidate Status", ""])
    for status, count in cand_counts.items():
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Gene x Sequence Status", ""])
    for _, rec in gene_seq.iterrows():
        lines.append(f"- {rec['human_gene_symbol']} / {rec['week4_sequence_status']}: {rec['rows']}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This matrix requires local reciprocal DIAMOND support for both rescue/unresolved transposon rows and direct NCBI Gene transposon rows. It is the strictest Week 4 transposon orthology matrix so far.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
