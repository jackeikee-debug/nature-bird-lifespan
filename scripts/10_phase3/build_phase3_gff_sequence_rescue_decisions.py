"""Build strict rescue decisions from Phase 3 GFF-linked protein validation."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def join_unique(values: pd.Series) -> str:
    vals = [str(v) for v in values if str(v) and str(v) != "nan"]
    return ";".join(sorted(set(vals)))


def best_decision(calls: set[str]) -> tuple[str, str, bool, bool]:
    if "reciprocal_same_gene_supported" in calls:
        return (
            "gff_sequence_strict_rescue",
            "At least one reciprocal DIAMOND/BLASTP validation returned the same human gene as the top reciprocal hit.",
            True,
            True,
        )
    if "same_gene_forward_supported_reciprocal_weak" in calls:
        return (
            "gff_sequence_forward_supported_manual_review",
            "Same-gene forward hit passed relaxed thresholds but reciprocal evidence is incomplete.",
            True,
            False,
        )
    if "weak_same_gene_support" in calls:
        return (
            "gff_sequence_weak_same_gene_not_strict",
            "A same-gene hit exists but does not meet reciprocal/coverage thresholds for strict scoring.",
            True,
            False,
        )
    if "no_same_gene_reference_hit" in calls:
        return (
            "reject_gff_sequence_no_same_gene_reference",
            "Candidate protein best hits do not include a same-gene human reference; likely wrong gene/paralog or incomplete reference.",
            False,
            False,
        )
    return (
        "gff_sequence_not_validated",
        "No DIAMOND/BLASTP validation record was available for this fetched GFF protein.",
        True,
        False,
    )


def decision_rank(decision: str) -> int:
    return {
        "gff_sequence_strict_rescue": 4,
        "gff_sequence_forward_supported_manual_review": 3,
        "gff_sequence_weak_same_gene_not_strict": 2,
        "gff_sequence_not_validated": 1,
        "reject_gff_sequence_no_same_gene_reference": 0,
    }.get(decision, 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=pathlib.Path, required=True)
    parser.add_argument("--validation", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")
    validation = pd.read_csv(args.validation, sep="\t", dtype=str).fillna("")
    validation["pident_num"] = pd.to_numeric(validation.get("pident", ""), errors="coerce")
    validation["qcovhsp_num"] = pd.to_numeric(validation.get("qcovhsp", ""), errors="coerce")
    validation["scovhsp_num"] = pd.to_numeric(validation.get("scovhsp", ""), errors="coerce")

    agg = (
        validation.groupby(["scientific_name", "human_gene_symbol", "candidate_protein_accession"], as_index=False)
        .agg(
            validation_tools=("tool", join_unique),
            validation_calls=("validation_call", join_unique),
            top_reference_genes=("top_reference_gene", join_unique),
            reciprocal_best_genes=("reciprocal_best_gene", join_unique),
            max_pident=("pident_num", "max"),
            max_qcovhsp=("qcovhsp_num", "max"),
            max_scovhsp=("scovhsp_num", "max"),
        )
    )
    merged = metadata.merge(
        agg,
        left_on=["scientific_name", "human_gene_symbol", "fetched_accession"],
        right_on=["scientific_name", "human_gene_symbol", "candidate_protein_accession"],
        how="left",
    )
    for col in ["validation_tools", "validation_calls", "top_reference_genes", "reciprocal_best_genes"]:
        merged[col] = merged[col].fillna("")
    for col in ["max_pident", "max_qcovhsp", "max_scovhsp"]:
        merged[col] = merged[col].fillna("")

    decisions = []
    for _, row in merged.iterrows():
        calls = {x for x in str(row["validation_calls"]).split(";") if x}
        decision, reason, coverage, strict = best_decision(calls)
        decisions.append((decision, reason, coverage, strict))
    merged["phase3_gff_sequence_decision"] = [x[0] for x in decisions]
    merged["phase3_gff_sequence_reason"] = [x[1] for x in decisions]
    merged["can_count_as_rescued_for_coverage_after_gff_sequence"] = [str(x[2]) for x in decisions]
    merged["can_count_as_strict_sequence_after_gff_sequence"] = [str(x[3]) for x in decisions]

    keep_cols = [
        "phase3_batch_id",
        "scientific_name",
        "clade",
        "flight_status",
        "species_taxid",
        "best_assembly_accession",
        "genome_analysis_tier",
        "maintenance_module",
        "human_gene_symbol",
        "gene_family_risk",
        "previous_phase3_rescue_decision",
        "gff_rescue_call",
        "gff_rescue_interpretation",
        "matched_id",
        "matched_gene",
        "protein_id",
        "gff_mrna_id",
        "gff_cds_id",
        "cds_product",
        "partial",
        "sequence_fetch_status",
        "protein_length",
        "validation_tools",
        "validation_calls",
        "top_reference_genes",
        "reciprocal_best_genes",
        "max_pident",
        "max_qcovhsp",
        "max_scovhsp",
        "phase3_gff_sequence_decision",
        "phase3_gff_sequence_reason",
        "can_count_as_rescued_for_coverage_after_gff_sequence",
        "can_count_as_strict_sequence_after_gff_sequence",
    ]
    for col in keep_cols:
        if col not in merged.columns:
            merged[col] = ""
    output = merged[keep_cols].copy()
    output["_decision_rank"] = output["phase3_gff_sequence_decision"].map(decision_rank).fillna(0)
    output["_max_pident_num"] = pd.to_numeric(output["max_pident"], errors="coerce").fillna(-1)
    output["_protein_length_num"] = pd.to_numeric(output["protein_length"], errors="coerce").fillna(-1)
    output = (
        output.sort_values(
            ["scientific_name", "human_gene_symbol", "_decision_rank", "_max_pident_num", "_protein_length_num"],
            ascending=[True, True, False, False, False],
        )
        .groupby(["scientific_name", "human_gene_symbol"], as_index=False)
        .first()
        .drop(columns=["_decision_rank", "_max_pident_num", "_protein_length_num"])
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)

    summary = (
        output.groupby(["scientific_name", "phase3_gff_sequence_decision"], as_index=False)
        .agg(rows=("human_gene_symbol", "count"))
        .sort_values(["scientific_name", "phase3_gff_sequence_decision"])
    )
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary_output, sep="\t", index=False)

    decision_counts = output["phase3_gff_sequence_decision"].value_counts().sort_index()
    strict_rows = int((output["can_count_as_strict_sequence_after_gff_sequence"] == "True").sum())
    coverage_rows = int((output["can_count_as_rescued_for_coverage_after_gff_sequence"] == "True").sum())
    species_strict = (
        output[output["can_count_as_strict_sequence_after_gff_sequence"] == "True"]
        .groupby("scientific_name")["human_gene_symbol"]
        .nunique()
        .sort_index()
    )
    lines = [
        "# Phase 3 GFF Sequence Rescue Decision Report",
        "",
        f"GFF-linked protein rows assessed: {len(output)}",
        f"Rows retaining annotation/coverage rescue after sequence screen: {coverage_rows}",
        f"Rows upgraded to strict sequence rescue: {strict_rows}",
        "",
        "## Decision Counts",
    ]
    for decision, count in decision_counts.items():
        lines.append(f"- {decision}: {count}")
    lines.extend(["", "## Strict Rescue by Species"])
    for species, count in species_strict.items():
        lines.append(f"- {species}: {count}")
    if species_strict.empty:
        lines.append("- none: 0")
    lines.extend(
        [
            "",
            "## Interpretation",
            "Rows upgraded here have both assembly/GFF annotation support and reciprocal same-gene sequence support. Weak or no-same-gene rows remain useful for auditing annotation bias, but should not enter the strict score.",
            "",
            "## Outputs",
            f"- decisions: `{args.output}`",
            f"- summary: `{args.summary_output}`",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
