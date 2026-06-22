"""Build strict decisions for assembly CDS translation rescue candidates."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def join_unique(values: pd.Series) -> str:
    vals = [str(v) for v in values if str(v) and str(v) != "nan"]
    return ";".join(sorted(set(vals)))


def best_decision(calls: set[str]) -> tuple[str, str, bool]:
    if "reciprocal_same_gene_supported" in calls:
        return (
            "cds_translation_strict_rescue",
            "Assembly CDS translation passed reciprocal same-gene DIAMOND/BLASTP validation.",
            True,
        )
    if "same_gene_forward_supported_reciprocal_weak" in calls:
        return (
            "cds_translation_forward_supported_manual_review",
            "Assembly CDS translation has same-gene forward support but reciprocal evidence is incomplete.",
            False,
        )
    if "weak_same_gene_support" in calls:
        return (
            "cds_translation_weak_same_gene_not_strict",
            "Assembly CDS translation has same-gene support below strict reciprocal thresholds.",
            False,
        )
    if "no_same_gene_reference_hit" in calls:
        return (
            "reject_cds_translation_no_same_gene_reference",
            "Assembly CDS translation best hits do not include a same-gene human reference.",
            False,
        )
    return (
        "cds_translation_not_validated",
        "No reciprocal validation record was available for this assembly CDS translation.",
        False,
    )


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
        decisions.append(best_decision(calls))
    merged["phase3_cds_translation_decision"] = [x[0] for x in decisions]
    merged["phase3_cds_translation_reason"] = [x[1] for x in decisions]
    merged["can_count_as_strict_sequence_after_cds_translation"] = [str(x[2]) for x in decisions]

    keep_cols = [
        "scientific_name",
        "clade",
        "flight_status",
        "species_taxid",
        "best_assembly_accession",
        "genome_analysis_tier",
        "maintenance_module",
        "human_gene_symbol",
        "gene_family_risk",
        "gff_rescue_call",
        "matched_id",
        "matched_gene",
        "cds_local_path",
        "cds_locus_tag",
        "cds_protein_id",
        "cds_product",
        "cds_partial",
        "fetched_accession",
        "sequence_fetch_status",
        "protein_length",
        "validation_tools",
        "validation_calls",
        "top_reference_genes",
        "reciprocal_best_genes",
        "max_pident",
        "max_qcovhsp",
        "max_scovhsp",
        "phase3_cds_translation_decision",
        "phase3_cds_translation_reason",
        "can_count_as_strict_sequence_after_cds_translation",
    ]
    for col in keep_cols:
        if col not in merged.columns:
            merged[col] = ""
    output = merged[keep_cols].copy()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    summary = (
        output.groupby(["human_gene_symbol", "phase3_cds_translation_decision"], as_index=False)
        .agg(rows=("scientific_name", "count"))
        .sort_values(["human_gene_symbol", "phase3_cds_translation_decision"])
    )
    summary.to_csv(args.summary_output, sep="\t", index=False)

    counts = output["phase3_cds_translation_decision"].value_counts().sort_index()
    strict = output[output["phase3_cds_translation_decision"] == "cds_translation_strict_rescue"]
    lines = [
        "# Phase 3 CDS Translation Rescue Decision Report",
        "",
        f"CDS translation candidate rows assessed: {len(output)}",
        f"Rows upgraded to strict local CDS translation rescue: {len(strict)}",
        "",
        "## Decision Counts",
    ]
    for decision, count in counts.items():
        lines.append(f"- {decision}: {count}")
    lines.extend(["", "## Strict Rescue Rows"])
    for _, row in strict.sort_values(["scientific_name", "human_gene_symbol"]).iterrows():
        lines.append(f"- {row['scientific_name']} / {row['human_gene_symbol']}: {row['fetched_accession']}")
    if strict.empty:
        lines.append("- none: 0")
    lines.extend(
        [
            "",
            "## Interpretation",
            "Strict rows here are local assembly CDS translations, so they are stronger than external-only UniProt rescue and can be treated as assembly-supported sequence rescue, while still flagged separately from GFF protein_id rescue.",
            "",
            "## Outputs",
            f"- decisions: `{args.output}`",
            f"- summary: `{args.summary_output}`",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
