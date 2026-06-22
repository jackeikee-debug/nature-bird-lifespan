"""Build decisions for Phase 3 external UniProt sequence rescue."""

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
            "uniprot_full_length_strict_rescue",
            "UniProt target-species sequence passed reciprocal same-gene DIAMOND/BLASTP validation.",
            True,
        )
    if "same_gene_forward_supported_reciprocal_weak" in calls:
        return (
            "uniprot_full_length_forward_supported_manual_review",
            "UniProt target-species sequence has same-gene forward support but reciprocal evidence is incomplete.",
            False,
        )
    if "weak_same_gene_support" in calls:
        return (
            "uniprot_full_length_weak_same_gene_not_strict",
            "UniProt target-species sequence has same-gene support below strict reciprocal thresholds.",
            False,
        )
    if "no_same_gene_reference_hit" in calls:
        return (
            "reject_uniprot_full_length_no_same_gene_reference",
            "UniProt candidate best hits do not include a same-gene human reference.",
            False,
        )
    return (
        "uniprot_full_length_not_validated",
        "No reciprocal validation record was available for this UniProt candidate.",
        False,
    )


def decision_rank(decision: str) -> int:
    return {
        "uniprot_full_length_strict_rescue": 4,
        "uniprot_full_length_forward_supported_manual_review": 3,
        "uniprot_full_length_weak_same_gene_not_strict": 2,
        "uniprot_full_length_not_validated": 1,
        "reject_uniprot_full_length_no_same_gene_reference": 0,
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
    prior_validation_cols = {
        "validation_tools": "previous_gff_validation_tools",
        "validation_calls": "previous_gff_validation_calls",
        "top_reference_genes": "previous_gff_top_reference_genes",
        "reciprocal_best_genes": "previous_gff_reciprocal_best_genes",
        "max_pident": "previous_gff_max_pident",
        "max_qcovhsp": "previous_gff_max_qcovhsp",
        "max_scovhsp": "previous_gff_max_scovhsp",
    }
    metadata = metadata.rename(columns={old: new for old, new in prior_validation_cols.items() if old in metadata.columns})
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
    merged["phase3_uniprot_sequence_decision"] = [x[0] for x in decisions]
    merged["phase3_uniprot_sequence_reason"] = [x[1] for x in decisions]
    merged["can_count_as_strict_sequence_after_uniprot_rescue"] = [str(x[2]) for x in decisions]
    merged["phase3_uniprot_evidence_scope"] = "external_uniprot_sequence_sensitivity"

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
        "family_domain_review_class",
        "phase3_uniprot_rescue_scope",
        "phase3_uniprot_rescue_reason",
        "uniprot_crosscheck_call",
        "uniprot_accession_top",
        "fetched_accession",
        "fasta_header",
        "sequence_fetch_status",
        "protein_length",
        "validation_tools",
        "validation_calls",
        "top_reference_genes",
        "reciprocal_best_genes",
        "max_pident",
        "max_qcovhsp",
        "max_scovhsp",
        "phase3_uniprot_sequence_decision",
        "phase3_uniprot_sequence_reason",
        "can_count_as_strict_sequence_after_uniprot_rescue",
        "phase3_uniprot_evidence_scope",
    ]
    for col in keep_cols:
        if col not in merged.columns:
            merged[col] = ""
    output = merged[keep_cols].copy()
    output["_decision_rank"] = output["phase3_uniprot_sequence_decision"].map(decision_rank).fillna(0)
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
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)

    summary = (
        output.groupby(["family_domain_review_class", "human_gene_symbol", "phase3_uniprot_sequence_decision"], as_index=False)
        .agg(rows=("scientific_name", "count"), species=("scientific_name", "nunique"))
        .sort_values(["family_domain_review_class", "human_gene_symbol", "phase3_uniprot_sequence_decision"])
    )
    summary.to_csv(args.summary_output, sep="\t", index=False)

    decision_counts = output["phase3_uniprot_sequence_decision"].value_counts().sort_index()
    strict = output[output["phase3_uniprot_sequence_decision"] == "uniprot_full_length_strict_rescue"]
    lines = [
        "# Phase 3 Full-Length UniProt Sequence Rescue Decision Report",
        "",
        f"UniProt candidate rows assessed: {len(output)}",
        f"Rows upgraded in external UniProt sensitivity branch: {len(strict)}",
        "",
        "## Decision Counts",
    ]
    for decision, count in decision_counts.items():
        lines.append(f"- {decision}: {count}")
    lines.extend(["", "## Strict Rescue by Gene"])
    for gene, count in strict.groupby("human_gene_symbol")["scientific_name"].nunique().sort_index().items():
        lines.append(f"- {gene}: {count}")
    if strict.empty:
        lines.append("- none: 0")
    lines.extend(
        [
            "",
            "## Interpretation",
            "Strict rows here have target-species UniProt support plus reciprocal same-gene sequence support. This remains a sensitivity branch because the sequence source is external UniProt rather than the local assembly/GFF rescue path.",
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
