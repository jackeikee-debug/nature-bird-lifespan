"""Summarize family/domain risks for Phase 3 priority-1 rescue rows."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


PARALOG_FAMILIES = {"DNMT1", "DNMT3A", "DNMT3B", "MBD2", "MBD3", "SETDB2"}


def num(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def classify(row: pd.Series) -> tuple[str, str, str]:
    gene = row.get("human_gene_symbol", "")
    route = row.get("next_rescue_route", "")
    decision = row.get("phase3_gff_sequence_decision", "")
    validation = row.get("validation_calls", "")
    top_refs = row.get("top_reference_genes", "")
    reciprocal = row.get("reciprocal_best_genes", "")
    protein_len = num(row.get("protein_length", ""))
    scov = num(row.get("max_scovhsp", ""))
    partial = str(row.get("partial", "")).lower() == "true"
    uniprot_call = row.get("uniprot_crosscheck_call", "")

    if row.get("can_count_as_strict_sequence_after_gff_sequence", "") == "True":
        return (
            "strict_sequence_supported",
            "Validated by reciprocal same-gene sequence search; UniProt/orthogroup evidence can be used as independent support.",
            "retain_as_strict_rescue",
        )

    if route == "external_db_or_assembly_reannotation":
        if uniprot_call in {"uniprot_gene_exact_support", "uniprot_symbol_like_support", "uniprot_product_support"}:
            return (
                "gff_absent_but_uniprot_supported",
                "Local GFF has no target hit, but UniProt suggests an annotation exists elsewhere.",
                "fetch_uniprot_sequence_or_updated_annotation_then_reciprocal_validate",
            )
        return (
            "gff_absent_external_only",
            "No local GFF hit and no strong UniProt target support; OrthoDB can only provide family-level context.",
            "query_updated_assembly_annotation_or_expand_to_orthodb_family_level",
        )

    if route == "cds_or_genome_sequence_extraction":
        return (
            "gff_annotation_without_protein_id",
            "GFF target annotation exists but no protein_id was available for sequence validation.",
            "extract_cds_or_genomic_region_and_translate_before_reciprocal_validation",
        )

    if gene == "DNMT3B" and "DNMT3A" in top_refs:
        return (
            "dnmt3b_dnmt3a_crosshit_probable_paralog",
            "DNMT3B candidate best human reference is DNMT3A; this is a high-risk DNMT-family ambiguity.",
            "do_not_score_as_strict; require domain_architecture_and_gene_tree_review",
        )

    if gene == "DNMT1" and decision in {"gff_sequence_weak_same_gene_not_strict", "gff_sequence_not_validated"}:
        if partial or protein_len < 700 or scov < 50:
            return (
                "dnmt1_partial_fragment_weak_support",
                "DNMT1 same-gene signal exists but candidate is partial/short or has low subject coverage.",
                "extract_full_length_cds_or_updated_protein_then_revalidate",
            )
        return (
            "dnmt1_manual_sequence_review",
            "DNMT1 is not strict despite same-gene evidence; inspect reciprocal hit table and isoform choice.",
            "manual_reciprocal_and_isoform_review",
        )

    if gene == "MBD2" and decision == "gff_sequence_weak_same_gene_not_strict":
        if "MBD3" in reciprocal or scov < 50 or protein_len < 150:
            return (
                "mbd2_mbd3_short_fragment_ambiguity",
                "MBD2 candidate is short/partial and MBD-family reciprocal evidence is ambiguous.",
                "do_not_score_as_strict; extract_full_length_sequence_and_run_MBD_domain_gene_tree_review",
            )
        return (
            "mbd2_manual_family_review",
            "MBD2 weak support requires family-aware review before scoring.",
            "manual_MBD_family_review",
        )

    if gene == "MBD3" and decision == "gff_sequence_weak_same_gene_not_strict":
        return (
            "mbd3_partial_fragment_weak_support",
            "MBD3 same-gene signal is weak, usually because the candidate is short/partial.",
            "extract_full_length_sequence_and_revalidate",
        )

    if route == "probable_paralog_or_reference_gap":
        return (
            "probable_paralog_or_reference_gap",
            "Best hits do not support the requested same-gene relationship.",
            "do_not_score_as_strict_unless_orthogroup_and_domain_review_resolves",
        )

    if route == "manual_reciprocal_review" or "forward_supported" in validation:
        return (
            "manual_reciprocal_review",
            "Forward support exists but reciprocal support is incomplete.",
            "inspect_all_hits_and_reference_isoforms",
        )

    if gene in PARALOG_FAMILIES:
        return (
            "paralog_family_manual_review",
            "Non-strict row belongs to a paralog-prone family and should not be treated as absence.",
            "family_specific_domain_or_gene_tree_review",
        )

    return (
        "non_strict_general_followup",
        "Non-strict row outside the main paralog-prone families.",
        "use_external_annotation_or_sequence_extraction_before_scoring",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--route-audit", type=pathlib.Path, required=True)
    parser.add_argument("--sequence-decisions", type=pathlib.Path, required=True)
    parser.add_argument("--uniprot", type=pathlib.Path, required=True)
    parser.add_argument("--crossdb", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    route = pd.read_csv(args.route_audit, sep="\t", dtype=str).fillna("")
    decisions = pd.read_csv(args.sequence_decisions, sep="\t", dtype=str).fillna("")
    uniprot = pd.read_csv(args.uniprot, sep="\t", dtype=str).fillna("")
    crossdb = pd.read_csv(args.crossdb, sep="\t", dtype=str).fillna("")
    keys = ["scientific_name", "human_gene_symbol"]
    merged = route.merge(
        decisions[
            keys
            + [
                "protein_length",
                "partial",
                "cds_product",
                "max_pident",
                "max_qcovhsp",
                "max_scovhsp",
                "phase3_gff_sequence_decision",
                "can_count_as_strict_sequence_after_gff_sequence",
            ]
        ],
        on=keys,
        how="left",
        suffixes=("", "_decision"),
    )
    for col in [
        "protein_length",
        "partial",
        "cds_product",
        "max_pident",
        "max_qcovhsp",
        "max_scovhsp",
        "phase3_gff_sequence_decision",
        "can_count_as_strict_sequence_after_gff_sequence",
    ]:
        alt = f"{col}_decision"
        if alt in merged.columns:
            merged[col] = merged[col].where(merged[col] != "", merged[alt])
    merged = merged.merge(
        uniprot[keys + ["uniprot_crosscheck_call", "uniprot_accession_top", "uniprot_gene_names_top", "uniprot_protein_name_top"]],
        on=keys,
        how="left",
    )
    merged = merged.merge(
        crossdb[keys + ["orthodb_call", "orthodb_group_id", "orthodb_group_name", "orthodb_level", "ensembl_species_call", "oma_call"]],
        on=keys,
        how="left",
    ).fillna("")

    classified = merged.apply(classify, axis=1, result_type="expand")
    merged["family_domain_review_class"] = classified[0]
    merged["family_domain_review_reason"] = classified[1]
    merged["recommended_domain_next_step"] = classified[2]

    keep = [
        "scientific_name",
        "human_gene_symbol",
        "next_rescue_route",
        "phase3_gff_sequence_decision",
        "can_count_as_strict_sequence_after_gff_sequence",
        "family_domain_review_class",
        "family_domain_review_reason",
        "recommended_domain_next_step",
        "protein_id",
        "protein_length",
        "partial",
        "cds_product",
        "validation_calls",
        "top_reference_genes",
        "reciprocal_best_genes",
        "max_pident",
        "max_qcovhsp",
        "max_scovhsp",
        "uniprot_crosscheck_call",
        "uniprot_accession_top",
        "orthodb_call",
        "orthodb_group_id",
        "orthodb_level",
        "oma_call",
        "ensembl_species_call",
    ]
    for col in keep:
        if col not in merged.columns:
            merged[col] = ""
    out = merged[keep].sort_values(["family_domain_review_class", "scientific_name", "human_gene_symbol"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    summary = (
        out.groupby(["family_domain_review_class", "recommended_domain_next_step"], as_index=False)
        .agg(rows=("human_gene_symbol", "count"), species=("scientific_name", "nunique"), genes=("human_gene_symbol", "nunique"))
        .sort_values(["family_domain_review_class", "recommended_domain_next_step"])
    )
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary_output, sep="\t", index=False)

    counts = out["family_domain_review_class"].value_counts().sort_index()
    lines = [
        "# Phase 3 Family/Domain Review Report",
        "",
        f"Rows reviewed: {len(out)}",
        "",
        "## Review Classes",
    ]
    for klass, count in counts.items():
        lines.append(f"- {klass}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "The dominant non-strict risks are not random missingness. DNMT3B rows frequently cross-hit DNMT3A, MBD2 rows are often short partial fragments with MBD-family ambiguity, and weak DNMT1 rows are mostly partial fragments. These rows should not be scored as absence or strict rescue until full-length sequence extraction, domain architecture checks, or gene-tree review resolves the ambiguity.",
            "",
            "## Outputs",
            f"- row table: `{args.output}`",
            f"- summary: `{args.summary_output}`",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
