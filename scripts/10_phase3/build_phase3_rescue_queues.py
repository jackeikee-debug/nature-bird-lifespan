"""Build Phase 3 orthology and coverage rescue queues."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


FINAL_VARIANT = "phase2_W3_full_background_sensitivity"
TARGET_MODULES = {
    "transposon_repeat_suppression",
    "chromatin_repression_heterochromatin",
}
MISSING_STATUSES = {
    "w2_expansion_no_ncbi_gene_candidate",
    "w3_expansion_no_ncbi_gene_candidate",
    "priority1_expansion_no_ncbi_gene_candidate",
    "not_found",
    "not_found_after_diamond_validation",
    "",
}


def boolish(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def priority_reason(row: pd.Series) -> str:
    reasons = []
    if row["maintenance_module"] == "transposon_repeat_suppression":
        reasons.append("target_transposon_module")
    elif row["maintenance_module"] == "chromatin_repression_heterochromatin":
        reasons.append("target_chromatin_competing_module")
    if row["clade"] == "Aves":
        reasons.append("bird_enriched_claim_relevant")
    if row["low_transposon_coverage_species"]:
        reasons.append("low_transposon_coverage_species")
    if row["gene_family_risk"] in {"module_high_priority", "high_paralog_risk"}:
        reasons.append("high_priority_or_paralog_risk")
    if row["v2_scoring_group"] in {"crossdb_confirm", "domain_supported_paralog_guard"}:
        reasons.append("needs_crossdb_or_domain_confirmation")
    return ";".join(reasons) if reasons else "background_rescue"


def rescue_route(row: pd.Series) -> str:
    gene = row["human_gene_symbol"]
    module = row["maintenance_module"]
    if gene.startswith("TDRD"):
        return "NCBI_Protein_or_UniProt_then_DIAMOND_reciprocal_plus_Tudor_domain_guard"
    if gene in {"PIWIL1", "PIWIL2", "PIWIL3", "PIWIL4", "DDX4", "MOV10L1"}:
        return "NCBI_Protein_or_UniProt_then_DIAMOND_reciprocal_plus_family_specific_domain_check"
    if module == "transposon_repeat_suppression":
        return "NCBI_Gene_rescue_then_UniProt_OrthoDB_OMA_crosscheck"
    if module == "chromatin_repression_heterochromatin":
        return "NCBI_Gene_rescue_then_OMA_OrthoDB_Ensembl_Compara_crosscheck"
    return "NCBI_Gene_rescue_then_crossdb_review"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", type=pathlib.Path, required=True)
    parser.add_argument("--matrix", type=pathlib.Path, required=True)
    parser.add_argument("--eligibility", type=pathlib.Path, required=True)
    parser.add_argument("--low-coverage-output", type=pathlib.Path, required=True)
    parser.add_argument("--orthology-queue-output", type=pathlib.Path, required=True)
    parser.add_argument("--priority-summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--coverage-threshold", type=float, default=0.5)
    args = parser.parse_args()

    scores = pd.read_csv(args.scores, sep="\t")
    matrix = pd.read_csv(args.matrix, sep="\t")
    eligibility = pd.read_csv(args.eligibility, sep="\t")

    scores = scores[scores["score_variant"].eq(FINAL_VARIANT)].copy()
    coverage_rows = []
    module_names = sorted(scores["maintenance_module_v2"].unique())
    for _, row in scores.iterrows():
        cov = pd.to_numeric(pd.Series([row["coverage_fraction"]]), errors="coerce").iloc[0]
        coverage_rows.append(
            {
                "scientific_name": row["scientific_name"],
                "clade": row["clade"],
                "flight_status": row["flight_status"],
                "genome_analysis_tier": row["genome_analysis_tier"],
                "maintenance_module": row["maintenance_module_v2"],
                "coverage_fraction": cov,
                "confidence_weighted_score": row["confidence_weighted_score"],
                "genes_total": row["genes_total"],
                "below_threshold": bool(cov < args.coverage_threshold),
                "coverage_deficit_to_threshold": max(0.0, args.coverage_threshold - cov),
            }
        )
    coverage = pd.DataFrame(coverage_rows)
    low_cov = coverage[coverage["below_threshold"]].copy()
    low_cov["rescue_priority"] = low_cov.apply(
        lambda r: 1
        if r["maintenance_module"] == "transposon_repeat_suppression" and r["clade"] == "Aves"
        else 2
        if r["maintenance_module"] in TARGET_MODULES
        else 3,
        axis=1,
    )
    low_cov = low_cov.sort_values(
        ["rescue_priority", "coverage_fraction", "maintenance_module", "clade", "scientific_name"]
    )

    low_trans_species = set(
        low_cov.loc[
            low_cov["maintenance_module"].eq("transposon_repeat_suppression"),
            "scientific_name",
        ]
    )

    elig_cols = [
        "human_gene_symbol",
        "maintenance_module_v2",
        "submodule_v2",
        "v2_scoring_group",
        "strict_score_allowed",
        "sensitivity_score_allowed",
        "absence_scoring_allowed",
        "gene_family_risk",
        "claim_use",
    ]
    elig = eligibility[[c for c in elig_cols if c in eligibility.columns]].copy()
    matrix = matrix.merge(
        elig,
        left_on=["human_gene_symbol", "maintenance_module"],
        right_on=["human_gene_symbol", "maintenance_module_v2"],
        how="left",
    )
    matrix["final_candidate_status"] = matrix["final_candidate_status"].fillna("")
    matrix["is_missing_or_unresolved"] = matrix["final_candidate_status"].isin(MISSING_STATUSES)
    matrix["low_transposon_coverage_species"] = matrix["scientific_name"].isin(low_trans_species)
    matrix["strict_score_allowed_bool"] = matrix.get("strict_score_allowed", False).apply(boolish)
    matrix["sensitivity_score_allowed_bool"] = matrix.get("sensitivity_score_allowed", False).apply(boolish)
    matrix["absence_scoring_allowed_bool"] = matrix.get("absence_scoring_allowed", False).apply(boolish)

    rescue = matrix[
        matrix["is_missing_or_unresolved"]
        & (
            matrix["maintenance_module"].isin(TARGET_MODULES)
            | matrix["low_transposon_coverage_species"]
            | matrix["strict_score_allowed_bool"]
            | matrix["sensitivity_score_allowed_bool"]
        )
    ].copy()
    rescue["priority_reason"] = rescue.apply(priority_reason, axis=1)
    rescue["recommended_rescue_route"] = rescue.apply(rescue_route, axis=1)
    rescue["rescue_priority"] = rescue.apply(
        lambda r: 1
        if r["maintenance_module"] == "transposon_repeat_suppression"
        and r["clade"] == "Aves"
        and r["low_transposon_coverage_species"]
        else 2
        if r["maintenance_module"] == "transposon_repeat_suppression"
        else 3
        if r["maintenance_module"] == "chromatin_repression_heterochromatin"
        else 4,
        axis=1,
    )
    rescue_cols = [
        "rescue_priority",
        "priority_reason",
        "recommended_rescue_route",
        "scientific_name",
        "clade",
        "flight_status",
        "species_taxid",
        "best_assembly_accession",
        "genome_analysis_tier",
        "maintenance_module",
        "human_gene_symbol",
        "submodule_v2",
        "v2_scoring_group",
        "gene_family_risk",
        "final_candidate_status",
        "final_candidate_source",
        "final_candidate_confidence",
        "ortholog_query_status",
        "ortholog_gene_id",
        "ortholog_gene_symbol",
        "strict_score_allowed",
        "sensitivity_score_allowed",
        "absence_scoring_allowed",
        "claim_use",
    ]
    rescue = rescue[[c for c in rescue_cols if c in rescue.columns]].sort_values(
        ["rescue_priority", "clade", "scientific_name", "maintenance_module", "human_gene_symbol"]
    )

    module_summary = coverage.groupby("maintenance_module", as_index=False).agg(
        species=("scientific_name", "nunique"),
        mean_coverage=("coverage_fraction", "mean"),
        min_coverage=("coverage_fraction", "min"),
        low_coverage_species=("below_threshold", "sum"),
    )
    rescue_summary = rescue.groupby("maintenance_module", as_index=False).agg(
        rescue_rows=("human_gene_symbol", "count"),
        rescue_genes=("human_gene_symbol", "nunique"),
        rescue_species=("scientific_name", "nunique"),
    )
    summary = module_summary.merge(rescue_summary, on="maintenance_module", how="left").fillna(
        {"rescue_rows": 0, "rescue_genes": 0, "rescue_species": 0}
    )
    summary["phase3_priority"] = summary["maintenance_module"].apply(
        lambda m: "primary" if m == "transposon_repeat_suppression" else "secondary" if m in TARGET_MODULES else "background"
    )
    summary = summary.sort_values(["phase3_priority", "mean_coverage"])

    args.low_coverage_output.parent.mkdir(parents=True, exist_ok=True)
    args.orthology_queue_output.parent.mkdir(parents=True, exist_ok=True)
    args.priority_summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    low_cov.to_csv(args.low_coverage_output, sep="\t", index=False)
    rescue.to_csv(args.orthology_queue_output, sep="\t", index=False)
    summary.to_csv(args.priority_summary_output, sep="\t", index=False)

    priority1 = rescue[rescue["rescue_priority"].eq(1)]
    low_trans = low_cov[low_cov["maintenance_module"].eq("transposon_repeat_suppression")]
    lines = [
        "# Phase 3 Rescue Queue Report",
        "",
        "## Scope",
        "",
        f"Final score variant: `{FINAL_VARIANT}`",
        f"Low-coverage threshold: {args.coverage_threshold}",
        "",
        "## Coverage Summary",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            "- {module}: mean coverage={mean:.3f}, min={minv:.3f}, low-coverage species={low}, rescue rows={rows}".format(
                module=row["maintenance_module"],
                mean=row["mean_coverage"],
                minv=row["min_coverage"],
                low=int(row["low_coverage_species"]),
                rows=int(row["rescue_rows"]),
            )
        )
    lines.extend(
        [
            "",
            "## Priority 1 Rescue Target",
            "",
            f"Priority-1 transposon/bird low-coverage rows: {len(priority1)}",
            f"Low-coverage transposon species: {low_trans['scientific_name'].nunique()}",
            "",
            "Top low-coverage transposon species:",
            "",
        ]
    )
    for _, row in low_trans.head(12).iterrows():
        lines.append(
            f"- {row['scientific_name']} ({row['clade']}): coverage={row['coverage_fraction']:.3f}"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Phase 3 should first rescue transposon/repeat rows in low-coverage bird species, then cross-check chromatin/repression rows as the closest competing module. No low-coverage row should be interpreted as a true absence until sequence or cross-database confirmation is complete.",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.low_coverage_output}, {args.orthology_queue_output}, and {args.report}")


if __name__ == "__main__":
    main()
