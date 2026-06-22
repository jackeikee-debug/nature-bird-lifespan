"""Build strict v2 scoring eligibility table for the expanded maintenance panel."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


STRICT_READY_CLASSES = {"observed_high_coverage_seed"}


def score_group(row: pd.Series) -> str:
    external_decision = str(row.get("p2_2_external_decision", ""))
    triage_class = str(row.get("triage_class", ""))
    feasibility = str(row.get("orthology_feasibility_class", ""))
    if external_decision == "exclude_from_strict_scoring_pending_domain_validation":
        return "exclude_from_absence_scoring"
    if external_decision == "manual_review_then_crossdb_confirmation":
        return "domain_required"
    if external_decision == "can_move_to_crossdb_confirmation":
        return "crossdb_confirm"
    if triage_class in {"pilot_symbol_stable", "pilot_mostly_stable"}:
        return "crossdb_confirm"
    if triage_class == "pilot_mixed_requires_crossdb":
        return "domain_required"
    if feasibility in STRICT_READY_CLASSES:
        return "strict_ready"
    if feasibility == "new_high_priority_validation_required":
        return "domain_required"
    if feasibility == "standard_mapping_candidate":
        return "standard_mapping_pending"
    if feasibility == "observed_seed_needs_review":
        return "sensitivity_only"
    return "review_required"


def absence_allowed(group: str) -> bool:
    return group not in {
        "exclude_from_absence_scoring",
        "domain_required",
        "standard_mapping_pending",
        "review_required",
    }


def strict_allowed(group: str) -> bool:
    return group == "strict_ready"


def sensitivity_allowed(group: str) -> bool:
    return group in {
        "strict_ready",
        "crossdb_confirm",
        "sensitivity_only",
        "standard_mapping_pending",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthology-audit", type=pathlib.Path, required=True)
    parser.add_argument("--pilot-triage", type=pathlib.Path, required=True)
    parser.add_argument("--external-status", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--module-summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--blocked-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    audit = pd.read_csv(args.orthology_audit, sep="\t")
    triage = pd.read_csv(args.pilot_triage, sep="\t")
    external = pd.read_csv(args.external_status, sep="\t")

    triage_cols = [
        "human_gene_symbol",
        "triage_class",
        "candidate_fraction",
        "strict_panel_ready",
        "recommended_next_action",
    ]
    external_cols = [
        "human_gene_symbol",
        "p2_2_external_decision",
        "strict_scoring_allowed_now",
        "symbol_like_hits",
        "broad_hits",
        "no_hits",
        "domain_validation_rule",
    ]
    table = audit.merge(triage[triage_cols], on="human_gene_symbol", how="left")
    table = table.merge(external[external_cols], on="human_gene_symbol", how="left")

    table["v2_scoring_group"] = table.apply(score_group, axis=1)
    table["strict_score_allowed"] = table["v2_scoring_group"].map(strict_allowed)
    table["sensitivity_score_allowed"] = table["v2_scoring_group"].map(sensitivity_allowed)
    table["absence_scoring_allowed"] = table["v2_scoring_group"].map(absence_allowed)
    table["next_validation_step"] = table["v2_scoring_group"].map(
        {
            "strict_ready": "can_enter_strict_v2_seed_supported_score",
            "crossdb_confirm": "run_cross_database_or_sequence_confirmation_before_strict_upgrade",
            "domain_required": "run_domain_or_family_level_validation_before_any_absence_call",
            "exclude_from_absence_scoring": "exclude_from_strict_absence_scoring_until_external_domain_support",
            "standard_mapping_pending": "run_standard_ncbi_or_crossdb_mapping_before_scoring",
            "sensitivity_only": "retain_only_in_sensitivity_models",
            "review_required": "manual_review_required",
        }
    )
    table["claim_use"] = table["v2_scoring_group"].map(
        {
            "strict_ready": "main_strict_score",
            "crossdb_confirm": "not_main_until_confirmed",
            "domain_required": "not_main_domain_review",
            "exclude_from_absence_scoring": "blocked_no_absence_claim",
            "standard_mapping_pending": "pending_not_main",
            "sensitivity_only": "sensitivity_only",
            "review_required": "blocked_review",
        }
    )

    sort_cols = [
        "module_order",
        "maintenance_module_v2",
        "v2_scoring_group",
        "submodule_v2",
        "gene_order_within_module",
        "human_gene_symbol",
    ]
    table = table.sort_values(sort_cols)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, sep="\t", index=False)

    module_summary = (
        table.groupby(["maintenance_module_v2", "v2_scoring_group"], as_index=False)
        .agg(genes=("human_gene_symbol", "nunique"))
        .pivot(index="maintenance_module_v2", columns="v2_scoring_group", values="genes")
        .fillna(0)
        .astype(int)
        .reset_index()
    )
    module_summary["total_genes"] = module_summary.drop(columns=["maintenance_module_v2"]).sum(axis=1)
    args.module_summary_output.parent.mkdir(parents=True, exist_ok=True)
    module_summary.to_csv(args.module_summary_output, sep="\t", index=False)

    blocked = table[
        table["v2_scoring_group"].isin(
            {"exclude_from_absence_scoring", "domain_required", "review_required"}
        )
    ].copy()
    args.blocked_output.parent.mkdir(parents=True, exist_ok=True)
    blocked.to_csv(args.blocked_output, sep="\t", index=False)

    group_counts = table["v2_scoring_group"].value_counts().sort_index()
    strict_count = int(table["strict_score_allowed"].sum())
    sensitivity_count = int(table["sensitivity_score_allowed"].sum())
    absence_count = int(table["absence_scoring_allowed"].sum())
    lines = [
        "# Phase 2 Strict v2 Scoring Eligibility Report",
        "",
        "## Summary",
        "",
        f"Genes assessed: {len(table)}",
        f"Strict-score allowed now: {strict_count}",
        f"Sensitivity-score allowed now: {sensitivity_count}",
        f"Absence-scoring allowed now: {absence_count}",
        "",
        "## Scoring Groups",
        "",
    ]
    for group, count in group_counts.items():
        lines.append(f"- {group}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Only strict_ready genes can enter the current strict v2 score without further validation. crossdb_confirm genes are promising but must be upgraded through external database or sequence confirmation before becoming strict. domain_required and exclude_from_absence_scoring genes must not be counted as biological absences.",
            "",
            "## Immediate Next Step",
            "",
            "Run cross-database confirmation for crossdb_confirm genes, then build strict and sensitivity ortholog matrices separately.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
