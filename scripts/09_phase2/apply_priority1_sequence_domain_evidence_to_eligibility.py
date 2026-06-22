"""Overlay priority-1 sequence/domain evidence onto strict v2 eligibility."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


STRICT_SEQUENCE_GROUP = "strict_sequence_supported"
DOMAIN_PARALOG_GROUP = "domain_supported_paralog_guard"
DOMAIN_MANUAL_GROUP = "domain_supported_manual_upgrade_candidate"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-eligibility", type=pathlib.Path, required=True)
    parser.add_argument("--diamond-strict", type=pathlib.Path, required=True)
    parser.add_argument("--domain-summary", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--module-summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--blocked-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    table = pd.read_csv(args.base_eligibility, sep="\t")
    strict = pd.read_csv(args.diamond_strict, sep="\t")
    domain = pd.read_csv(args.domain_summary, sep="\t")

    strict_genes = set(strict["human_gene_symbol"].dropna())
    domain_supported = domain[
        domain["sensitivity_upgrade_allowed"].astype(bool)
        & domain["all_species_domain_rule_passed"].astype(bool)
    ].copy()
    domain_paralog_genes = set(
        domain_supported[
            domain_supported["domain_confirmation_recommendation"]
            == "domain_supported_keep_paralog_guard"
        ]["human_gene_symbol"]
    )
    domain_manual_genes = set(
        domain_supported[
            domain_supported["domain_confirmation_recommendation"]
            == "domain_supported_manual_queue_upgrade_candidate"
        ]["human_gene_symbol"]
    )

    table["priority1_sequence_evidence"] = "none"
    table["priority1_domain_evidence"] = "none"
    table["evidence_overlay_applied"] = False
    table["base_v2_scoring_group"] = table["v2_scoring_group"]
    table["base_next_validation_step"] = table["next_validation_step"]
    table["base_claim_use"] = table["claim_use"]

    strict_mask = table["human_gene_symbol"].isin(strict_genes)
    table.loc[strict_mask, "v2_scoring_group"] = STRICT_SEQUENCE_GROUP
    table.loc[strict_mask, "strict_score_allowed"] = True
    table.loc[strict_mask, "sensitivity_score_allowed"] = True
    table.loc[strict_mask, "absence_scoring_allowed"] = False
    table.loc[
        strict_mask, "next_validation_step"
    ] = "enter_strict_v2_sequence_supported_score_no_absence_claim"
    table.loc[strict_mask, "claim_use"] = "main_strict_score_sequence_supported"
    table.loc[
        strict_mask, "priority1_sequence_evidence"
    ] = "NCBI_Gene_plus_DIAMOND_BLASTP_reciprocal_same_gene_support"
    table.loc[strict_mask, "evidence_overlay_applied"] = True

    paralog_mask = table["human_gene_symbol"].isin(domain_paralog_genes)
    table.loc[paralog_mask, "v2_scoring_group"] = DOMAIN_PARALOG_GROUP
    table.loc[paralog_mask, "strict_score_allowed"] = False
    table.loc[paralog_mask, "sensitivity_score_allowed"] = True
    table.loc[paralog_mask, "absence_scoring_allowed"] = False
    table.loc[
        paralog_mask, "next_validation_step"
    ] = "retain_sensitivity_domain_supported_but_require_tree_or_HMM_for_strict_TDRD_paralog_resolution"
    table.loc[paralog_mask, "claim_use"] = "domain_supported_sensitivity_only_paralog_guard"
    table.loc[paralog_mask, "priority1_domain_evidence"] = "InterProScan_domain_rule_passed"
    table.loc[paralog_mask, "evidence_overlay_applied"] = True

    manual_mask = table["human_gene_symbol"].isin(domain_manual_genes)
    table.loc[manual_mask, "v2_scoring_group"] = DOMAIN_MANUAL_GROUP
    table.loc[manual_mask, "strict_score_allowed"] = False
    table.loc[manual_mask, "sensitivity_score_allowed"] = True
    table.loc[manual_mask, "absence_scoring_allowed"] = False
    table.loc[
        manual_mask, "next_validation_step"
    ] = "retain_sensitivity_domain_supported_manual_queue_pending_stricter_orthology_evidence"
    table.loc[manual_mask, "claim_use"] = "domain_supported_sensitivity_only"
    table.loc[manual_mask, "priority1_domain_evidence"] = "InterProScan_domain_rule_passed"
    table.loc[manual_mask, "evidence_overlay_applied"] = True

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
        (~table["strict_score_allowed"].astype(bool))
        & (~table["sensitivity_score_allowed"].astype(bool))
    ].copy()
    args.blocked_output.parent.mkdir(parents=True, exist_ok=True)
    blocked.to_csv(args.blocked_output, sep="\t", index=False)

    group_counts = table["v2_scoring_group"].value_counts().sort_index()
    lines = [
        "# Phase 2 Sequence-Updated Strict v2 Scoring Eligibility Report",
        "",
        "## Summary",
        "",
        f"Genes assessed: {table['human_gene_symbol'].nunique()}",
        f"Strict-score allowed now: {int(table['strict_score_allowed'].sum())}",
        f"Sensitivity-score allowed now: {int(table['sensitivity_score_allowed'].sum())}",
        f"Absence-scoring allowed now: {int(table['absence_scoring_allowed'].sum())}",
        f"DIAMOND/BLAST strict-upgrade genes overlaid: {len(strict_genes)}",
        f"Domain-supported paralog-guard genes overlaid: {len(domain_paralog_genes)}",
        f"Domain-supported manual-queue genes overlaid: {len(domain_manual_genes)}",
        "",
        "## Overlay Groups",
        "",
    ]
    for group, count in group_counts.items():
        lines.append(f"- {group}: {count}")
    lines.extend(
        [
            "",
            "## Strict Sequence-Supported Upgrades",
            "",
        ]
    )
    for gene in sorted(strict_genes):
        lines.append(f"- {gene}")
    lines.extend(["", "## Domain-Supported Sensitivity Genes", ""])
    for gene in sorted(domain_paralog_genes | domain_manual_genes):
        lines.append(f"- {gene}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The 7 DIAMOND/BLAST reciprocal-support genes can enter the strict v2 presence score, but they remain protected from absence claims. TDRD-family domain-supported genes remain under paralog guard and should be used only in sensitivity analyses until tree/HMM discrimination resolves exact paralogs.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
