"""Export strict, sensitivity, and absence-eligible v2 scoring gene panels."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


PANEL_COLUMNS = [
    "human_gene_symbol",
    "maintenance_module_v2",
    "submodule_v2",
    "module_order",
    "gene_order_within_module",
    "seed_status",
    "inclusion_tier",
    "gene_family_risk",
    "orthology_feasibility_class",
    "triage_class",
    "p2_2_external_decision",
    "v2_scoring_group",
    "next_validation_step",
    "claim_use",
]


def write_panel(table: pd.DataFrame, flag: str, output: pathlib.Path) -> pd.DataFrame:
    panel = table.loc[table[flag].astype(bool), PANEL_COLUMNS].copy()
    panel = panel.sort_values(
        [
            "module_order",
            "maintenance_module_v2",
            "submodule_v2",
            "gene_order_within_module",
            "human_gene_symbol",
        ]
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    panel.to_csv(output, sep="\t", index=False)
    return panel


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eligibility", type=pathlib.Path, required=True)
    parser.add_argument("--strict-output", type=pathlib.Path, required=True)
    parser.add_argument("--sensitivity-output", type=pathlib.Path, required=True)
    parser.add_argument("--absence-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    table = pd.read_csv(args.eligibility, sep="\t")
    strict = write_panel(table, "strict_score_allowed", args.strict_output)
    sensitivity = write_panel(table, "sensitivity_score_allowed", args.sensitivity_output)
    absence = write_panel(table, "absence_scoring_allowed", args.absence_output)

    lines = [
        "# Phase 2 v2 Gene Scoring Panels Report",
        "",
        "## Summary",
        "",
        f"Eligibility genes assessed: {table['human_gene_symbol'].nunique()}",
        f"Strict v2 panel genes: {strict['human_gene_symbol'].nunique()}",
        f"Sensitivity v2 panel genes: {sensitivity['human_gene_symbol'].nunique()}",
        f"Absence-allowed v2 panel genes: {absence['human_gene_symbol'].nunique()}",
        "",
        "## Use",
        "",
        "The strict panel is the only panel currently eligible for main v2 scoring. The sensitivity panel can be used for exploratory or robustness analyses. The absence-allowed panel defines genes for which missing ortholog calls may be interpreted after additional row-level evidence checks; domain-required and excluded genes remain protected from absence claims.",
        "",
        "## Inputs",
        "",
        f"- eligibility: {args.eligibility.as_posix()}",
    ]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.strict_output}")


if __name__ == "__main__":
    main()
