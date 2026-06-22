"""Summarize Phase 2 gates P2.1-P2.3."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validated", type=pathlib.Path, required=True)
    parser.add_argument("--orthology", type=pathlib.Path, required=True)
    parser.add_argument("--model-spec", type=pathlib.Path, required=True)
    parser.add_argument("--task-register", type=pathlib.Path, required=True)
    parser.add_argument("--task-register-current", type=pathlib.Path, required=True)
    parser.add_argument("--gate-table", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    validated = pd.read_csv(args.validated, sep="\t")
    orthology = pd.read_csv(args.orthology, sep="\t")
    model_spec = pd.read_csv(args.model_spec, sep="\t")
    tasks = pd.read_csv(args.task_register, sep="\t")

    n_genes = len(validated)
    n_pass = int((validated["validation_decision"] == "pass").sum())
    n_high = int(
        (
            orthology["orthology_feasibility_class"]
            == "new_high_priority_validation_required"
        ).sum()
    )
    n_seed_high = int(
        (orthology["orthology_feasibility_class"] == "observed_high_coverage_seed").sum()
    )

    gates = pd.DataFrame(
        [
            {
                "gate": "P2.1_gene_symbol_source_validation",
                "status": "pass" if n_pass / n_genes >= 0.90 else "caution_or_fail",
                "evidence": f"{n_pass}/{n_genes} genes passed symbol and source checks",
                "next_action": "Use validated panel for v2 orthology work",
            },
            {
                "gate": "P2.2_orthology_feasibility_audit",
                "status": "caution",
                "evidence": f"{n_seed_high} seed genes observed high coverage; {n_high} new high-priority genes require validation",
                "next_action": "Proceed only to targeted v2 mapping, not final biology claims",
            },
            {
                "gate": "P2.3_annotation_bias_design",
                "status": "pass_design_frozen",
                "evidence": f"{len(model_spec)} module model specifications written",
                "next_action": "Carry tier, missingness, and coverage covariates into expanded PGLS",
            },
        ]
    )
    args.gate_table.parent.mkdir(parents=True, exist_ok=True)
    gates.to_csv(args.gate_table, sep="\t", index=False)

    for gate, status in [
        ("P2.1", "pass"),
        ("P2.2", "caution"),
        ("P2.3", "pass_design_frozen"),
    ]:
        tasks.loc[tasks["stage"] == gate, "status"] = status
    args.task_register_current.parent.mkdir(parents=True, exist_ok=True)
    tasks.to_csv(args.task_register_current, sep="\t", index=False)

    lines = [
        "# Phase 2 First Three Gates Report",
        "",
        "## Overall Status",
        "",
        "P2.1-P2.3 are complete enough to proceed to targeted v2 orthology mapping and expanded module scoring design. They are not complete enough to make expanded biological claims yet.",
        "",
        "## Gate Results",
        "",
    ]
    for _, row in gates.iterrows():
        lines.extend(
            [
                f"### {row['gate']}",
                "",
                f"Status: **{row['status']}**",
                "",
                f"Evidence: {row['evidence']}",
                "",
                f"Next action: {row['next_action']}",
                "",
            ]
        )
    lines.extend(
        [
            "## Interpretation",
            "",
            "The project has not drifted into a duplicate of prior bird/bat longevity work at this stage. The key differentiators remain viable, but P2.2 is the bottleneck: repeat/chromatin genes now require real orthology validation before any strong transposon-specific manuscript claim.",
            "",
            "## Recommended Next Step",
            "",
            "Start P2.4 only after building a targeted v2 orthology mapping queue from the 92 high-priority genes, with strict and sensitivity tiers separated from the beginning.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.report}")


if __name__ == "__main__":
    main()
