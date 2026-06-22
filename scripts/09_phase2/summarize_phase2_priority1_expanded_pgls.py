"""Write a Phase 2-specific summary for priority1-expanded PGLS results."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgls", type=pathlib.Path, required=True)
    parser.add_argument("--scores", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    pgls = pd.read_csv(args.pgls, sep="\t")
    scores = pd.read_csv(args.scores, sep="\t")
    trans = pgls[pgls["maintenance_module"] == "transposon_repeat_suppression"].copy()
    keep_models = {"mass_clade_module", "residual_module", "pgls_clade_residual_module"}
    trans = trans[trans["model"].isin(keep_models)].copy()
    for col in ["module_estimate", "module_p", "module_p_bh_by_variant_model", "n", "lambda"]:
        trans[col] = pd.to_numeric(trans[col], errors="coerce")
    trans = trans.sort_values(["score_variant", "model"])
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    trans.to_csv(args.summary_output, sep="\t", index=False)

    score_trans = scores[scores["maintenance_module_v2"] == "transposon_repeat_suppression"].copy()
    score_trans["coverage_fraction"] = pd.to_numeric(score_trans["coverage_fraction"], errors="coerce")
    score_summary = score_trans.groupby("score_variant", as_index=False).agg(
        genes_total=("genes_total", "first"),
        mean_coverage=("coverage_fraction", "mean"),
        min_coverage=("coverage_fraction", "min"),
        max_coverage=("coverage_fraction", "max"),
    )

    lines = [
        "# Phase 2 Priority1-Expanded Transposon PGLS Summary",
        "",
        "## Score Coverage",
        "",
    ]
    for _, row in score_summary.iterrows():
        lines.append(
            f"- {row['score_variant']}: genes_total={int(row['genes_total'])}, "
            f"mean_coverage={row['mean_coverage']:.3f}, "
            f"range={row['min_coverage']:.3f}-{row['max_coverage']:.3f}"
        )
    lines.extend(["", "## Transposon Main Effects", ""])
    for _, row in trans.iterrows():
        lines.append(
            f"- {row['score_variant']} / {row['model']}: "
            f"estimate={row['module_estimate']:.4f}, p={row['module_p']:.4g}, "
            f"BH={row['module_p_bh_by_variant_model']:.4g}, n={int(row['n'])}, "
            f"lambda={row['lambda']:.3f}"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The transposon/repeat-suppression score remains positive after expanding from the original 41-gene strict-ready baseline to the 48-gene sequence-updated strict panel. The priority1-domain sensitivity variant is also positive, but remains sensitivity-only because TDRD paralogs are protected by the paralog guard and absence claims are not allowed for the newly expanded rows.",
            "",
            "This is an interim P2.4 result: it validates the priority-1 expansion path, but it is not yet the full 200-gene sensitivity analysis.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.report}")


if __name__ == "__main__":
    main()
