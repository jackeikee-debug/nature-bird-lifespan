"""Summarize Phase 2 W2-expanded PGLS results."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


TARGET_MODULES = {
    "transposon_repeat_suppression",
    "chromatin_repression_heterochromatin",
}
KEY_MODELS = {
    "mass_clade_module",
    "residual_module",
    "pgls_clade_residual_module",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgls", type=pathlib.Path, required=True)
    parser.add_argument("--scores", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    pgls = pd.read_csv(args.pgls, sep="\t")
    scores = pd.read_csv(args.scores, sep="\t")
    for col in ["module_estimate", "module_p", "module_p_bh_by_variant_model", "n", "lambda"]:
        pgls[col] = pd.to_numeric(pgls[col], errors="coerce")
    target = pgls[
        pgls["maintenance_module"].isin(TARGET_MODULES) & pgls["model"].isin(KEY_MODELS)
    ].copy()
    target = target.sort_values(["maintenance_module", "score_variant", "model"])
    valid_target = target.dropna(subset=["module_estimate", "module_p", "n"]).copy()
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    target.to_csv(args.summary_output, sep="\t", index=False)

    score_target = scores[scores["maintenance_module_v2"].isin(TARGET_MODULES)].copy()
    score_target["coverage_fraction"] = pd.to_numeric(score_target["coverage_fraction"], errors="coerce")
    coverage = score_target.groupby(["score_variant", "maintenance_module_v2"], as_index=False).agg(
        genes_total=("genes_total", "first"),
        mean_coverage=("coverage_fraction", "mean"),
        min_coverage=("coverage_fraction", "min"),
        max_coverage=("coverage_fraction", "max"),
    )

    errors = pgls[pgls["error"].fillna("").astype(str) != ""]
    lines = [
        "# Phase 2 W2-Expanded PGLS Summary",
        "",
        "## Coverage",
        "",
    ]
    for _, row in coverage.iterrows():
        lines.append(
            f"- {row['score_variant']} / {row['maintenance_module_v2']}: "
            f"genes_total={int(row['genes_total'])}, mean_coverage={row['mean_coverage']:.3f}, "
            f"range={row['min_coverage']:.3f}-{row['max_coverage']:.3f}"
        )
    lines.extend(["", "## Key Module Effects", ""])
    for _, row in valid_target.iterrows():
        lines.append(
            f"- {row['score_variant']} / {row['maintenance_module']} / {row['model']}: "
            f"estimate={row['module_estimate']:.4f}, p={row['module_p']:.4g}, "
            f"BH={row['module_p_bh_by_variant_model']:.4g}, n={int(row['n'])}, "
            f"lambda={row['lambda']:.3f}"
        )
    skipped = target[target["n"].isna() | target["module_p"].isna()]
    if not skipped.empty:
        lines.extend(["", "## Skipped Non-informative Module Models", ""])
        for _, row in skipped.iterrows():
            lines.append(
                f"- {row['score_variant']} / {row['maintenance_module']} / {row['model']}: "
                f"not modelled because this module is absent or non-informative in that score variant "
                f"({row['error']})"
            )
    lines.extend(["", "## Model Errors", ""])
    if errors.empty:
        lines.append("- None")
    else:
        for _, row in errors.iterrows():
            lines.append(
                f"- {row['score_variant']} / {row['maintenance_module']} / {row['model']}: {row['error']}"
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The W2-expanded variant adds crossdb_confirm chromatin/repression and transposon/repeat genes to the priority1-expanded matrix. Its key use is to test whether the transposon signal survives when chromatin repression becomes an explicit competing module.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.report}")


if __name__ == "__main__":
    main()
