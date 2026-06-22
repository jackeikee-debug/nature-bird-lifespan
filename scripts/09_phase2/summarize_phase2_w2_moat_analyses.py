"""Summarize the three Phase 2 W2 moat analyses."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def fmt(value: object, digits: int = 4) -> str:
    try:
        val = float(value)
    except (TypeError, ValueError):
        return "NA"
    if pd.isna(val):
        return "NA"
    return f"{val:.{digits}g}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation", type=pathlib.Path, required=True)
    parser.add_argument("--random", type=pathlib.Path, required=True)
    parser.add_argument("--submodule", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    ann = pd.read_csv(args.annotation, sep="\t")
    rnd = pd.read_csv(args.random, sep="\t")
    sub = pd.read_csv(args.submodule, sep="\t")
    for frame in [ann, rnd, sub]:
        for col in frame.columns:
            if col in {"p", "p_bh", "estimate", "observed_r", "observed_r2", "empirical_p_r_greater", "empirical_p_r2_greater"}:
                frame[col] = pd.to_numeric(frame[col], errors="coerce")

    lines = [
        "# Phase 2 W2 Moat Analyses Summary",
        "",
        "## Overall Interpretation",
        "",
        "The three moat analyses support a cautious repeat/chromatin-axis interpretation, not a transposon-specific independent-effect claim. The signal survives genome-tier adjustment and is stronger than W2 random maintenance-gene sets by correlation/R2, but it is highly entangled with module coverage and distributed across related repeat/chromatin submodules.",
        "",
        "## 1. Annotation-Bias Models",
        "",
    ]
    for module in ["transposon_repeat_suppression", "chromatin_repression_heterochromatin"]:
        rows = ann[ann["maintenance_module"] == module].set_index("model")
        base = rows.loc["residual_base"]
        tier = rows.loc["residual_tier"]
        coverage = rows.loc["residual_coverage"]
        lines.append(
            f"- {module}: base PGLS residual estimate={fmt(base['estimate'])}, p={fmt(base['p'])}; "
            f"tier-adjusted estimate={fmt(tier['estimate'])}, p={fmt(tier['p'])}; "
            f"coverage-adjusted estimate={fmt(coverage['estimate'])}, p={fmt(coverage['p'])}."
        )
    lines.extend(
        [
            "",
            "Interpretation: genome tier does not explain away the signal, but module coverage does. Because the score is confidence-weighted presence divided by gene count, coverage adjustment is an intentionally severe test and indicates that annotation/observability remains the leading vulnerability.",
            "",
            "## 2. Matched Random Gene-Set Tests",
            "",
        ]
    )
    for _, row in rnd.iterrows():
        lines.append(
            f"- {row['target_module']}: observed r={fmt(row['observed_r'])}, R2={fmt(row['observed_r2'])}; "
            f"empirical p(r)={fmt(row['empirical_p_r_greater'])}, p(R2)={fmt(row['empirical_p_r2_greater'])}; "
            f"target mean coverage={fmt(row['target_mean_gene_coverage'])}, nearest-random mean coverage={fmt(row['matched_mean_coverage'])}."
        )
    lines.extend(
        [
            "",
            "Interpretation: target repeat/chromatin modules outperform nearest matched W2 random maintenance-gene sets by correlation/R2. However, the non-repeat/chromatin W2 background has higher coverage, so this is supportive but not yet a final matched-control proof.",
            "",
            "## 3. Submodule Split",
            "",
        ]
    )
    sub_valid = sub[(sub["error"].fillna("") == "") & sub["submodule_v2"].notna()].copy()
    sub_valid = sub_valid[~sub_valid["submodule_v2"].isin(["transposon_repeat_suppression", "chromatin_repression_heterochromatin"])]
    sub_valid = sub_valid.sort_values("p").head(10)
    for _, row in sub_valid.iterrows():
        lines.append(
            f"- {row['submodule_v2']} / {row['model']}: estimate={fmt(row['estimate'])}, "
            f"p={fmt(row['p'])}, BH={fmt(row['p_bh'])}."
        )
    lines.extend(
        [
            "",
            "Interpretation: signals are distributed across piRNA/repeat-control and chromatin-repression submodules. This pattern favors a broader repeat/chromatin maintenance axis over a single-gene or single-submodule artifact.",
            "",
            "## Decision",
            "",
            "**Caution-go.** Continue Phase 2, but manuscript language should remain `repeat/chromatin maintenance axis` until coverage-matched controls can be improved with a larger W3/full200 gene universe or independent orthology evidence.",
        ]
    )

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.report}")


if __name__ == "__main__":
    main()
