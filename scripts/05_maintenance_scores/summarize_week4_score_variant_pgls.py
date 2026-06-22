"""Summarize Week 4 score-variant PGLS robustness."""

from __future__ import annotations

import argparse
import math
import pathlib

import pandas as pd


DEFAULT_PGLS = pathlib.Path("results/tables/week4_score_variant_pgls.tsv")
DEFAULT_SCORES = pathlib.Path("data/processed/maintenance_scores_week4_variants.tsv")
DEFAULT_OUTPUT = pathlib.Path("results/tables/week4_score_variant_pgls_summary.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week4_score_variant_pgls_summary_report.md")

KEY_MODEL = "mass_clade_module"
FOUNDATION_VARIANTS = ["all_validated", "no_protein_rescue", "ncbi_only"]


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def classify(row: dict[str, object]) -> str:
    foundation_positive = int(row["foundation_positive_significant"])
    all_p = safe_float(row.get("all_validated_p"))
    ncbi_p = safe_float(row.get("ncbi_only_p"))
    no_protein_p = safe_float(row.get("no_protein_rescue_p"))
    high_mean = safe_float(row.get("high_confidence_only_mean_score"))
    if foundation_positive == 3:
        if high_mean < 0.05:
            return "source_controlled_high_conf_sparse"
        return "robust"
    if foundation_positive == 2 and all_p < 0.05 and (ncbi_p < 0.05 or no_protein_p < 0.05):
        return "moderate_source_controlled"
    if foundation_positive >= 1 and all_p < 0.05:
        return "rescue_sensitive"
    return "exploratory"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgls", type=pathlib.Path, default=DEFAULT_PGLS)
    parser.add_argument("--scores", type=pathlib.Path, default=DEFAULT_SCORES)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    pgls = pd.read_csv(args.pgls, sep="\t")
    scores = pd.read_csv(args.scores, sep="\t")

    means = (
        scores.assign(mean_score=pd.to_numeric(scores["confidence_weighted_score"], errors="coerce"))
        .groupby(["score_variant", "maintenance_module"], as_index=False)["mean_score"]
        .mean()
    )
    key = pgls[(pgls["model"] == KEY_MODEL) & (pgls["error"].fillna("") == "")].copy()
    key["module_p"] = pd.to_numeric(key["module_p"], errors="coerce")
    key["module_estimate"] = pd.to_numeric(key["module_estimate"], errors="coerce")
    key["module_p_bh_by_variant_model"] = pd.to_numeric(
        key["module_p_bh_by_variant_model"], errors="coerce"
    )

    rows = []
    for module in sorted(key["maintenance_module"].unique()):
        row: dict[str, object] = {"maintenance_module": module, "model": KEY_MODEL}
        foundation_hits = 0
        for variant in sorted(key["score_variant"].unique()):
            sub = key[(key["maintenance_module"] == module) & (key["score_variant"] == variant)]
            mean_sub = means[(means["maintenance_module"] == module) & (means["score_variant"] == variant)]
            if sub.empty:
                continue
            rec = sub.iloc[0]
            row[f"{variant}_estimate"] = rec["module_estimate"]
            row[f"{variant}_p"] = rec["module_p"]
            row[f"{variant}_bh"] = rec["module_p_bh_by_variant_model"]
            if not mean_sub.empty:
                row[f"{variant}_mean_score"] = mean_sub.iloc[0]["mean_score"]
            if (
                variant in FOUNDATION_VARIANTS
                and rec["module_estimate"] > 0
                and rec["module_p_bh_by_variant_model"] < 0.05
            ):
                foundation_hits += 1
        row["foundation_positive_significant"] = foundation_hits
        row["week4_priority"] = classify(row)
        rows.append(row)

    out = pd.DataFrame(rows)
    priority_order = {
        "robust": 0,
        "source_controlled_high_conf_sparse": 1,
        "moderate_source_controlled": 2,
        "rescue_sensitive": 3,
        "exploratory": 4,
    }
    out["_priority_order"] = out["week4_priority"].map(priority_order).fillna(9)
    out = out.sort_values(["_priority_order", "all_validated_p", "maintenance_module"]).drop(
        columns=["_priority_order"]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    top = out.head(7)
    lines = [
        "# Week 4 Score-Variant PGLS Summary Report",
        "",
        f"Key model: `{KEY_MODEL}`.",
        "",
        "Foundation variants are `all_validated`, `no_protein_rescue`, and `ncbi_only`. A module is source-controlled when the positive mass+clade PGLS signal survives these variants after BH correction within each variant/model.",
        "",
        "## Module Priorities",
        "",
    ]
    for _, rec in top.iterrows():
        lines.append(
            f"- {rec['maintenance_module']}: {rec['week4_priority']}; "
            f"foundation_hits={int(rec['foundation_positive_significant'])}; "
            f"all_validated_p={rec.get('all_validated_p', math.nan):.4g}; "
            f"ncbi_only_p={rec.get('ncbi_only_p', math.nan):.4g}; "
            f"no_protein_rescue_p={rec.get('no_protein_rescue_p', math.nan):.4g}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The strongest Week 4 result is that transposon suppression remains positive and significant in the source-control variants, including `ncbi_only`. This weakens the concern that the Week 3 signal is purely caused by GFF/protein rescue.",
            "",
            "`high_confidence_only` is not a decisive negative control yet because high-confidence scores are nearly all zero under the current NCBI/GFF candidate schema. It should be replaced later by reciprocal best-hit, OMA/OrthoDB/Ensembl agreement, or curated orthology confidence rather than simple source labels.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
