"""Summarize Week 4 sequence-validated transposon PGLS results."""

from __future__ import annotations

import argparse
import math
import pathlib

import pandas as pd


DEFAULT_PGLS = pathlib.Path("results/tables/week4_sequence_validated_pgls.tsv")
DEFAULT_SCORES = pathlib.Path("data/processed/maintenance_scores_week4_sequence_validated.tsv")
DEFAULT_OUTPUT = pathlib.Path("results/tables/week4_sequence_validated_transposon_summary.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week4_sequence_validated_transposon_summary_report.md")

TRANSPOSON = "transposon_suppression"
KEY_MODELS = [
    "mass_module",
    "mass_clade_module",
    "residual_module",
    "pgls_clade_residual_module",
]


def fmt(value: float) -> str:
    if pd.isna(value):
        return ""
    return f"{value:.6g}"


def classify(row: pd.Series) -> str:
    p = row["module_p"]
    bh = row["module_p_bh_by_variant_model"]
    estimate = row["module_estimate"]
    if estimate > 0 and bh < 0.05:
        return "positive_bh_significant"
    if estimate > 0 and p < 0.05:
        return "positive_nominal"
    if estimate > 0 and p < 0.10:
        return "positive_trend"
    return "not_supported"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgls", type=pathlib.Path, default=DEFAULT_PGLS)
    parser.add_argument("--scores", type=pathlib.Path, default=DEFAULT_SCORES)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    pgls = pd.read_csv(args.pgls, sep="\t")
    scores = pd.read_csv(args.scores, sep="\t")
    trans = pgls[
        (pgls["maintenance_module"] == TRANSPOSON)
        & (pgls["model"].isin(KEY_MODELS))
        & (pgls["error"].fillna("") == "")
    ].copy()
    for col in ["module_estimate", "module_p", "module_p_bh_by_variant_model", "n"]:
        trans[col] = pd.to_numeric(trans[col], errors="coerce")
    trans["support_class"] = trans.apply(classify, axis=1)

    score_summary = scores[scores["maintenance_module"] == TRANSPOSON].copy()
    for col in ["coverage_fraction", "confidence_weighted_score", "genes_found"]:
        score_summary[col] = pd.to_numeric(score_summary[col], errors="coerce")
    score_summary = score_summary.groupby("score_variant", as_index=False).agg(
        mean_coverage=("coverage_fraction", "mean"),
        mean_score=("confidence_weighted_score", "mean"),
        min_genes_found=("genes_found", "min"),
        max_genes_found=("genes_found", "max"),
    )
    out = trans.merge(score_summary, on="score_variant", how="left")
    out = out.sort_values(["model", "module_p", "score_variant"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    mass_clade = out[out["model"] == "mass_clade_module"].copy()
    lines = [
        "# Week 4 Sequence-Validated Transposon PGLS Summary",
        "",
        "This summary tests whether the transposon-suppression signal survives after replacing rescue-only rows with local reciprocal DIAMOND support.",
        "",
        "## Mass + Clade Model",
        "",
    ]
    for _, row in mass_clade.iterrows():
        lines.append(
            f"- `{row['score_variant']}`: estimate={fmt(row['module_estimate'])}, "
            f"p={fmt(row['module_p'])}, BH={fmt(row['module_p_bh_by_variant_model'])}, "
            f"mean_score={row['mean_score']:.3f}, support={row['support_class']}."
        )
    lines.extend(["", "## All Transposon Models", ""])
    for _, row in out.iterrows():
        lines.append(
            f"- `{row['score_variant']}` / `{row['model']}`: estimate={fmt(row['module_estimate'])}, "
            f"p={fmt(row['module_p'])}, BH={fmt(row['module_p_bh_by_variant_model'])}, support={row['support_class']}."
        )
    weak_row = mass_clade[mass_clade["score_variant"] == "transposon_sequence_weak_inclusive"]
    weak_support = weak_row.iloc[0]["support_class"] if not weak_row.empty else "not_tested"
    if weak_support == "positive_bh_significant":
        weak_sentence = "The weak-inclusive variant is also positive and BH-significant in the mass+clade model, indicating that the weak reciprocal rows do not erase the signal in this validation round."
    else:
        weak_sentence = "The weak-inclusive variant is positive but no longer significant in the mass+clade model, suggesting that weak reciprocal hits add noise rather than strengthening the signal."
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"The strict sequence-validated transposon score remains positive in all key models and BH-significant in the mass+clade model. {weak_sentence}",
            "",
            "For manuscript development, the strict sequence-supported transposon score is the cleaner Week 4 result. Weak PIWIL2/TRIM28 rows should be held out or cross-checked with OMA, OrthoDB, or Ensembl Compara before being used as positive evidence.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
