"""Scan associations between maintenance module scores and lifespan residuals."""

from __future__ import annotations

import argparse
import csv
import math
import pathlib

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf


DEFAULT_SCORES = pathlib.Path("data/processed/maintenance_scores_primary_wide.tsv")
DEFAULT_TRAITS = pathlib.Path("data/processed/pgls_trait_table.tsv")
DEFAULT_PGLS = pathlib.Path("data/processed/pgls_first_pass_residuals.tsv")
DEFAULT_OUTPUT = pathlib.Path("results/tables/maintenance_lifespan_association_primary.tsv")
DEFAULT_MERGED = pathlib.Path("data/processed/maintenance_lifespan_primary.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/maintenance_lifespan_association_primary_report.md")


OUTCOMES = [
    "lifespan_residual_log10",
    "pgls_model_b_mass_flight_residual",
    "pgls_model_c_mass_clade_residual",
]


def score_columns(df: pd.DataFrame) -> list[str]:
    return [col for col in df.columns if col.endswith("_score")]


def safe_corr(x: pd.Series, y: pd.Series, method: str) -> tuple[float, float]:
    mask = x.notna() & y.notna()
    if mask.sum() < 5 or x[mask].nunique() < 2 or y[mask].nunique() < 2:
        return math.nan, math.nan
    if method == "pearson":
        val = stats.pearsonr(x[mask], y[mask])
    else:
        val = stats.spearmanr(x[mask], y[mask])
    return float(val.statistic), float(val.pvalue)


def fit_ols(df: pd.DataFrame, outcome: str, predictor: str) -> tuple[float, float, float, int]:
    cols = [outcome, predictor, "clade", "flight_status", "log10_body_mass_g"]
    data = df[cols].dropna()
    if len(data) < 10 or data[predictor].nunique() < 2:
        return math.nan, math.nan, math.nan, len(data)
    model = smf.ols(f"{outcome} ~ {predictor} + log10_body_mass_g + C(clade)", data=data).fit()
    return (
        float(model.params.get(predictor, math.nan)),
        float(model.pvalues.get(predictor, math.nan)),
        float(model.rsquared_adj),
        int(model.nobs),
    )


def build_merged(scores: pd.DataFrame, traits: pd.DataFrame, pgls: pd.DataFrame) -> pd.DataFrame:
    keep_traits = [
        "scientific_name",
        "opentree_tip_label",
        "body_mass_g",
        "max_lifespan_years",
        "lifespan_residual_log10",
        "lifespan_residual_ratio",
        "data_quality",
        "sample_size",
        "specimen_origin",
    ]
    merged = scores.merge(traits[keep_traits], on="scientific_name", how="left")
    merged["log10_body_mass_g"] = np.log10(pd.to_numeric(merged["body_mass_g"], errors="coerce"))
    merged = merged.merge(pgls, on="opentree_tip_label", how="left")
    return merged


def association_rows(df: pd.DataFrame) -> list[dict[str, str]]:
    rows = []
    for predictor in score_columns(df):
        module = predictor.removesuffix("_score")
        for outcome in OUTCOMES:
            if outcome not in df.columns:
                continue
            x = pd.to_numeric(df[predictor], errors="coerce")
            y = pd.to_numeric(df[outcome], errors="coerce")
            pearson_r, pearson_p = safe_corr(x, y, "pearson")
            spearman_rho, spearman_p = safe_corr(x, y, "spearman")
            beta, beta_p, adj_r2, n_ols = fit_ols(df.assign(**{predictor: x, outcome: y}), outcome, predictor)
            n_pairwise = int((x.notna() & y.notna()).sum())
            rows.append(
                {
                    "maintenance_module": module,
                    "predictor": predictor,
                    "outcome": outcome,
                    "n_pairwise": str(n_pairwise),
                    "pearson_r": f"{pearson_r:.6f}" if not math.isnan(pearson_r) else "",
                    "pearson_p": f"{pearson_p:.6g}" if not math.isnan(pearson_p) else "",
                    "spearman_rho": f"{spearman_rho:.6f}" if not math.isnan(spearman_rho) else "",
                    "spearman_p": f"{spearman_p:.6g}" if not math.isnan(spearman_p) else "",
                    "ols_n": str(n_ols),
                    "ols_beta": f"{beta:.6f}" if not math.isnan(beta) else "",
                    "ols_p": f"{beta_p:.6g}" if not math.isnan(beta_p) else "",
                    "ols_adj_r2": f"{adj_r2:.6f}" if not math.isnan(adj_r2) else "",
                }
            )
    rows.sort(key=lambda row: float(row["spearman_p"] or "inf"))
    return rows


def write_report(path: pathlib.Path, rows: list[dict[str, str]], merged: pd.DataFrame) -> None:
    top = rows[:10]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Maintenance-Lifespan Association Primary Report",
                "",
                f"Species in merged table: {len(merged)}",
                f"Association tests: {len(rows)}",
                "",
                "## Top Spearman Associations",
                *[
                    f"- {row['maintenance_module']} -> {row['outcome']}: rho={row['spearman_rho']}, p={row['spearman_p']}, n={row['n_pairwise']}"
                    for row in top
                ],
                "",
                "## Interpretation",
                "This is a first-pass non-phylogenetic screen on the primary genome panel. Positive hits are hypotheses for Week 4/5 PGLS or phylogenetic mixed-model testing, not final evidence.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", type=pathlib.Path, default=DEFAULT_SCORES)
    parser.add_argument("--traits", type=pathlib.Path, default=DEFAULT_TRAITS)
    parser.add_argument("--pgls", type=pathlib.Path, default=DEFAULT_PGLS)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--merged-output", type=pathlib.Path, default=DEFAULT_MERGED)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    scores = pd.read_csv(args.scores, sep="\t")
    traits = pd.read_csv(args.traits, sep="\t")
    pgls = pd.read_csv(args.pgls, sep="\t")
    merged = build_merged(scores, traits, pgls)
    rows = association_rows(merged)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    args.merged_output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.merged_output, sep="\t", index=False)
    write_report(args.report, rows, merged)
    print(f"Wrote {args.output}, {args.merged_output}, and {args.report}")


if __name__ == "__main__":
    main()
