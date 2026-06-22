"""Plot Phase 3 combined rescue model impact and coverage boundary."""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.pyplot as plt
import pandas as pd


def read(path: pathlib.Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t").fillna("")


def to_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def p_label(value: float) -> str:
    if pd.isna(value):
        return "P=NA"
    if value < 0.001:
        return f"P={value:.1e}"
    return f"P={value:.3g}"


def write_report(
    output: pathlib.Path,
    impact: pd.DataFrame,
    highcov: pd.DataFrame,
    bird: pd.DataFrame,
    png: pathlib.Path,
    pdf: pathlib.Path,
) -> None:
    residual = impact[
        (impact["comparison_layer"] == "all_module_pgls")
        & (impact["model"] == "residual_module")
        & (impact["metric"] == "transposon_rank_estimate_p")
    ].copy()
    baseline = residual[residual["dataset"] == "baseline"].iloc[0]
    rescued = residual[residual["dataset"] == "rescued"].iloc[0]

    highcov_resid = highcov[
        (highcov["model"] == "pgls_clade_residual_transposon")
        & (highcov["term"] == "transposon_repeat_suppression_score")
    ].copy()
    bird_interaction = bird[
        (bird["model"] == "residual_bird_interaction")
        & (bird["term"] == "transposon_repeat_suppression_score:bird_statusbird")
    ].copy()
    bird_p = float(bird_interaction["p"].iloc[0]) if len(bird_interaction) else float("nan")

    lines = [
        "# Phase 3 Model Impact Figure Report",
        "",
        "## Main All-Module Residual Model",
        "",
        f"- Baseline: rank {baseline['rank']}, estimate = {float(baseline['estimate']):.4f}, P = {float(baseline['p']):.4g}.",
        f"- Combined rescued: rank {rescued['rank']}, estimate = {float(rescued['estimate']):.4f}, P = {float(rescued['p']):.4g}.",
        f"- Mean transposon/repeat coverage: {float(baseline['mean_transposon_coverage']):.3f} to {float(rescued['mean_transposon_coverage']):.3f}.",
        "",
        "## Coverage Boundary",
        "",
    ]
    for _, row in highcov_resid.iterrows():
        lines.append(
            f"- {row['subset']}: n = {int(row['n'])}, birds = {int(row['birds'])}, "
            f"estimate = {float(row['estimate']):.4f}, P = {float(row['p']):.4g}."
        )
    lines.extend(
        [
            "",
            "## Bird Interaction",
            "",
            f"- Residual-model bird interaction P = {bird_p:.4g}.",
            "",
            "## Interpretation",
            "",
            "The combined rescue preserves the positive transposon/repeat residual association, but the signal remains weak in high-coverage and birds-only high-coverage subsets. This figure should be used as a claim-boundary plot rather than as evidence for a robust bird-specific mechanism.",
            "",
            "## Outputs",
            f"- PNG: `{png}`",
            f"- PDF: `{pdf}`",
        ]
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot(
    impact: pd.DataFrame,
    highcov: pd.DataFrame,
    bird: pd.DataFrame,
    output_png: pathlib.Path,
    output_pdf: pathlib.Path,
) -> None:
    impact = impact.copy()
    highcov = highcov.copy()
    bird = bird.copy()

    for col in ["estimate", "p", "mean_transposon_coverage"]:
        if col in impact.columns:
            impact[col] = to_num(impact[col])
    for col in ["estimate", "se", "p", "n", "birds"]:
        if col in highcov.columns:
            highcov[col] = to_num(highcov[col])
    for col in ["estimate", "se", "p"]:
        if col in bird.columns:
            bird[col] = to_num(bird[col])

    main = impact[
        (impact["comparison_layer"] == "all_module_pgls")
        & (impact["metric"] == "transposon_rank_estimate_p")
        & (impact["model"].isin(["residual_module", "mass_clade_module"]))
    ].copy()
    main["model_label"] = main["model"].map(
        {
            "residual_module": "Residual model",
            "mass_clade_module": "Mass + clade model",
        }
    )
    main["dataset_label"] = main["dataset"].map({"baseline": "Baseline", "rescued": "Validation overlay"})

    subset_order = [
        "all_primary",
        "transposon_coverage_ge_0_25",
        "transposon_coverage_ge_0_50",
        "transposon_coverage_ge_0_70",
        "all_module_coverage_ge_0_50",
        "birds_only_transposon_coverage_ge_0_50",
    ]
    subset_labels = {
        "all_primary": "All primary",
        "transposon_coverage_ge_0_25": "Transposon coverage >= 0.25",
        "transposon_coverage_ge_0_50": "Transposon coverage >= 0.50",
        "transposon_coverage_ge_0_70": "Transposon coverage >= 0.70",
        "all_module_coverage_ge_0_50": "All-module coverage >= 0.50",
        "birds_only_transposon_coverage_ge_0_50": "Birds only, coverage >= 0.50",
    }
    forest = highcov[
        (highcov["model"] == "pgls_clade_residual_transposon")
        & (highcov["term"] == "transposon_repeat_suppression_score")
        & (highcov["subset"].isin(subset_order))
    ].copy()
    forest["subset"] = pd.Categorical(forest["subset"], categories=subset_order, ordered=True)
    forest = forest.sort_values("subset")
    forest["label"] = forest["subset"].map(subset_labels)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), gridspec_kw={"width_ratios": [1, 1.45]})
    ax0, ax1 = axes

    colors = {"Baseline": "#7A869A", "Validation overlay": "#2E7D5B"}
    x_positions = {"Residual model": 0, "Mass + clade model": 1}
    offsets = {"Baseline": -0.14, "Validation overlay": 0.14}
    for _, row in main.iterrows():
        x = x_positions[row["model_label"]] + offsets[row["dataset_label"]]
        ax0.scatter(x, row["estimate"], s=86, color=colors[row["dataset_label"]], zorder=3)
        text_dx = -0.03 if row["dataset_label"] == "Baseline" else 0.03
        text_ha = "right" if row["dataset_label"] == "Baseline" else "left"
        ax0.text(
            x + text_dx,
            row["estimate"] + 0.018,
            p_label(row["p"]),
            ha=text_ha,
            va="bottom",
            fontsize=8,
            rotation=0,
        )
    for dataset, color in colors.items():
        ax0.scatter([], [], s=86, color=color, label=dataset)
    ax0.set_xticks([0, 1])
    ax0.set_xticklabels(["Residual\nmodel", "Mass + clade\nmodel"])
    ax0.set_ylabel("Transposon/repeat score estimate")
    ax0.set_title("A. Main model after validation")
    ax0.set_xlim(-0.65, 1.42)
    ax0.set_ylim(0, max(0.45, float(main["estimate"].max()) + 0.08))
    ax0.legend(frameon=False, loc="upper left")
    ax0.grid(axis="y", color="#D8DEE9", linewidth=0.8, alpha=0.8)
    ax0.spines[["top", "right"]].set_visible(False)

    y = list(range(len(forest)))[::-1]
    estimates = forest["estimate"].astype(float).to_list()
    ses = forest["se"].astype(float).to_list()
    lower = [est - 1.96 * se for est, se in zip(estimates, ses, strict=False)]
    upper = [est + 1.96 * se for est, se in zip(estimates, ses, strict=False)]
    point_colors = ["#2E7D5B" if p < 0.05 else "#BF616A" for p in forest["p"].astype(float)]
    ax1.axvline(0, color="#4C566A", linewidth=1)
    for yi, est, lo, hi, pval, color, n, birds in zip(
        y,
        estimates,
        lower,
        upper,
        forest["p"].astype(float),
        point_colors,
        forest["n"].astype(int),
        forest["birds"].astype(int),
        strict=False,
    ):
        ax1.plot([lo, hi], [yi, yi], color=color, linewidth=2)
        ax1.scatter(est, yi, s=72, color=color, zorder=3)
        ax1.text(hi + 0.05, yi, f"{p_label(pval)}; n={n}, birds={birds}", va="center", fontsize=8)
    ax1.set_yticks(y)
    ax1.set_yticklabels(forest["label"])
    ax1.set_xlabel("Estimate with approximate 95% CI")
    ax1.set_title("B. Coverage-filtered residual models")
    ax1.grid(axis="x", color="#D8DEE9", linewidth=0.8, alpha=0.8)
    ax1.spines[["top", "right"]].set_visible(False)

    interaction = bird[
        (bird["model"] == "residual_bird_interaction")
        & (bird["term"] == "transposon_repeat_suppression_score:bird_statusbird")
    ]
    if len(interaction):
        pval = float(interaction["p"].iloc[0])
        ax1.text(
            0.02,
            -0.18,
            f"Formal bird interaction: {p_label(pval)}",
            transform=ax1.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            color="#4C566A",
        )

    fig.suptitle("Sequence validation and coverage sensitivity", y=1.02, fontsize=14)
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--impact", type=pathlib.Path, required=True)
    parser.add_argument("--high-coverage", type=pathlib.Path, required=True)
    parser.add_argument("--bird-models", type=pathlib.Path, required=True)
    parser.add_argument("--output-png", type=pathlib.Path, required=True)
    parser.add_argument("--output-pdf", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    impact = read(args.impact)
    highcov = read(args.high_coverage)
    bird = read(args.bird_models)
    plot(impact, highcov, bird, args.output_png, args.output_pdf)
    write_report(args.report, impact, highcov, bird, args.output_png, args.output_pdf)


if __name__ == "__main__":
    main()
