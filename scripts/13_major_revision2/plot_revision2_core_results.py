"""Plot study-effort/completeness sensitivity and degradation simulation results."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"


def plot_effort_busco() -> None:
    models = pd.read_csv(TABLES / "jme_revision2_effort_busco_pgls.tsv", sep="\t")
    predictors = ["transposon_repeat_suppression_score", "shared_module_pc1"]
    labels = {
        "base": "Mass + clade",
        "anage_effort": "+ AnAge references",
        "pubmed_effort": "+ PubMed records",
        "dual_effort": "+ both effort proxies",
        "busco_complete_case": "+ BUSCO (complete cases)",
        "busco_and_dual_effort": "+ BUSCO + effort",
        "assembly_and_dual_effort": "+ assembly metrics + effort",
    }
    order = list(labels)
    colors = {predictors[0]: "#236A73", predictors[1]: "#C45B3C"}
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.1), sharey=True)
    for axis, predictor, title in zip(axes, predictors, ["Transposon/repeat score", "Shared six-module axis (PC1)"]):
        subset = models[models.predictor.eq(predictor)].set_index("model").loc[order].reset_index()
        y = np.arange(len(subset))[::-1]
        axis.axvline(0, color="#777777", linewidth=0.8)
        axis.errorbar(
            subset.estimate_per_predictor_sd, y,
            xerr=[subset.estimate_per_predictor_sd - subset.conf_low, subset.conf_high - subset.estimate_per_predictor_sd],
            fmt="o", color=colors[predictor], ecolor=colors[predictor], capsize=3, markersize=5,
        )
        for x, yy, n in zip(subset.estimate_per_predictor_sd, y, subset.n):
            axis.text(x + 0.006, yy + 0.16, f"n={int(n)}", fontsize=7, color="#555555")
        axis.set_yticks(y, [labels[item] for item in order], fontsize=8.5)
        axis.set_title(title, loc="left", fontsize=11, weight="bold")
        axis.set_xlabel("PGLS estimate per predictor SD (95% CI)")
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Sensitivity to phenotype study effort and genome completeness", x=0.08, ha="left", fontsize=14, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / "jme_revision2_effort_busco_forest.png", dpi=500, bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / "jme_revision2_effort_busco_forest.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_degradation() -> None:
    summary = pd.read_csv(TABLES / "jme_revision2_degradation_summary.tsv", sep="\t")
    colors = {"random": "#236A73", "low_research_attention": "#D8A23A", "outcome_related": "#C45B3C"}
    labels = {"random": "Random loss", "low_research_attention": "Low-attention loss", "outcome_related": "Outcome-related loss"}
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.4))
    null = summary[summary.world.eq("null_independent")]
    confounded = summary[summary.world.eq("attention_confounded")]
    planted = summary[summary.world.eq("planted")]
    for mechanism in colors:
        a = null[null.mechanism.eq(mechanism)].sort_values("loss_rate")
        c = confounded[confounded.mechanism.eq(mechanism)].sort_values("loss_rate")
        b = planted[planted.mechanism.eq(mechanism)].sort_values("loss_rate")
        axes[0].plot(100 * a.loss_rate, a["rejected_p_lt_0.05"], marker="o", color=colors[mechanism], label=labels[mechanism])
        axes[1].plot(100 * c.loss_rate, c["rejected_p_lt_0.05"], marker="o", color=colors[mechanism])
        axes[2].plot(100 * b.loss_rate, b.estimated_beta, marker="o", color=colors[mechanism])
    axes[0].axhline(0.05, color="#777777", linestyle="--", linewidth=0.9)
    axes[0].set_title("Independent null", loc="left", weight="bold")
    axes[0].set_ylabel("False positive rate")
    axes[1].axhline(0.05, color="#777777", linestyle="--", linewidth=0.9)
    adjusted = confounded[confounded.mechanism.eq("low_research_attention")].sort_values("loss_rate")
    axes[1].plot(100 * adjusted.loss_rate, adjusted["adjusted_rejected_p_lt_0.05"], marker="s",
                 color="#222222", linestyle="--", label="Low-attention + effort adjustment")
    axes[1].set_title("Attention-confounded null", loc="left", weight="bold")
    axes[1].set_ylabel("False positive rate")
    axes[2].axhline(0.10, color="#777777", linestyle="--", linewidth=0.9)
    axes[2].set_title("Planted biological effect", loc="left", weight="bold")
    axes[2].set_ylabel("Mean estimated effect")
    for axis in axes:
        axis.set_xlabel("Masked observed rows (%)")
        axis.set_ylim(bottom=min(-0.02, axis.get_ylim()[0]))
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    axes[1].legend(frameon=False, fontsize=7.5, loc="lower right")
    fig.suptitle("Annotation-matrix degradation changes comparative inference", x=0.07, ha="left", fontsize=14, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(FIGURES / "jme_revision2_degradation_simulation.png", dpi=500, bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / "jme_revision2_degradation_simulation.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    plot_effort_busco()
    plot_degradation()
    print("Wrote revision 2 core figures.")
