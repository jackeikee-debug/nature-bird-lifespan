"""Build revised main figures for the second JME major revision."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

INK = "#17212B"
MUTED = "#5D6B78"
TEAL = "#236A73"
ORANGE = "#C45B3C"
GOLD = "#D8A23A"
BLUE = "#3478A5"
LIGHT = "#F2F5F6"


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.png", dpi=500, bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def box(ax, x, y, w, h, label, title, lines, color, title_fontsize=10.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.018",
                                facecolor="white", edgecolor=color, linewidth=1.8))
    ax.text(x + 0.025, y + h - 0.04, label, color=color, fontsize=8.5, weight="bold", va="top")
    ax.text(x + 0.025, y + h - 0.09, title, color=INK, fontsize=title_fontsize, weight="bold", va="top")
    for i, line in enumerate(lines):
        ax.text(x + 0.025, y + h - 0.145 - i * 0.033, line, color=MUTED, fontsize=7.3, va="top")


def workflow() -> None:
    fig, ax = plt.subplots(figsize=(11.2, 7.0))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.025, 0.97, "Testing observation bias in comparative genomics of lifespan",
            fontsize=17, weight="bold", color=INK, va="top")
    ax.text(0.025, 0.925, "The design separates biological association, phenotype study effort, genome completeness, and orthology support",
            fontsize=9.8, color=MUTED, va="top")
    box(ax, 0.035, 0.565, 0.28, 0.285, "DATA", "Comparative case study",
        ["68 amniotes: 31 birds, 11 bats,", "18 other mammals, 8 reptiles", "Maximum lifespan, mass, clade", "200 genes in six modules"], BLUE)
    box(ax, 0.36, 0.565, 0.28, 0.285, "HYPOTHESES", "Sources of error",
        ["H1: observation loss lowers scores", "H2: informative loss biases effects", "H3: orthology errors remain after", "assembly-level quality control"], ORANGE)
    box(ax, 0.685, 0.565, 0.28, 0.285, "EMPIRICAL TESTS", "Measured covariates",
        ["BUSCO completeness (49 species)", "Assembly tier and contig N50", "AnAge reference count", "Exact-name PubMed record count"], TEAL)
    box(ax, 0.18, 0.225, 0.28, 0.265, "INTERVENTION", "Annotation-matrix degradation",
        ["31 high-quality genomes", "5%, 10%, 20%, 30% masking", "Random, low-attention, and", "outcome-related observation loss"], GOLD,
        title_fontsize=8.7)
    box(ax, 0.54, 0.225, 0.28, 0.265, "ESTIMANDS", "Quantified performance",
        ["Score RMSE and effect bias", "95% interval coverage", "False positive rate", "Power for a planted effect"], TEAL)
    for start, end in [((.315,.705),(.36,.705)),((.64,.705),(.685,.705)),((.825,.57),(.72,.485)),((.46,.36),(.54,.36))]:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=13, color="#87939D", linewidth=1.4))
    ax.add_patch(FancyBboxPatch((0.035, 0.055), 0.93, 0.11, boxstyle="round,pad=0.012,rounding_size=0.018",
                                facecolor=LIGHT, edgecolor="#C7D0D6", linewidth=1.0))
    ax.text(0.06, 0.125, "SECONDARY SEQUENCE FOLLOW-UP", color=ORANGE, fontsize=9.2, weight="bold", va="center")
    ax.text(0.385, 0.125, "140 low-observability avian rows classified by reciprocal sequence support; SAMHD1 retained as exploratory",
            color=INK, fontsize=7.8, va="center")
    ax.text(0.385, 0.085, "No claim of a causal longevity pathway, bird-specific mechanism, flight convergence, or positive selection",
            color=MUTED, fontsize=7.8, va="center")
    save(fig, "jme_revision2_figure1_study_design")


def shared_axis() -> None:
    all_models = pd.read_csv(TABLES / "jme_revision2_effort_busco_pgls.tsv", sep="\t")
    corr = pd.read_csv(TABLES / "jme_revision2_module_correlations.tsv", sep="\t")
    modules = list(dict.fromkeys(corr.module_1))
    short = {
        "DNA_repair_replication_stress": "DNA repair",
        "cancer_surveillance_senescence": "Cancer control",
        "chromatin_repression_heterochromatin": "Chromatin",
        "inflammation_innate_immune_restraint": "Inflammation",
        "proteostasis_autophagy_mitophagy": "Proteostasis",
        "transposon_repeat_suppression": "Repeat control",
    }
    matrix = corr.pivot(index="module_1", columns="module_2", values="spearman_rho").loc[modules, modules]
    effects = all_models[
        (all_models["model"] == "base")
        & (all_models["predictor"] != "shared_module_pc1")
    ].copy()
    effects["maintenance_module"] = effects["predictor"].str.replace("_score$", "", regex=True)
    effects = effects.set_index("maintenance_module").loc[modules].reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.9), gridspec_kw={"width_ratios": [1.05, 1]})
    y = np.arange(len(effects))[::-1]
    axes[0].axvline(0, color="#777", linewidth=.8)
    axes[0].errorbar(effects.estimate_per_predictor_sd, y,
                     xerr=[effects.estimate_per_predictor_sd-effects.conf_low, effects.conf_high-effects.estimate_per_predictor_sd],
                     fmt="o", color=TEAL, ecolor=TEAL, capsize=3)
    axes[0].set_yticks(y, [short[m] for m in effects.maintenance_module])
    axes[0].set_xlabel("PGLS estimate per score SD (95% CI)")
    axes[0].set_title("a  Concordant module associations", loc="left", weight="bold")
    axes[0].spines[["top","right"]].set_visible(False)
    im = axes[1].imshow(matrix.values, vmin=.7, vmax=1, cmap="YlGnBu")
    axes[1].set_xticks(range(len(modules)), [short[m] for m in modules], rotation=45, ha="right", fontsize=8)
    axes[1].set_yticks(range(len(modules)), [short[m] for m in modules], fontsize=8)
    axes[1].set_title("b  Pairwise Spearman correlations", loc="left", weight="bold")
    for i in range(len(modules)):
        for j in range(len(modules)):
            axes[1].text(j, i, f"{matrix.iloc[i,j]:.2f}", ha="center", va="center",
                         fontsize=7.5, color="white" if matrix.iloc[i,j] > .87 else INK)
    cbar = fig.colorbar(im, ax=axes[1], fraction=.046, pad=.04)
    cbar.set_label("Spearman rho")
    fig.suptitle("Six maintenance modules primarily represent one shared score axis", x=.06, ha="left", fontsize=14, weight="bold")
    fig.text(.06, .015, "Pairwise rho = 0.742-0.923; principal component 1 explains 99.01% of standardized module-score variance.", fontsize=9, color=MUTED)
    fig.tight_layout(rect=[0, .05, 1, .92])
    save(fig, "jme_revision2_figure2_shared_module_axis")


def source_tables() -> None:
    pd.DataFrame([
        ("species_total", 68), ("birds", 31), ("bats", 11), ("other_mammals", 18),
        ("reptiles", 8), ("genes", 200), ("modules", 6), ("busco_available", 49),
        ("high_quality_simulation_species", 31), ("sequence_followup_rows", 140),
    ], columns=["quantity", "value"]).to_csv(TABLES / "jme_revision2_figure1_design_counts.tsv", sep="\t", index=False)

    all_models = pd.read_csv(TABLES / "jme_revision2_effort_busco_pgls.tsv", sep="\t")
    base = all_models[
        (all_models["model"] == "base")
        & (all_models["predictor"] != "shared_module_pc1")
    ].copy()
    base["maintenance_module"] = base["predictor"].str.replace("_score$", "", regex=True)
    base.to_csv(TABLES / "jme_revision2_base_module_models.tsv", sep="\t", index=False)

    legacy = pd.read_csv(TABLES / "jme_table1_gene_matrix_summary.tsv", sep="\t")
    label_to_key = {
        "DNA repair / replication stress": "DNA_repair_replication_stress",
        "Cancer surveillance / senescence": "cancer_surveillance_senescence",
        "Chromatin repression / heterochromatin": "chromatin_repression_heterochromatin",
        "Inflammation / immune restraint": "inflammation_innate_immune_restraint",
        "Proteostasis / autophagy / mitophagy": "proteostasis_autophagy_mitophagy",
        "Transposon / repeat suppression": "transposon_repeat_suppression",
    }
    legacy["maintenance_module"] = legacy["module"].map(label_to_key)
    stats = base[["maintenance_module", "estimate_per_predictor_sd", "conf_low", "conf_high", "p", "q_within_model"]]
    table1 = legacy.drop(columns=["estimate", "conf_low", "conf_high", "p", "p_bh_by_model"]).merge(
        stats, on="maintenance_module", how="left", validate="one_to_one"
    )
    table1 = table1.rename(columns={"estimate_per_predictor_sd": "estimate", "q_within_model": "p_bh_by_model"})
    table1.to_csv(TABLES / "jme_revision2_table1_gene_matrix_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    workflow(); shared_axis(); source_tables()
    print("Wrote revised JME main figures 1 and 2.")
