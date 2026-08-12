#!/usr/bin/env python
"""Build the combined annotation-degradation and sequence-classification figure."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

INK = "#24343B"
MUTED = "#6C7A81"
GRID = "#DDE3E6"
MECHANISM_COLORS = {
    "random": "#236A73",
    "low_research_attention": "#D8A23A",
    "outcome_related": "#C45B3C",
}
MECHANISM_LABELS = {
    "random": "Random loss",
    "low_research_attention": "Low-attention loss",
    "outcome_related": "Outcome-related loss",
}


def style_axis(axis: plt.Axes) -> None:
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", color=GRID, linewidth=0.75, zorder=0)
    axis.tick_params(colors=INK, labelsize=8)


def add_panel_label(axis: plt.Axes, letter: str, title: str) -> None:
    axis.set_title(f"{letter}  {title}", loc="left", fontsize=10.2, weight="bold", color=INK, pad=8)


def build_figure() -> None:
    summary = pd.read_csv(TABLES / "jme_revision2_degradation_summary.tsv", sep="\t")
    waterfall = pd.read_csv(TABLES / "phase3_evidence_waterfall_counts.tsv", sep="\t")
    by_gene = pd.read_csv(TABLES / "phase3_evidence_ladder_by_gene.tsv", sep="\t")

    fig = plt.figure(figsize=(13.2, 9.2), facecolor="white")
    grid = fig.add_gridspec(2, 6, height_ratios=[0.96, 1.08], hspace=0.48, wspace=0.7)
    axes_top = [fig.add_subplot(grid[0, 0:2]), fig.add_subplot(grid[0, 2:4]), fig.add_subplot(grid[0, 4:6])]
    ax_waterfall = fig.add_subplot(grid[1, 0:2])
    ax_gene = fig.add_subplot(grid[1, 2:6])

    worlds = ["null_independent", "attention_confounded", "planted"]
    titles = ["Independent null", "Attention-confounded null", "Planted biological effect"]
    letters = ["a", "b", "c"]
    for axis, world, title, letter in zip(axes_top, worlds, titles, letters):
        subset = summary[summary["world"].eq(world)]
        for mechanism, color in MECHANISM_COLORS.items():
            rows = subset[subset["mechanism"].eq(mechanism)].sort_values("loss_rate")
            metric = "estimated_beta" if world == "planted" else "rejected_p_lt_0.05"
            axis.plot(
                rows["loss_rate"] * 100,
                rows[metric],
                marker="o",
                markersize=4.2,
                linewidth=1.7,
                color=color,
                label=MECHANISM_LABELS[mechanism],
                zorder=3,
            )
        if world == "attention_confounded":
            rows = subset[subset["mechanism"].eq("low_research_attention")].sort_values("loss_rate")
            axis.plot(
                rows["loss_rate"] * 100,
                rows["adjusted_rejected_p_lt_0.05"],
                marker="s",
                markersize=3.8,
                linewidth=1.3,
                color="#222222",
                linestyle="--",
                label="Low-attention + effort adjustment",
                zorder=3,
            )
        if world == "planted":
            axis.axhline(0.10, color="#767676", linestyle="--", linewidth=1.0, zorder=1)
            axis.set_ylabel("Mean estimated effect", color=INK)
        else:
            axis.axhline(0.05, color="#767676", linestyle="--", linewidth=1.0, zorder=1)
            axis.set_ylabel("False positive rate", color=INK)
        axis.set_xlabel("Masked observed rows (%)", color=INK)
        axis.set_xticks([5, 10, 20, 30])
        add_panel_label(axis, letter, title)
        style_axis(axis)
    axes_top[0].legend(frameon=False, fontsize=7.2, loc="upper left")
    axes_top[1].legend(frameon=False, fontsize=6.8, loc="center right")

    keep = [
        "input_priority1_rows",
        "gff_annotation_rescue",
        "local_gff_protein_strict",
        "local_cds_translation_strict",
        "external_uniprot_strict",
        "partial_or_family_not_scoreable",
        "not_found_or_unresolved",
    ]
    labels = {
        "input_priority1_rows": "Input rows",
        "gff_annotation_rescue": "GFF annotation",
        "local_gff_protein_strict": "Local GFF protein\n+ reciprocal support",
        "local_cds_translation_strict": "Local CDS translation\n+ reciprocal support",
        "external_uniprot_strict": "External UniProt\n+ reciprocal support",
        "partial_or_family_not_scoreable": "Partial / family\nambiguous",
        "not_found_or_unresolved": "Unresolved",
    }
    checkpoint_colors = ["#8A989F", "#3F7CAC", "#2A7886", "#4B8F77", "#735D9C", "#C49A3A", "#B75B55"]
    wf = waterfall.set_index("waterfall_step").reindex(keep).dropna(subset=["rows"]).reset_index()
    y = np.arange(len(wf))
    ax_waterfall.barh(y, wf["rows"], color=checkpoint_colors[: len(wf)], height=0.67, zorder=3)
    ax_waterfall.set_yticks(y, [labels[value] for value in wf["waterfall_step"]], fontsize=7.6)
    ax_waterfall.invert_yaxis()
    ax_waterfall.set_xlabel("Species-gene rows", color=INK)
    ax_waterfall.set_xlim(0, 158)
    for index, value in enumerate(wf["rows"]):
        ax_waterfall.text(value + 2.2, index, str(int(value)), va="center", fontsize=8, weight="bold", color=INK)
    add_panel_label(ax_waterfall, "d", "Sequence-evidence checkpoints")
    style_axis(ax_waterfall)
    ax_waterfall.spines["left"].set_visible(False)
    ax_waterfall.tick_params(axis="y", length=0)

    status_map = {
        "scoreable_strict": "Strict local",
        "scoreable_sensitivity": "Strict external",
        "not_scoreable_not_absence": "Partial/family ambiguous",
        "not_scoreable_review": "Review",
        "not_scoreable_unknown": "Unresolved",
    }
    by_gene["display_class"] = by_gene["phase3_scoring_status"].map(status_map)
    pivot = by_gene.pivot_table(
        index="human_gene_symbol", columns="display_class", values="rows", aggfunc="sum", fill_value=0
    )
    order = ["Strict local", "Strict external", "Partial/family ambiguous", "Review", "Unresolved"]
    pivot = pivot.reindex(columns=order, fill_value=0)
    stack_colors = ["#2A7886", "#735D9C", "#C49A3A", "#987560", "#B75B55"]
    x = np.arange(len(pivot))
    bottom = np.zeros(len(pivot))
    for class_name, color in zip(order, stack_colors):
        values = pivot[class_name].to_numpy()
        ax_gene.bar(x, values, bottom=bottom, label=class_name, color=color, width=0.72, zorder=3)
        bottom += values
    ax_gene.set_xticks(x, pivot.index, rotation=35, ha="right", fontsize=8)
    ax_gene.set_ylabel("Rows", color=INK)
    ax_gene.set_ylim(0, max(15, bottom.max() + 1))
    add_panel_label(ax_gene, "e", "Final evidence classes by gene")
    style_axis(ax_gene)
    ax_gene.legend(frameon=False, fontsize=7.1, loc="upper center", bbox_to_anchor=(0.5, -0.23), ncol=3)

    fig.suptitle(
        "Annotation loss and orthology uncertainty affect different stages of inference",
        x=0.055,
        y=0.985,
        ha="left",
        fontsize=14.2,
        weight="bold",
        color=INK,
    )
    fig.text(
        0.055,
        0.948,
        "Panels a-c intervene on the gene matrix; panels d-e classify sequence evidence in 140 low-observability avian rows.",
        fontsize=8.8,
        color=MUTED,
    )
    fig.subplots_adjust(left=0.16, right=0.985, top=0.90, bottom=0.16)
    fig.savefig(FIGURES / "jme_revision2_figure5_degradation_sequence.png", dpi=500, bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / "jme_revision2_figure5_degradation_sequence.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    wf.to_csv(TABLES / "jme_revision2_figure5_sequence_checkpoints.tsv", sep="\t", index=False)
    pivot.reset_index().to_csv(TABLES / "jme_revision2_figure5_sequence_classes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    build_figure()
    print("Wrote combined Figure 5 and source-data tables.")
