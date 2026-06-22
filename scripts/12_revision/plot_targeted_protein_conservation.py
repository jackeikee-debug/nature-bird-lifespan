#!/usr/bin/env python3
"""Plot targeted protein-sequence coverage and PGLS sensitivity results."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


ROWS = Path("results/tables/targeted_protein_conservation_rows.tsv")
MODELS = Path("results/tables/targeted_protein_conservation_pgls.tsv")
OUTPUT = Path("results/figures/targeted_protein_conservation_sensitivity")

GENE_ORDER = ["DNMT1", "DNMT3A", "DNMT3B", "HELLS", "MBD2", "MBD3", "MORC3", "SAMHD1", "SETDB2", "UHRF1"]
CLADE_ORDER = ["Aves", "Mammalia_Chiroptera", "Mammalia_nonChiroptera", "Reptilia"]
CLADE_LABELS = ["Birds", "Bats", "Other mammals", "Reptiles"]


def main() -> None:
    rows = pd.read_csv(ROWS, sep="\t")
    models = pd.read_csv(MODELS, sep="\t")
    rows["sequence_available"] = rows["sequence_available"].astype(bool)

    availability = (
        rows.groupby(["human_gene_symbol", "clade"])["sequence_available"]
        .agg(["sum", "count", "mean"])
        .reset_index()
    )
    fraction = availability.pivot(index="human_gene_symbol", columns="clade", values="mean").reindex(
        index=GENE_ORDER, columns=CLADE_ORDER
    )
    counts = availability.pivot(index="human_gene_symbol", columns="clade", values="sum").reindex(
        index=GENE_ORDER, columns=CLADE_ORDER
    )
    totals = availability.pivot(index="human_gene_symbol", columns="clade", values="count").reindex(
        index=GENE_ORDER, columns=CLADE_ORDER
    )

    gene_models = models.loc[models["test_family"].eq("gene_level_primary")].copy()
    all_models = gene_models.loc[gene_models["scope"].eq("all_species")].set_index("human_gene_symbol").reindex(GENE_ORDER)
    bird_models = gene_models.loc[gene_models["scope"].eq("aves_only")].set_index("human_gene_symbol").reindex(GENE_ORDER)
    aggregate = models.loc[
        models["test_family"].eq("module_aggregate")
        & models["predictor"].isin(["aggregate_identity_gene_clade_z", "mean_length_completeness"])
    ].copy()
    aggregate_order = ["all_species_min5", "all_species_min8", "aves_min5", "aves_min8"]
    aggregate["scope"] = pd.Categorical(aggregate["scope"], aggregate_order, ordered=True)
    aggregate = aggregate.sort_values(["scope", "predictor"])

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "pdf.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(11.2, 8.4), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.45], height_ratios=[1.25, 1.0])

    ax_a = fig.add_subplot(grid[0, 0])
    cmap = LinearSegmentedColormap.from_list("coverage", ["#F2F0EA", "#7AA6A1", "#1D5A62"])
    image = ax_a.imshow(fraction.values, vmin=0, vmax=1, cmap=cmap, aspect="auto")
    ax_a.set_xticks(range(len(CLADE_LABELS)), CLADE_LABELS, rotation=30, ha="right")
    ax_a.set_yticks(range(len(GENE_ORDER)), GENE_ORDER)
    ax_a.tick_params(length=0)
    for row_index in range(len(GENE_ORDER)):
        for column_index in range(len(CLADE_ORDER)):
            value = fraction.iloc[row_index, column_index]
            text_color = "white" if value >= 0.72 else "#202020"
            ax_a.text(
                column_index,
                row_index,
                f"{int(counts.iloc[row_index, column_index])}/{int(totals.iloc[row_index, column_index])}",
                ha="center",
                va="center",
                color=text_color,
                fontsize=8,
            )
    colorbar = fig.colorbar(image, ax=ax_a, fraction=0.047, pad=0.03)
    colorbar.set_label("Sequence availability")
    ax_a.set_title("a  Protein sequence coverage", loc="left", fontweight="bold")

    ax_b = fig.add_subplot(grid[:, 1])
    y = np.arange(len(GENE_ORDER))[::-1]
    offsets = {"All species": 0.16, "Birds only": -0.16}
    colors = {"All species": "#2D718E", "Birds only": "#C2573F"}
    for label, table in [("All species", all_models), ("Birds only", bird_models)]:
        estimate = table["estimate_per_sd"].to_numpy(float)
        low = table["conf_low"].to_numpy(float)
        high = table["conf_high"].to_numpy(float)
        ax_b.errorbar(
            estimate,
            y + offsets[label],
            xerr=np.vstack([estimate - low, high - estimate]),
            fmt="o",
            ms=5,
            capsize=2,
            lw=1.2,
            color=colors[label],
            label=label,
            zorder=3,
        )
    ax_b.axvline(0, color="#555555", lw=0.9, ls="--")
    ax_b.set_xlim(-0.18, 0.38)
    ax_b.set_yticks(y, GENE_ORDER)
    ax_b.set_xlabel("PGLS effect on lifespan residual per 1 SD\n(identity x human-reference coverage)")
    ax_b.set_title("b  Gene-level sequence-conservation tests", loc="left", fontweight="bold")
    ax_b.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.075), ncol=2)
    ax_b.grid(axis="x", color="#D9D9D9", lw=0.6)
    ax_b.spines[["top", "right"]].set_visible(False)
    for row_y, gene in zip(y, GENE_ORDER):
        q_all = all_models.loc[gene, "q"]
        q_bird = bird_models.loc[gene, "q"]
        right_edge = max(all_models.loc[gene, "conf_high"], bird_models.loc[gene, "conf_high"])
        ax_b.text(right_edge + 0.012, row_y, f"q={q_all:.2f} / {q_bird:.2f}", va="center", fontsize=7.5, color="#555555")
    ax_b.text(0.99, 1.01, "q: all species / birds only", transform=ax_b.transAxes, ha="right", va="bottom", fontsize=7.5, color="#555555")

    ax_c = fig.add_subplot(grid[1, 0])
    scope_labels = {
        "all_species_min5": "All species, >=5 genes",
        "all_species_min8": "All species, >=8 genes",
        "aves_min5": "Birds, >=5 genes",
        "aves_min8": "Birds, >=8 genes",
    }
    metric_labels = {
        "aggregate_identity_gene_clade_z": "Within-clade conservation",
        "mean_length_completeness": "Length completeness",
    }
    positions = np.arange(len(aggregate))[::-1]
    for pos, (_, row) in zip(positions, aggregate.iterrows()):
        color = "#C2573F" if str(row["scope"]).startswith("aves") else "#2D718E"
        ax_c.errorbar(
            row["estimate_per_sd"],
            pos,
            xerr=[[row["estimate_per_sd"] - row["conf_low"]], [row["conf_high"] - row["estimate_per_sd"]]],
            fmt="o",
            ms=4.5,
            capsize=2,
            color=color,
        )
    labels = [f"{scope_labels[str(row.scope)]}\n{metric_labels[row.predictor]}" for _, row in aggregate.iterrows()]
    ax_c.set_yticks(positions, labels, fontsize=7.5)
    ax_c.axvline(0, color="#555555", lw=0.9, ls="--")
    ax_c.set_xlabel("PGLS effect per 1 SD")
    ax_c.set_title("c  Module-level sensitivity", loc="left", fontweight="bold")
    ax_c.grid(axis="x", color="#D9D9D9", lw=0.6)
    ax_c.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        "Targeted protein conservation does not independently explain lifespan residuals",
        x=0.02,
        ha="left",
        fontsize=13,
        fontweight="bold",
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUTPUT.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {OUTPUT.with_suffix('.png')} and {OUTPUT.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
