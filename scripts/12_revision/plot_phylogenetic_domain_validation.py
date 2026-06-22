#!/usr/bin/env python3
"""Plot clade, ambiguity, and SAMHD1 domain-aware validation results."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


OUTPUT = Path("results/figures/phylogenetic_domain_validation")
CLADE_COLORS = {
    "Aves": "#C5523C",
    "Mammalia_Chiroptera": "#7B5AA6",
    "Mammalia_nonChiroptera": "#2F718E",
    "Reptilia": "#4F8A5B",
}


def forest(ax, data, y, color, label=None, marker="o"):
    estimate = data["estimate"].to_numpy(float)
    low = data["low"].to_numpy(float)
    high = data["high"].to_numpy(float)
    ax.errorbar(
        estimate,
        y,
        xerr=np.vstack([estimate - low, high - estimate]),
        fmt=marker,
        ms=5,
        capsize=2,
        lw=1.2,
        color=color,
        label=label,
        zorder=3,
    )


def main() -> None:
    clade = pd.read_csv("results/tables/communications_biology_clade_specific_datelife_pgls.tsv", sep="\t")
    ambiguity = pd.read_csv("results/tables/ambiguous_row_scoring_strategy_pgls.tsv", sep="\t")
    samhd1 = pd.read_csv("results/tables/samhd1_domain_robustness_pgls.tsv", sep="\t")
    domain_rows = pd.read_csv("results/tables/targeted_domain_conservation_rows.tsv", sep="\t")
    traits = pd.read_csv(
        "data/processed/maintenance_lifespan_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv", sep="\t"
    )
    traits = traits.loc[traits["score_variant"].eq("phase2_W3_full_background_sensitivity")].drop_duplicates(
        "scientific_name"
    )

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.7), constrained_layout=True)

    ax = axes[0, 0]
    scope_order = [
        "all_primary", "birds_only", "bats_only", "nonflying_mammals_only",
        "mammals_all", "reptiles_only", "no_birds",
    ]
    scope_labels = ["All species", "Birds", "Bats", "Other mammals", "All mammals", "Reptiles", "Without birds"]
    x = clade.loc[
        clade["model"].eq("clade_adjusted_residual") & clade["subset"].isin(scope_order)
    ].set_index("subset").reindex(scope_order)
    plot_data = pd.DataFrame(
        {"estimate": x["estimate_per_score_sd"], "low": x["conf_low"], "high": x["conf_high"]}
    )
    y = np.arange(len(scope_order))[::-1]
    colors = ["#333333", CLADE_COLORS["Aves"], CLADE_COLORS["Mammalia_Chiroptera"], CLADE_COLORS["Mammalia_nonChiroptera"], "#486B8A", CLADE_COLORS["Reptilia"], "#777777"]
    for row_index, color in enumerate(colors):
        forest(ax, plot_data.iloc[[row_index]], np.array([y[row_index]]), color)
    ax.set_yticks(y, scope_labels)
    ax.axvline(0, color="#555555", ls="--", lw=0.9)
    ax.set_xlabel("PGLS effect per score SD")
    ax.set_title("a  Clade-specific lifespan-residual models", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#DDDDDD", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[0, 1]
    strategy_order = ["ambiguous_as_missing", "ambiguous_excluded", "ambiguous_present_like_0.5"]
    strategy_labels = ["Missing", "Excluded", "Present-like (0.5)"]
    amb = ambiguity.loc[
        ambiguity["model"].eq("residual_score") & ambiguity["ambiguity_strategy"].isin(strategy_order)
    ]
    y = np.arange(len(strategy_order))[::-1]
    for label, scope, color, offset in [
        ("All species", "all_species", "#2F718E", 0.12),
        ("Birds only", "birds_only", "#C5523C", -0.12),
    ]:
        table = amb.loc[amb["scope"].eq(scope)].set_index("ambiguity_strategy").reindex(strategy_order)
        plot_data = pd.DataFrame(
            {"estimate": table["estimate_per_score_sd"], "low": table["conf_low"], "high": table["conf_high"]}
        )
        forest(ax, plot_data, y + offset, color, label)
    ax.set_yticks(y, strategy_labels)
    ax.axvline(0, color="#555555", ls="--", lw=0.9)
    ax.set_xlabel("PGLS effect per score SD")
    ax.set_title("b  Ambiguous-row encoding", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2)
    ax.grid(axis="x", color="#DDDDDD", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1, 0]
    sam_scopes = ["all_species", "domain_coverage_ge_0.8", "birds_only", "nonbirds"]
    sam_labels = ["All species", "Domain coverage >=0.8", "Birds", "Non-birds"]
    predictors = [
        ("domain_minus_nondomain_product", "Domain - non-domain", "#C5523C", 0.18),
        ("domain_identity_coverage_product", "Pfam domains", "#2F718E", 0.0),
        ("whole_protein_identity_coverage_product", "Whole protein", "#777777", -0.18),
    ]
    y = np.arange(len(sam_scopes))[::-1]
    for predictor, label, color, offset in predictors:
        table = samhd1.loc[
            samhd1["predictor"].eq(predictor) & samhd1["scope"].isin(sam_scopes)
        ].set_index("scope").reindex(sam_scopes)
        plot_data = pd.DataFrame(
            {"estimate": table["estimate_per_sd"], "low": table["conf_low"], "high": table["conf_high"]}
        )
        forest(ax, plot_data, y + offset, color, label)
    ax.set_yticks(y, sam_labels)
    ax.axvline(0, color="#555555", ls="--", lw=0.9)
    ax.set_xlabel("SAMHD1 PGLS effect per metric SD")
    ax.set_title("c  SAMHD1 domain-specific robustness", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3)
    ax.grid(axis="x", color="#DDDDDD", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1, 1]
    sam_rows = domain_rows.loc[
        domain_rows["human_gene_symbol"].eq("SAMHD1") & domain_rows["sequence_available"].astype(bool)
    ].merge(
        traits[["scientific_name", "pgls_model_c_mass_clade_residual"]],
        on="scientific_name",
        how="left",
    )
    for clade_name, group in sam_rows.groupby("clade"):
        ax.scatter(
            group["domain_minus_nondomain_product"],
            group["pgls_model_c_mass_clade_residual"],
            s=28,
            color=CLADE_COLORS[clade_name],
            alpha=0.82,
            edgecolor="white",
            linewidth=0.4,
            label={
                "Aves": "Birds",
                "Mammalia_Chiroptera": "Bats",
                "Mammalia_nonChiroptera": "Other mammals",
                "Reptilia": "Reptiles",
            }[clade_name],
        )
    ax.axhline(0, color="#999999", lw=0.7)
    ax.set_xlabel("SAMHD1 domain - non-domain conservation")
    ax.set_ylabel("Lifespan residual")
    ax.set_title("d  SAMHD1 species distribution", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="lower right")
    ax.grid(color="#E2E2E2", lw=0.5)
    ax.spines[["top", "right"]].set_visible(False)

    inset = inset_axes(ax, width="52%", height="18%", loc="upper left", borderpad=1.0)
    inset.plot([1, 626], [0, 0], color="#444444", lw=2)
    inset.add_patch(plt.Rectangle((42, -0.18), 66, 0.36, color="#79A6A1"))
    inset.add_patch(plt.Rectangle((164, -0.18), 64, 0.36, color="#C5523C"))
    inset.text(75, 0, "SAM", ha="center", va="center", fontsize=7)
    inset.text(196, 0, "HD", ha="center", va="center", fontsize=7, color="white")
    inset.set_xlim(1, 626)
    inset.set_ylim(-0.35, 0.35)
    inset.set_xticks([1, 626])
    inset.set_yticks([])
    inset.tick_params(axis="x", labelsize=6, length=2)
    inset.spines[["top", "left", "right"]].set_visible(False)

    fig.suptitle(
        "Phylogenetic and domain-aware sensitivity defines the repeat-control claim boundary",
        x=0.01,
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
