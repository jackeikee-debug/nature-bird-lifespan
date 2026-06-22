#!/usr/bin/env python
"""Plot PGLS residual and leave-one-species influence diagnostics."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
INFLUENCE = ROOT / "results/tables/communications_biology_species_influence_pgls.tsv"
DIAGNOSTICS = ROOT / "results/tables/communications_biology_pgls_residual_diagnostics.tsv"
OUT_PNG = ROOT / "results/figures/communications_biology_species_influence_diagnostics.png"
OUT_PDF = ROOT / "results/figures/communications_biology_species_influence_diagnostics.pdf"

COLORS = {
    "Aves": "#247B6B",
    "Mammalia_Chiroptera": "#6B5FA7",
    "Mammalia_nonChiroptera": "#3D6C93",
    "Reptilia": "#B56A35",
}


def main() -> None:
    influence = pd.read_csv(INFLUENCE, sep="\t")
    diagnostics = pd.read_csv(DIAGNOSTICS, sep="\t")
    focal = influence.loc[
        (influence["model"] == "mass_clade_module") & influence["error"].fillna("").eq("")
    ].copy()
    focal["abs_shift"] = focal["standardized_estimate_shift"].abs()
    top = focal.nlargest(12, "abs_shift").sort_values("loo_estimate")

    fig, axes = plt.subplots(1, 2, figsize=(14.2, 6.3), gridspec_kw={"width_ratios": [1.0, 1.12]})
    ax = axes[0]
    for clade, group in diagnostics.groupby("clade"):
        ax.scatter(
            group["fitted_log10_lifespan"], group["standardized_residual"],
            s=38 + 260 * group["cook_distance_gls_approx"].clip(upper=0.25),
            color=COLORS.get(clade, "#666666"), alpha=0.82, edgecolor="white", linewidth=0.5,
            label=clade.replace("Mammalia_", ""),
        )
    ax.axhline(0, color="#777777", lw=1)
    ax.axhline(2, color="#B24745", lw=1, ls="--")
    ax.axhline(-2, color="#B24745", lw=1, ls="--")
    offsets = {
        "Anolis carolinensis": (8, 8),
        "Candoia aspera": (8, -20),
        "Pedionomus torquatus": (8, 8),
        "Saccopteryx bilineata": (-82, 8),
        "Heterocephalus glaber": (8, 5),
    }
    for _, row in diagnostics.nlargest(5, "cook_distance_gls_approx").iterrows():
        offset = offsets.get(row["scientific_name"], (6, 6))
        ax.annotate(
            row["scientific_name"],
            (row["fitted_log10_lifespan"], row["standardized_residual"]),
            xytext=offset, textcoords="offset points", fontsize=7.5,
            arrowprops={"arrowstyle": "-", "color": "#777777", "lw": 0.6},
        )
    ax.set_xlabel("Fitted log10 maximum lifespan")
    ax.set_ylabel("Standardized PGLS residual")
    ax.set_title("A. Residual and leverage screen", loc="left", weight="bold")
    ax.grid(color="#E5E7EB", lw=0.7)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    y = np.arange(len(top))
    baseline = float(top["baseline_estimate"].iloc[0])
    baseline_se = float(top["baseline_se"].iloc[0])
    ax.axvspan(baseline - 1.96 * baseline_se, baseline + 1.96 * baseline_se, color="#DCEBE7", alpha=0.75)
    ax.axvline(baseline, color="#247B6B", lw=1.5, label="Full-panel estimate")
    for pos, (_, row) in zip(y, top.iterrows()):
        color = COLORS.get(row["clade"], "#666666")
        ax.errorbar(
            row["loo_estimate"], pos,
            xerr=[[row["loo_estimate"] - row["loo_conf_low"]], [row["loo_conf_high"] - row["loo_estimate"]]],
            fmt="o", color=color, ecolor=color, capsize=3, elinewidth=1.5, markersize=5.5,
        )
    ax.set_yticks(y, [x.replace("_", " ") for x in top["scientific_name"]])
    ax.set_xlabel("Estimate after omitting one species (95% CI)")
    ax.set_title("B. Twelve largest estimate shifts", loc="left", weight="bold")
    ax.grid(axis="x", color="#E5E7EB", lw=0.7)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0, labelsize=8.5)

    fig.suptitle("Species-level influence diagnostics for the transposon/repeat PGLS model", x=0.07, ha="left", fontsize=15, weight="bold")
    fig.text(0.07, 0.925, "Point size in panel A reflects the GLS Cook-distance approximation; all 68 omissions are reported in Supplementary Data.", fontsize=9, color="#555555")
    fig.subplots_adjust(left=0.07, right=0.98, top=0.86, bottom=0.12, wspace=0.62)
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
