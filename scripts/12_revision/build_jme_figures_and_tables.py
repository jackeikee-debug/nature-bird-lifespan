#!/usr/bin/env python
"""Build Journal of Molecular Evolution figures and summary tables."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

COLORS = {
    "ink": "#18212B",
    "muted": "#607080",
    "blue": "#2F6F8F",
    "teal": "#2A8C82",
    "green": "#4D8B63",
    "orange": "#D06B4F",
    "gold": "#C69A3B",
    "purple": "#7D65A8",
    "rose": "#B95C6B",
    "light": "#EEF2F4",
}

CLADE_COLORS = {
    "Aves": COLORS["orange"],
    "Mammalia_Chiroptera": COLORS["purple"],
    "Mammalia_nonChiroptera": COLORS["blue"],
    "Reptilia": COLORS["green"],
}

CLADE_LABELS = {
    "Aves": "Birds",
    "Mammalia_Chiroptera": "Bats",
    "Mammalia_nonChiroptera": "Other mammals",
    "Reptilia": "Reptiles",
}


def save_figure(fig: plt.Figure, stem: str, dpi: int = 320) -> None:
    fig.savefig(FIGURES / f"{stem}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_summary_tables() -> None:
    reconciliation = pd.read_csv(
        TABLES / "communications_biology_gene_panel_reconciliation_by_module.tsv", sep="\t"
    )
    effects = pd.read_csv(TABLES / "communications_biology_module_effect_forest.tsv", sep="\t")
    effects = effects[["maintenance_module", "mean_coverage", "estimate", "conf_low", "conf_high", "p", "p_bh_by_model"]]
    table1 = reconciliation.merge(
        effects, left_on="maintenance_module_v2", right_on="maintenance_module", how="left"
    )
    labels = {
        "DNA_repair_replication_stress": "DNA repair / replication stress",
        "proteostasis_autophagy_mitophagy": "Proteostasis / autophagy / mitophagy",
        "inflammation_innate_immune_restraint": "Inflammation / immune restraint",
        "cancer_surveillance_senescence": "Cancer surveillance / senescence",
        "chromatin_repression_heterochromatin": "Chromatin repression / heterochromatin",
        "transposon_repeat_suppression": "Transposon / repeat suppression",
    }
    table1.insert(0, "module", table1["maintenance_module_v2"].map(labels))
    table1 = table1[
        [
            "module",
            "expanded_genes",
            "final_genes",
            "removed_genes",
            "domain_required",
            "exclude_from_absence_scoring",
            "mean_coverage",
            "estimate",
            "conf_low",
            "conf_high",
            "p",
            "p_bh_by_model",
        ]
    ]
    table1.to_csv(TABLES / "jme_table1_gene_matrix_summary.tsv", sep="\t", index=False)

    criteria = [
        {
            "evidence_class": "Strict local GFF protein",
            "required_evidence": "Target-species GFF-linked protein and reciprocal same-gene sequence support",
            "orthology_risk_handling": "Reject if a related paralog is the better reciprocal match",
            "scoring_use": "Primary strict presence",
        },
        {
            "evidence_class": "Strict local CDS translation",
            "required_evidence": "Target-species assembly CDS translated from a GFF-linked transcript and reciprocal same-gene support",
            "orthology_risk_handling": "Apply the same paralog and fragment checks as for proteins",
            "scoring_use": "Primary strict presence",
        },
        {
            "evidence_class": "Strict external UniProt",
            "required_evidence": "Target-species UniProt sequence and reciprocal same-gene support",
            "orthology_risk_handling": "Retain outside the local assembly evidence tier",
            "scoring_use": "External sensitivity presence",
        },
        {
            "evidence_class": "Partial or family ambiguous",
            "required_evidence": "Short fragment, conserved-domain-only hit, or unresolved same-family placement",
            "orthology_risk_handling": "Do not promote a family-level match to a gene-level ortholog",
            "scoring_use": "Not scoreable and not evidence of absence",
        },
        {
            "evidence_class": "Review",
            "required_evidence": "Weak, conflicting, or wrong-family evidence",
            "orthology_risk_handling": "Retain the reason for rejection at row level",
            "scoring_use": "Not scoreable",
        },
        {
            "evidence_class": "Unresolved or not found",
            "required_evidence": "No sequence-supported assignment after local and external review",
            "orthology_risk_handling": "Do not infer biological absence from failed observation",
            "scoring_use": "Not scoreable unless an independent absence rule applies",
        },
    ]
    pd.DataFrame(criteria).to_csv(TABLES / "jme_orthology_evidence_criteria.tsv", sep="\t", index=False)

    selection_order = [
        "DNMT1",
        "DNMT3A",
        "DNMT3B",
        "HELLS",
        "UHRF1",
        "SETDB2",
        "MBD2",
        "MBD3",
        "MORC3",
        "SAMHD1",
    ]
    eligibility = pd.read_csv(ROOT / "data/processed/phase2_strict_v2_scoring_eligibility.tsv", sep="\t")
    selection = eligibility[eligibility["human_gene_symbol"].isin(selection_order)].copy()
    selection["rescue_rank"] = selection["human_gene_symbol"].map(
        {gene: rank for rank, gene in enumerate(selection_order, start=1)}
    )
    selection["audited_low_coverage_bird_species"] = 14
    selection["audited_species_gene_rows"] = 14
    selection["pre_rescue_selection_basis"] = (
        "High orthology-validation priority in the transposon/repeat module; selected in a fixed rescue order for "
        "low-observability avian rows with module-high-priority or paralog/family risk"
    )
    selection["posthoc_protein_or_domain_result_used"] = False
    selection = selection.sort_values("rescue_rank")[[
        "rescue_rank",
        "human_gene_symbol",
        "submodule_v2",
        "orthology_validation_priority",
        "gene_family_risk",
        "v2_scoring_group",
        "claim_use",
        "audited_low_coverage_bird_species",
        "audited_species_gene_rows",
        "pre_rescue_selection_basis",
        "posthoc_protein_or_domain_result_used",
    ]]
    selection.to_csv(TABLES / "jme_audited_gene_selection.tsv", sep="\t", index=False)

    interpro_columns = [
        "sequence_id",
        "md5",
        "sequence_length",
        "analysis",
        "signature_accession",
        "signature_description",
        "start",
        "end",
        "score",
        "status",
        "date",
        "interpro_accession",
        "interpro_description",
        "go_terms",
        "pathway_terms",
    ]
    interpro = pd.read_csv(
        TABLES / "human_reference_10genes_interpro_pfam.tsv",
        sep="\t",
        header=None,
        names=interpro_columns,
    )
    interpro.to_csv(
        TABLES / "human_reference_10genes_interpro_pfam_with_header.tsv",
        sep="\t",
        index=False,
    )


def add_box(ax, xy, width, height, title, lines, color, step):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        facecolor="white",
        edgecolor=color,
        linewidth=1.8,
    )
    ax.add_patch(patch)
    ax.add_patch(Rectangle((x, y + height - 0.08), width, 0.08, facecolor=color, edgecolor="none"))
    ax.text(x + 0.025, y + height - 0.04, step, color="white", fontsize=9, weight="bold", va="center")
    ax.text(x + 0.025, y + height - 0.115, title, color=COLORS["ink"], fontsize=10.3, weight="bold", va="top")
    for i, line in enumerate(lines):
        ax.text(x + 0.025, y + height - 0.17 - i * 0.041, line, color=COLORS["muted"], fontsize=8.3, va="top")


def build_workflow_figure() -> None:
    fig, ax = plt.subplots(figsize=(10.6, 7.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(
        0.02,
        0.97,
        "Gene observability in comparative genomics of vertebrate lifespan evolution",
        fontsize=17,
        weight="bold",
        color=COLORS["ink"],
        va="top",
    )
    ax.text(
        0.02,
        0.925,
        "Study design separates trait association, observability, orthology support, and domain evidence",
        fontsize=10.2,
        color=COLORS["muted"],
        va="top",
    )

    boxes = [
        ((0.03, 0.585), "Comparative phenotype", ["Maximum lifespan", "Body mass and clade", "Residual adjusted for mass and clade"], COLORS["blue"], "1  TRAITS"),
        ((0.355, 0.585), "Genome maintenance matrix", ["68 vertebrate species", "236 gene design -> 200 scored", "Six prespecified modules"], COLORS["teal"], "2  MATRIX"),
        ((0.68, 0.585), "Phylogenetic association", ["DateLife calibrated OpenTree", "PGLS and phylogenetic signal", "Branch and clade sensitivity"], COLORS["purple"], "3  MODELS"),
        ((0.03, 0.235), "Gene observability", ["Genome tier and module coverage", "5,000 matched random gene sets", "High coverage sensitivity"], COLORS["gold"], "4  BIAS CONTROL"),
        ((0.355, 0.235), "Orthology evidence audit", ["140 avian rows with low coverage", "GFF/CDS + reciprocal sequence", "Strict, ambiguous, unresolved"], COLORS["orange"], "5  SEQUENCE AUDIT"),
        ((0.68, 0.235), "Protein and domain evolution", ["508 proteins; gene-wise MAFFT", "InterProScan/Pfam projection", "SAMHD1 SAM/HD hypothesis"], COLORS["green"], "6  DOMAIN TEST"),
    ]
    for xy, title, lines, color, step in boxes:
        add_box(ax, xy, 0.29, 0.27, title, lines, color, step)

    arrows = [
        ((0.32, 0.72), (0.355, 0.72)),
        ((0.645, 0.72), (0.68, 0.72)),
        ((0.825, 0.585), (0.825, 0.505)),
        ((0.68, 0.37), (0.645, 0.37)),
        ((0.355, 0.37), (0.32, 0.37)),
        ((0.175, 0.505), (0.175, 0.585)),
    ]
    for start, end in arrows:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=13, color="#89949E", linewidth=1.4))

    ax.add_patch(
        FancyBboxPatch(
            (0.03, 0.035),
            0.94,
            0.135,
            boxstyle="round,pad=0.012,rounding_size=0.018",
            facecolor=COLORS["light"],
            edgecolor="#C8D0D7",
            linewidth=1.0,
        )
    )
    ax.text(0.055, 0.132, "INFERENCE BOUNDARY", fontsize=9.7, weight="bold", color=COLORS["rose"], va="center")
    ax.text(
        0.27,
        0.132,
        "Broad maintenance/annotation axis; repeat-control as a worked stress test; exploratory SAMHD1 domain candidate",
        fontsize=9.3,
        color=COLORS["ink"],
        va="center",
    )
    ax.text(
        0.27,
        0.078,
        "No claim of pathway specificity, bird specific mechanism, flight convergence, positive selection, or altered enzyme activity",
        fontsize=8.9,
        color=COLORS["muted"],
        va="center",
    )
    save_figure(fig, "jme_figure1_study_design")


def build_module_forest() -> None:
    df = pd.read_csv(TABLES / "communications_biology_module_effect_forest.tsv", sep="\t")
    df = df.sort_values("estimate", ascending=True).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9.4, 5.6))
    y = np.arange(len(df))
    colors = [COLORS["orange"] if x == "transposon_repeat_suppression" else COLORS["teal"] for x in df["maintenance_module"]]
    for i, row in df.iterrows():
        ax.errorbar(
            row["estimate"],
            i,
            xerr=[[row["estimate"] - row["conf_low"]], [row["conf_high"] - row["estimate"]]],
            fmt="o",
            color=colors[i],
            markersize=5.5,
            elinewidth=1.8,
            capsize=3,
            zorder=3,
        )
    ax.axvline(0, color="#777777", linestyle="--", linewidth=1)
    ax.set_yticks(y, df["module_label"])
    ax.set_xlabel("Module-score effect on log10 maximum lifespan (95% CI)")
    ax.set_title("Broad genome-maintenance module effects in mass- and clade-adjusted lifespan models", loc="left", weight="bold", pad=13)
    ax.text(
        0,
        1.01,
        "DateLife-tree models of log10 lifespan adjusted for log10 body mass and four-level clade; n = 68",
        transform=ax.transAxes,
        fontsize=9,
        color=COLORS["muted"],
        va="bottom",
    )
    right = max(df["conf_high"].max() + 0.18, 0.72)
    ax.set_xlim(-0.03, right)
    for i, row in df.iterrows():
        ax.text(
            row["conf_high"] + 0.018,
            i,
            f"q={row['p_bh_by_model']:.3g};  n={int(row['genes_total'])} genes",
            va="center",
            fontsize=9.2,
            color=COLORS["muted"],
        )
    ax.grid(axis="x", color="#D9DEE2", linewidth=0.8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    save_figure(fig, "jme_figure2_module_forest", dpi=450)


def build_orthology_ladder() -> None:
    waterfall = pd.read_csv(TABLES / "phase3_evidence_waterfall_counts.tsv", sep="\t")
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
        "partial_or_family_not_scoreable": "Partial / paralog\nambiguous",
        "not_found_or_unresolved": "Unresolved",
    }
    class_colors = [COLORS["muted"], COLORS["blue"], COLORS["teal"], COLORS["green"], COLORS["purple"], COLORS["gold"], COLORS["rose"]]
    wf = waterfall.set_index("waterfall_step").loc[keep].reset_index()

    by_gene = pd.read_csv(TABLES / "phase3_evidence_ladder_by_gene.tsv", sep="\t")
    status_map = {
        "scoreable_strict": "Strict local",
        "scoreable_sensitivity": "Strict external",
        "not_scoreable_not_absence": "Partial/family ambiguous",
        "not_scoreable_review": "Review",
        "not_scoreable_unknown": "Unresolved",
    }
    by_gene["display_class"] = by_gene["phase3_scoring_status"].map(status_map)
    pivot = by_gene.pivot_table(index="human_gene_symbol", columns="display_class", values="rows", aggfunc="sum", fill_value=0)
    order = ["Strict local", "Strict external", "Partial/family ambiguous", "Review", "Unresolved"]
    pivot = pivot.reindex(columns=order, fill_value=0)
    stack_colors = [COLORS["teal"], COLORS["purple"], COLORS["gold"], "#9B7C65", COLORS["rose"]]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 5.7), gridspec_kw={"width_ratios": [1.0, 1.38]})
    y = np.arange(len(wf))
    ax1.barh(y, wf["rows"], color=class_colors, height=0.68)
    ax1.set_yticks(y, [labels[x] for x in wf["waterfall_step"]])
    ax1.invert_yaxis()
    ax1.set_xlabel("Species-gene rows")
    ax1.set_title("a  Orthology-evidence ladder", loc="left", weight="bold")
    for i, value in enumerate(wf["rows"]):
        ax1.text(value + 2, i, str(int(value)), va="center", fontsize=9, weight="bold", color=COLORS["ink"])
    ax1.set_xlim(0, 155)
    ax1.grid(axis="x", color="#E0E4E7", linewidth=0.8)
    ax1.spines[["top", "right", "left"]].set_visible(False)
    ax1.tick_params(axis="y", length=0)

    x = np.arange(len(pivot))
    bottom = np.zeros(len(pivot))
    for cls, color in zip(order, stack_colors):
        vals = pivot[cls].to_numpy()
        ax2.bar(x, vals, bottom=bottom, label=cls, color=color, width=0.72)
        bottom += vals
    ax2.set_xticks(x, pivot.index, rotation=45, ha="right")
    ax2.set_ylabel("Rows")
    ax2.set_title("b  Final evidence classes by gene", loc="left", weight="bold")
    ax2.set_ylim(0, max(15, bottom.max() + 1))
    ax2.grid(axis="y", color="#E0E4E7", linewidth=0.8)
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=3)
    fig.suptitle("Orthology and sequence-evidence audit of low-coverage avian repeat-control rows", x=0.02, ha="left", fontsize=15.5, weight="bold")
    fig.text(0.02, 0.925, "Counts in panel a are evidence-ladder checkpoints and are not mutually exclusive partitions.", fontsize=9, color=COLORS["muted"])
    fig.subplots_adjust(left=0.22, right=0.98, top=0.85, bottom=0.28, wspace=0.35)
    save_figure(fig, "jme_figure3_orthology_ladder")


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            elif name is not None:
                records[name].append(line)
    return {key: "".join(value) for key, value in records.items()}


def build_samhd1_figure() -> None:
    rows = pd.read_csv(TABLES / "targeted_domain_conservation_rows.tsv", sep="\t")
    rows = rows[
        (rows["human_gene_symbol"] == "SAMHD1")
        & rows["sequence_available"].astype(bool)
        & (rows["domain_reference_coverage"] >= 0.5)
    ].copy()
    traits = pd.read_csv(TABLES / "targeted_domain_conservation_species.tsv", sep="\t")
    rows = rows.merge(traits[["scientific_name", "pgls_model_c_mass_clade_residual"]], on="scientific_name", how="left")
    sequences = read_fasta(ROOT / "data" / "interim" / "protein_conservation" / "SAMHD1.aligned.faa")
    reference = sequences["REF_Homo_sapiens"]
    ref_columns = [i for i, aa in enumerate(reference) if aa != "-"]
    if len(ref_columns) != 626:
        raise ValueError(f"Expected 626 human SAMHD1 residues, found {len(ref_columns)}")

    clade_order = {name: i for i, name in enumerate(CLADE_COLORS)}
    rows["clade_order"] = rows["clade"].map(clade_order)
    rows = rows.sort_values(["clade_order", "domain_minus_nondomain_product", "scientific_name"]).reset_index(drop=True)

    matrix = []
    for _, row in rows.iterrows():
        seq = sequences[row["alignment_record_id"]]
        states = []
        for col in ref_columns:
            aa = seq[col]
            if aa == "-":
                states.append(0)
            elif aa == reference[col]:
                states.append(2)
            else:
                states.append(1)
        matrix.append(states)
    matrix = np.asarray(matrix)

    robustness = pd.read_csv(TABLES / "samhd1_domain_robustness_pgls.tsv", sep="\t")
    scopes = ["all_species", "domain_coverage_ge_0.8", "birds_only", "nonbirds"]
    predictors = ["domain_minus_nondomain_product", "domain_identity_coverage_product", "whole_protein_identity_coverage_product"]
    forest = robustness[robustness["scope"].isin(scopes) & robustness["predictor"].isin(predictors)].copy()
    scope_labels = {
        "all_species": "All species",
        "domain_coverage_ge_0.8": "Domain coverage >= 0.8",
        "birds_only": "Birds",
        "nonbirds": "Non-birds",
    }
    predictor_labels = {
        "domain_minus_nondomain_product": "Domain - non-domain",
        "domain_identity_coverage_product": "Pfam domains",
        "whole_protein_identity_coverage_product": "Whole protein",
    }
    predictor_colors = {
        "domain_minus_nondomain_product": COLORS["orange"],
        "domain_identity_coverage_product": COLORS["blue"],
        "whole_protein_identity_coverage_product": "#7B7B7B",
    }

    fig = plt.figure(figsize=(11.4, 9.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.16, 1.15, 1.0], hspace=0.34, wspace=0.28)
    ax_arch = fig.add_subplot(gs[0, :])
    ax_heat = fig.add_subplot(gs[1, :])
    ax_forest = fig.add_subplot(gs[2, 0])
    ax_scatter = fig.add_subplot(gs[2, 1])

    ax_arch.set_xlim(1, 626)
    ax_arch.set_ylim(0, 1)
    ax_arch.axis("off")
    ax_arch.plot([1, 626], [0.38, 0.38], color=COLORS["ink"], linewidth=3)
    for start, end, label, color in [(42, 107, "SAM", "#6AA49D"), (164, 227, "HD", COLORS["orange"])]:
        ax_arch.add_patch(Rectangle((start, 0.13), end - start + 1, 0.5, facecolor=color, edgecolor="none"))
        ax_arch.text((start + end) / 2, 0.38, label, ha="center", va="center", fontsize=9, weight="bold", color="white")
        ax_arch.text((start + end) / 2, 0.72, f"{start}-{end}", ha="center", va="bottom", fontsize=8, color=COLORS["muted"])
    ax_arch.text(1, 0.93, "a  Human SAMHD1 domain architecture (NP_056289.2; 626 aa)", fontsize=10.5, weight="bold", va="bottom")

    cmap = ListedColormap(["#E7EAED", "#E6B8A9", "#277D78"])
    ax_heat.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=2, extent=[1, 626, len(rows), 0])
    ax_heat.axvspan(42, 107, facecolor="#6AA49D", alpha=0.13, edgecolor="#4B837D", linewidth=1)
    ax_heat.axvspan(164, 227, facecolor=COLORS["orange"], alpha=0.12, edgecolor=COLORS["orange"], linewidth=1)
    ax_heat.set_xlabel("Human-reference amino-acid position")
    ax_heat.set_ylabel("")
    ax_heat.set_title("b  Residue identity projected onto human SAMHD1 coordinates", loc="left", fontsize=10.5, weight="bold")
    ax_heat.set_yticks([])
    boundaries = []
    start = 0
    for clade in CLADE_COLORS:
        n = int((rows["clade"] == clade).sum())
        if n:
            ax_heat.add_patch(Rectangle((-13, start), 7, n, transform=ax_heat.transData, facecolor=CLADE_COLORS[clade], clip_on=False))
            ax_heat.text(-18, start + n / 2, f"{CLADE_LABELS[clade]} ({n})", ha="right", va="center", fontsize=8.2, color=COLORS["ink"], clip_on=False)
            start += n
            boundaries.append(start)
    for boundary in boundaries[:-1]:
        ax_heat.axhline(boundary, color="white", linewidth=1.3)
    heat_legend = [
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#277D78", markeredgecolor="none", markersize=8, label="Identical"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#E6B8A9", markeredgecolor="none", markersize=8, label="Substitution"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#E7EAED", markeredgecolor="none", markersize=8, label="Gap"),
    ]
    ax_heat.legend(handles=heat_legend, frameon=False, ncol=3, loc="upper right", bbox_to_anchor=(1, 1.14), fontsize=8)

    offsets = {predictors[0]: -0.18, predictors[1]: 0, predictors[2]: 0.18}
    for i, scope in enumerate(scopes):
        subset = forest[forest["scope"] == scope]
        for predictor in predictors:
            row = subset[subset["predictor"] == predictor]
            if row.empty:
                continue
            row = row.iloc[0]
            y = len(scopes) - 1 - i + offsets[predictor]
            ax_forest.errorbar(
                row["estimate_per_sd"],
                y,
                xerr=[[row["estimate_per_sd"] - row["conf_low"]], [row["conf_high"] - row["estimate_per_sd"]]],
                fmt="o",
                color=predictor_colors[predictor],
                markersize=4.8,
                capsize=2.5,
                linewidth=1.4,
            )
    ax_forest.axvline(0, color="#666666", linestyle="--", linewidth=1)
    ax_forest.set_yticks(range(len(scopes)), [scope_labels[x] for x in reversed(scopes)])
    ax_forest.set_xlabel("PGLS effect per metric SD (95% CI)")
    ax_forest.set_title("c  Domain-specific robustness", loc="left", fontsize=10.5, weight="bold")
    ax_forest.grid(axis="x", color="#E0E4E7", linewidth=0.8)
    ax_forest.spines[["top", "right", "left"]].set_visible(False)
    ax_forest.tick_params(axis="y", length=0)
    legend = [Line2D([0], [0], marker="o", color=predictor_colors[p], label=predictor_labels[p], linewidth=1.4, markersize=5) for p in predictors]
    ax_forest.legend(handles=legend, frameon=False, fontsize=7.8, loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2)

    for clade, color in CLADE_COLORS.items():
        sub = rows[rows["clade"] == clade]
        ax_scatter.scatter(
            sub["domain_minus_nondomain_product"],
            sub["pgls_model_c_mass_clade_residual"],
            s=35,
            color=color,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.45,
            label=CLADE_LABELS[clade],
        )
    ax_scatter.axhline(0, color="#888888", linewidth=0.9)
    ax_scatter.axvline(0, color="#D0D4D7", linewidth=0.8)
    ax_scatter.set_xlabel("SAMHD1 domain - non-domain conservation")
    ax_scatter.set_ylabel("Lifespan residual")
    ax_scatter.set_title("d  Species distribution", loc="left", fontsize=10.5, weight="bold")
    ax_scatter.grid(color="#E0E4E7", linewidth=0.7)
    ax_scatter.spines[["top", "right"]].set_visible(False)
    ax_scatter.legend(frameon=False, fontsize=7.8, loc="lower right", ncol=2)
    ax_scatter.text(
        0.02,
        0.98,
        "DateLife PGLS: beta = 0.128 per SD\nP = 7.30 x 10^-5; joint q = 0.00146",
        transform=ax_scatter.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        color=COLORS["ink"],
        bbox={"facecolor": "white", "edgecolor": "#CDD3D8", "boxstyle": "round,pad=0.3", "alpha": 0.92},
    )

    fig.suptitle("SAMHD1 domain conservation emerges as an exploratory molecular-evolutionary hypothesis", x=0.03, y=0.995, ha="left", fontsize=15.2, weight="bold")
    fig.subplots_adjust(left=0.18, right=0.98, top=0.93, bottom=0.13)
    save_figure(fig, "jme_figure4_samhd1_domain_evolution", dpi=600)


def write_figure_source_data() -> None:
    manifest = [
        ("Figure_1_design_counts.tsv", TABLES / "jme_table1_gene_matrix_summary.tsv"),
        ("Figure_2_module_models.tsv", TABLES / "communications_biology_module_effect_forest.tsv"),
        ("Figure_3_evidence_ladder.tsv", TABLES / "phase3_evidence_waterfall_counts.tsv"),
        ("Figure_3_evidence_by_gene.tsv", TABLES / "phase3_evidence_ladder_by_gene.tsv"),
        ("Figure_4_samhd1_species.tsv", TABLES / "targeted_domain_conservation_rows.tsv"),
        ("Figure_4_samhd1_robustness.tsv", TABLES / "samhd1_domain_robustness_pgls.tsv"),
        ("Figure_4_samhd1_leave_one_species.tsv", TABLES / "samhd1_domain_leave_one_species_out.tsv"),
    ]
    with (TABLES / "jme_figure_source_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["submission_filename", "workspace_source"])
        for target, source in manifest:
            writer.writerow([target, str(source.relative_to(ROOT)).replace("\\", "/")])

    figure_map = [
        ("Figure 1", "Source_Data/Figure_1_design_counts.tsv"),
        ("Figure 2", "Source_Data/Figure_2_module_ranking.tsv; Supplementary_Data/Supplementary_Data_2_final_module_models.tsv"),
        ("Figure 3", "Source_Data/Figure_3_evidence_waterfall.tsv; Source_Data/Figure_3_evidence_by_gene.tsv"),
        ("Figure 4", "Source_Data/Figure_4_samhd1_species.tsv; Source_Data/Figure_4_samhd1_robustness.tsv; Source_Data/Figure_4_samhd1_leave_one_species.tsv"),
        ("Supplementary Figure 1", "Supplementary_Data/Supplementary_Data_1_species_traits.tsv; Supplementary_Data/Supplementary_Data_13_datelife_node_calibration_audit.tsv; Supplementary_Data/Supplementary_Data_14_datelife_species_coverage.tsv"),
        ("Supplementary Figure 2", "Supplementary_Data/Supplementary_Data_6_branch_length_sensitivity.tsv"),
        ("Supplementary Figure 3", "Supplementary_Data/Supplementary_Data_41_sequence_validation_model_impact.tsv"),
        ("Supplementary Figure 4a-d", "Supplementary_Data/Supplementary_Data_4_gene_tree_summary.tsv; Supplementary_Data/Supplementary_Data_5_trimmed_alignment_qc.tsv; Supplementary_Data/Gene_Trees/"),
        ("Supplementary Figure 4e-h", "Supplementary_Data/Supplementary_Data_4_gene_tree_summary.tsv; Supplementary_Data/Gene_Trees/"),
        ("Supplementary Figure 5", "Supplementary_Data/Supplementary_Data_11_species_influence_pgls.tsv"),
        ("Supplementary Figure 6", "Supplementary_Data/Supplementary_Data_15_targeted_protein_conservation_rows.tsv; Supplementary_Data/Supplementary_Data_16_targeted_protein_conservation_species.tsv; Supplementary_Data/Supplementary_Data_17_targeted_protein_conservation_pgls.tsv"),
        ("Supplementary Figure 7", "Supplementary_Data/Supplementary_Data_38_samhd1_alignment_position_qc.tsv; Supplementary_Data/Supplementary_Data_39_samhd1_alignment_species_qc.tsv; Supplementary_Data/Supplementary_Data_40_samhd1_alignment_sensitivity_pgls.tsv"),
    ]
    with (TABLES / "jme_figure_source_data_map.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["figure", "source_data_files"])
        writer.writerows(figure_map)


def main() -> None:
    build_summary_tables()
    build_workflow_figure()
    build_module_forest()
    build_orthology_ladder()
    build_samhd1_figure()
    write_figure_source_data()
    print("Built JME figures, summary tables, and source manifest.")


if __name__ == "__main__":
    main()
