"""Create a polished phylogenetic heatmap candidate for Figure 1.

The figure is intentionally descriptive: it combines the primary genome-panel
tree with lifespan residuals, transposon/repeat scores, coverage, genome tier,
and targeted sequence-validation evidence. It is not used as a statistical test.
"""

from __future__ import annotations

import argparse
import copy
import math
import pathlib
from dataclasses import dataclass

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Patch, Rectangle


DEFAULT_TREE = pathlib.Path("data/processed/phylogeny_inputs/opentree_induced_subtree.tre")
DEFAULT_LIFESPAN = pathlib.Path(
    "data/processed/maintenance_lifespan_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv"
)
DEFAULT_EVIDENCE = pathlib.Path("data/processed/phase3_evidence_ladder.tsv")
DEFAULT_TABLE = pathlib.Path("data/processed/figure1_phylo_heatmap_table.tsv")
DEFAULT_COMPACT_TABLE = pathlib.Path("data/processed/figure1_phylo_heatmap_compact_table.tsv")
DEFAULT_PNG = pathlib.Path("results/figures/figure1_phylogenetic_heatmap_polished.png")
DEFAULT_PDF = pathlib.Path("results/figures/figure1_phylogenetic_heatmap_polished.pdf")
DEFAULT_SVG = pathlib.Path("results/figures/figure1_phylogenetic_heatmap_polished.svg")
DEFAULT_FULL_PNG = pathlib.Path("results/figures/figureS1_phylogenetic_heatmap_full.png")
DEFAULT_FULL_PDF = pathlib.Path("results/figures/figureS1_phylogenetic_heatmap_full.pdf")
DEFAULT_FULL_SVG = pathlib.Path("results/figures/figureS1_phylogenetic_heatmap_full.svg")
DEFAULT_REPORT = pathlib.Path("results/reports/figure1_phylogenetic_heatmap_report.md")
DEFAULT_VARIANT = "phase2_W3_full_background_sensitivity"


CLADE_COLORS = {
    "Aves": "#0F766E",
    "Mammalia_Chiroptera": "#7C3AED",
    "Mammalia_nonChiroptera": "#475569",
    "Reptilia": "#B7791F",
    "Anchor": "#111827",
}

TIER_COLORS = {
    "tier1_refseq_annotated_chromosome": "#176B87",
    "tier1_refseq_annotated": "#2A9D8F",
    "tier2_annotated": "#E9A23B",
    "tier3_assembly_only": "#9CA3AF",
}

EVIDENCE_STRICT_LEVELS = {
    "local_gff_protein_strict",
    "local_cds_translation_strict",
    "external_uniprot_strict",
}

EVIDENCE_AMBIGUOUS_LEVELS = {
    "dnmt1_partial_fragment_not_absence",
    "dnmt1_longer_local_isoform_available",
    "mbd2_short_mbd_domain_ambiguous",
    "mbd3_partial_fragment_not_absence",
    "local_gff_forward_supported_manual_review",
    "local_gff_rejected_no_same_gene_reference",
    "gff_annotation_only_pending_sequence",
    "gff_probable_wrong_gene_family",
}

EVIDENCE_UNRESOLVED_LEVELS = {
    "not_found_in_local_gff",
}


@dataclass
class TreeLayout:
    x: dict[int, float]
    y: dict[int, float]
    terminals: list
    max_x: float


def safe_read_tsv(path: pathlib.Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def collapse_unary(clade) -> None:
    """Collapse one-child internal nodes after pruning an induced subtree."""
    for child in list(clade.clades):
        collapse_unary(child)

    while len(clade.clades) == 1:
        child = clade.clades[0]
        clade.name = child.name
        clade.clades = child.clades


def prune_tree_to_labels(tree, labels: set[str]):
    tree = copy.deepcopy(tree)
    changed = True
    while changed:
        changed = False
        for terminal in list(tree.get_terminals()):
            if terminal.name not in labels:
                tree.prune(terminal)
                changed = True
    collapse_unary(tree.root)
    tree.ladderize()
    return tree


def layout_tree(tree) -> TreeLayout:
    terminals = tree.get_terminals()
    y = {id(term): i for i, term in enumerate(terminals)}
    x: dict[int, float] = {}

    def assign_x(clade, depth: float) -> None:
        x[id(clade)] = depth
        for child in clade.clades:
            assign_x(child, depth + 1.0)

    def assign_internal_y(clade) -> float:
        if clade.is_terminal():
            return y[id(clade)]
        child_ys = [assign_internal_y(child) for child in clade.clades]
        y[id(clade)] = float(sum(child_ys) / len(child_ys))
        return y[id(clade)]

    assign_x(tree.root, 0.0)
    assign_internal_y(tree.root)
    max_x = max(x[id(term)] for term in terminals)
    return TreeLayout(x=x, y=y, terminals=terminals, max_x=max_x)


def species_label(name: str) -> str:
    parts = name.split()
    if len(parts) >= 2:
        return f"{parts[0][0]}. {parts[1]}"
    return name


def summarize_evidence(evidence: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for species, sub in evidence.groupby("scientific_name", dropna=False):
        strict = sub["phase3_evidence_level"].isin(EVIDENCE_STRICT_LEVELS).sum()
        ambiguous = sub["phase3_evidence_level"].isin(EVIDENCE_AMBIGUOUS_LEVELS).sum()
        unresolved = sub["phase3_evidence_level"].isin(EVIDENCE_UNRESOLVED_LEVELS).sum()
        total = len(sub)
        rows.append(
            {
                "scientific_name": species,
                "validation_rows": total,
                "strict_sequence_rows": strict,
                "ambiguous_or_review_rows": ambiguous,
                "unresolved_rows": unresolved,
                "strict_sequence_fraction": strict / total if total else np.nan,
                "unresolved_or_ambiguous_fraction": (ambiguous + unresolved) / total if total else np.nan,
                "targeted_sequence_validation": True,
            }
        )
    return pd.DataFrame(rows)


def build_plot_table(lifespan: pd.DataFrame, evidence: pd.DataFrame, variant: str) -> pd.DataFrame:
    base = lifespan[lifespan["score_variant"].eq(variant)].copy()
    if base.empty:
        raise ValueError(f"No rows found for score variant: {variant}")

    ev = summarize_evidence(evidence)
    merged = base.merge(ev, on="scientific_name", how="left")
    for col in [
        "validation_rows",
        "strict_sequence_rows",
        "ambiguous_or_review_rows",
        "unresolved_rows",
        "strict_sequence_fraction",
        "unresolved_or_ambiguous_fraction",
    ]:
        merged[col] = merged[col].astype(float)
    merged["targeted_sequence_validation"] = merged["targeted_sequence_validation"].fillna(False)
    merged["plot_label"] = merged["scientific_name"].map(species_label)
    return merged


def select_compact_table(table: pd.DataFrame) -> pd.DataFrame:
    """Select a manuscript-readable subset while retaining the validation core."""
    selected: set[str] = set(table.loc[table["targeted_sequence_validation"], "scientific_name"])
    anchors = {
        "Homo sapiens",
        "Mus musculus",
        "Heterocephalus glaber",
        "Gallus gallus",
        "Taeniopygia guttata",
        "Erythrura gouldiae",
        "Pteropus vampyrus",
        "Myotis brandtii",
        "Desmodus rotundus",
        "Alligator mississippiensis",
        "Crocodylus porosus",
        "Anolis carolinensis",
        "Chelonia mydas",
    }
    selected.update(set(table.loc[table["scientific_name"].isin(anchors), "scientific_name"]))

    quotas = {
        "Aves": 22,
        "Mammalia_Chiroptera": 7,
        "Mammalia_nonChiroptera": 8,
        "Reptilia": 8,
    }
    for clade, quota in quotas.items():
        current_n = int((table["scientific_name"].isin(selected) & table["clade"].eq(clade)).sum())
        need = max(0, quota - current_n)
        if need == 0:
            continue
        candidates = table[~table["scientific_name"].isin(selected) & table["clade"].eq(clade)].copy()
        if candidates.empty:
            continue
        candidates = candidates.sort_values(
            ["lifespan_residual_log10", "transposon_repeat_suppression_coverage"],
            ascending=[False, False],
        )
        selected.update(set(candidates.head(need)["scientific_name"]))

    return table[table["scientific_name"].isin(selected)].copy()


def draw_tree(ax, tree, layout: TreeLayout, terminal_x: float, row_lookup: dict[str, pd.Series]) -> None:
    def clade_color_for_terminal(term_name: str) -> str:
        row = row_lookup.get(term_name)
        if row is None:
            return "#6B7280"
        return CLADE_COLORS.get(str(row["clade"]), "#6B7280")

    def draw_clade(clade) -> None:
        x0 = layout.x[id(clade)]
        y0 = layout.y[id(clade)]
        if clade.clades:
            child_ys = [layout.y[id(child)] for child in clade.clades]
            ax.plot([x0, x0], [min(child_ys), max(child_ys)], color="#334155", lw=0.55, zorder=1)
            for child in clade.clades:
                x1 = layout.x[id(child)]
                y1 = layout.y[id(child)]
                ax.plot([x0, x1], [y1, y1], color="#334155", lw=0.55, zorder=1)
                draw_clade(child)
        else:
            color = clade_color_for_terminal(clade.name)
            ax.plot([x0, terminal_x], [y0, y0], color="#CBD5E1", lw=0.45, zorder=0)
            ax.scatter([terminal_x], [y0], s=22, color=color, edgecolor="white", linewidth=0.45, zorder=3)

    draw_clade(tree.root)


def add_heat_square(ax, x: float, y: float, value: float, cmap, norm, na_color="#E5E7EB") -> None:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        face = na_color
    else:
        face = cmap(norm(value))
    ax.add_patch(Rectangle((x - 0.38, y - 0.38), 0.76, 0.76, facecolor=face, edgecolor="white", lw=0.5))


def plot_figure(
    table: pd.DataFrame,
    tree_path: pathlib.Path,
    out_png: pathlib.Path,
    out_pdf: pathlib.Path,
    out_svg: pathlib.Path,
    title: str,
    subtitle: str,
    panel_label: str,
) -> None:
    labels = set(table["opentree_tip_label"].dropna())
    tree = Phylo.read(tree_path, "newick")
    tree = prune_tree_to_labels(tree, labels)
    layout = layout_tree(tree)

    order = [term.name for term in layout.terminals]
    row_lookup = table.set_index("opentree_tip_label", drop=False).to_dict(orient="index")
    table_ordered = pd.DataFrame([row_lookup[name] for name in order if name in row_lookup])

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    n = len(table_ordered)
    fig_h = max(10.6, n * 0.27)
    fig, ax = plt.subplots(figsize=(15.6, fig_h), dpi=300)
    terminal_x = layout.max_x + 0.4
    label_x = terminal_x + 0.55
    heat_start = terminal_x + 4.4
    heat_cols = [
        ("Residual\nlifespan", "lifespan_residual_log10", plt.cm.RdBu_r, TwoSlopeNorm(vmin=-0.45, vcenter=0, vmax=0.45)),
        ("Repeat\nscore", "transposon_repeat_suppression_score", plt.cm.YlGn, Normalize(vmin=0, vmax=1)),
        ("Repeat\ncoverage", "transposon_repeat_suppression_coverage", plt.cm.Blues, Normalize(vmin=0, vmax=1)),
        ("Strict seq.\nrows", "strict_sequence_rows", plt.cm.Greens, Normalize(vmin=0, vmax=10)),
        ("Ambig./\nunresolved", "unresolved_or_ambiguous_fraction", plt.cm.OrRd, Normalize(vmin=0, vmax=1)),
    ]

    draw_tree(ax, tree, layout, terminal_x, row_lookup)

    y_by_name = {term.name: layout.y[id(term)] for term in layout.terminals}

    # Subtle clade bands make the tree readable without replacing it with a legend-heavy plot.
    clade_ranges = []
    for clade, sub in table_ordered.groupby("clade", sort=False):
        ys = [y_by_name[name] for name in sub["opentree_tip_label"] if name in y_by_name]
        if ys:
            clade_ranges.append((clade, min(ys) - 0.42, max(ys) + 0.42))
    for clade, y0, y1 in clade_ranges:
        color = CLADE_COLORS.get(str(clade), "#6B7280")
        ax.add_patch(
            Rectangle((-0.12, y0), 0.09, y1 - y0, facecolor=color, edgecolor="none", alpha=0.95, zorder=4)
        )

    for _, row in table_ordered.iterrows():
        y = y_by_name[row["opentree_tip_label"]]
        clade_color = CLADE_COLORS.get(str(row["clade"]), "#6B7280")
        ax.text(
            label_x,
            y,
            row["plot_label"],
            va="center",
            ha="left",
            fontsize=6.2,
            color="#111827",
        )
        ax.add_patch(Rectangle((heat_start - 0.82, y - 0.38), 0.18, 0.76, facecolor=clade_color, edgecolor="none"))

        tier = str(row["genome_analysis_tier"])
        tier_color = TIER_COLORS.get(tier, "#9CA3AF")
        ax.add_patch(Rectangle((heat_start - 0.53, y - 0.38), 0.18, 0.76, facecolor=tier_color, edgecolor="none"))

        for i, (_, col, cmap, norm) in enumerate(heat_cols):
            na = "#F3F4F6" if col in {"strict_sequence_rows", "unresolved_or_ambiguous_fraction"} else "#E5E7EB"
            add_heat_square(ax, heat_start + i, y, row[col], cmap, norm, na_color=na)
            if col == "strict_sequence_rows" and bool(row["targeted_sequence_validation"]):
                ax.text(
                    heat_start + i,
                    y,
                    f"{int(row[col])}",
                    va="center",
                    ha="center",
                    fontsize=5.5,
                    color="#0F172A" if row[col] < 6 else "white",
                )

    title_y = -2.55
    subtitle_y = -1.55
    header_y = -0.72

    ax.text(0, title_y, panel_label, fontsize=15, fontweight="bold", ha="left", va="bottom")
    ax.text(
        0.6,
        title_y,
        title,
        fontsize=12,
        fontweight="bold",
        ha="left",
        va="bottom",
        color="#111827",
    )
    ax.text(
        0.6,
        subtitle_y,
        subtitle,
        fontsize=7.5,
        ha="left",
        va="bottom",
        color="#475569",
    )

    ax.text(heat_start - 0.73, header_y, "Clade", rotation=45, ha="right", va="bottom", fontsize=7)
    ax.text(heat_start - 0.44, header_y, "Genome\ntier", rotation=45, ha="right", va="bottom", fontsize=7)
    for i, (label, _, _, _) in enumerate(heat_cols):
        ax.text(heat_start + i, header_y, label, rotation=45, ha="right", va="bottom", fontsize=7)

    ax.set_ylim(n + 0.9, -3.0)
    ax.set_xlim(-0.2, heat_start + len(heat_cols) + 3.2)
    ax.axis("off")

    legend_x = heat_start + len(heat_cols) + 0.2
    ax.text(legend_x, 1.0, "Clade", fontsize=8, fontweight="bold", ha="left", va="center")
    for j, (label, color) in enumerate(
        [
            ("Birds", CLADE_COLORS["Aves"]),
            ("Bats", CLADE_COLORS["Mammalia_Chiroptera"]),
            ("Other mammals", CLADE_COLORS["Mammalia_nonChiroptera"]),
            ("Reptiles", CLADE_COLORS["Reptilia"]),
        ]
    ):
        yy = 2.1 + j * 1.0
        ax.scatter([legend_x + 0.1], [yy], s=28, color=color, edgecolor="white", linewidth=0.45)
        ax.text(legend_x + 0.35, yy, label, fontsize=7, va="center", ha="left")

    ax.text(legend_x, 7.0, "Genome tier", fontsize=8, fontweight="bold", ha="left", va="center")
    tier_labels = [
        ("Tier 1 chr.", TIER_COLORS["tier1_refseq_annotated_chromosome"]),
        ("Tier 1 annot.", TIER_COLORS["tier1_refseq_annotated"]),
        ("Tier 2 annot.", TIER_COLORS["tier2_annotated"]),
        ("Other/held", TIER_COLORS["tier3_assembly_only"]),
    ]
    for j, (label, color) in enumerate(tier_labels):
        yy = 8.1 + j * 1.0
        ax.add_patch(Rectangle((legend_x, yy - 0.28), 0.22, 0.56, facecolor=color, edgecolor="none"))
        ax.text(legend_x + 0.35, yy, label, fontsize=7, va="center", ha="left")

    ax.text(legend_x, 13.2, "Heatmap scales", fontsize=8, fontweight="bold", ha="left", va="center")
    scale_notes = [
        ("Residual lifespan", "blue = lower, red = higher"),
        ("Repeat score / coverage", "light = low, dark = high"),
        ("Strict sequence rows", "0 to 10 targeted rows"),
        ("Ambig./unresolved", "white = low, orange = high"),
    ]
    for j, (label, note) in enumerate(scale_notes):
        yy = 14.25 + j * 1.0
        ax.text(legend_x, yy, label, fontsize=6.8, ha="left", va="center", color="#111827")
        ax.text(legend_x, yy + 0.38, note, fontsize=6.1, ha="left", va="center", color="#64748B")

    fig.subplots_adjust(left=0.025, right=0.985, top=0.985, bottom=0.025)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    fig.savefig(out_svg)
    plt.close(fig)


def write_report(full_table: pd.DataFrame, compact_table: pd.DataFrame, report: pathlib.Path, variant: str) -> None:
    targeted = int(full_table["targeted_sequence_validation"].sum())
    clade_counts = full_table["clade"].value_counts().to_dict()
    compact_clade_counts = compact_table["clade"].value_counts().to_dict()
    tier_counts = full_table["genome_analysis_tier"].value_counts().to_dict()
    rows = [
        "# Figure 1 Phylogenetic Heatmap Report",
        "",
        "This figure is a polished manuscript candidate that combines the primary genome-panel topology with lifespan residuals, transposon/repeat scores, coverage, genome tier, and targeted sequence-validation evidence.",
        "",
        "## Inputs",
        "",
        f"- score variant: `{variant}`",
        f"- compact Figure 1 species plotted: {len(compact_table)}",
        f"- full supplementary species plotted: {len(full_table)}",
        f"- targeted sequence-validation species: {targeted}",
        f"- compact clade counts: {compact_clade_counts}",
        f"- full clade counts: {clade_counts}",
        f"- genome tier counts: {tier_counts}",
        "",
        "## Interpretation",
        "",
        "The compact figure is descriptive and should be used as a visual overview, not as a standalone statistical test. It is suitable as a Figure 1 candidate because it links the comparative species framework to the sequence-audited transposon/repeat module. The full 68-species version is better suited to supplementary display.",
        "",
        "## Claim Boundary",
        "",
        "The tree uses the OpenTree synthetic topology and cladogram-depth display. Strong evolutionary timing claims still require dated or TimeTree-informed branch lengths.",
    ]
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", type=pathlib.Path, default=DEFAULT_TREE)
    parser.add_argument("--lifespan", type=pathlib.Path, default=DEFAULT_LIFESPAN)
    parser.add_argument("--evidence", type=pathlib.Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--score-variant", default=DEFAULT_VARIANT)
    parser.add_argument("--table", type=pathlib.Path, default=DEFAULT_TABLE)
    parser.add_argument("--compact-table", type=pathlib.Path, default=DEFAULT_COMPACT_TABLE)
    parser.add_argument("--png", type=pathlib.Path, default=DEFAULT_PNG)
    parser.add_argument("--pdf", type=pathlib.Path, default=DEFAULT_PDF)
    parser.add_argument("--svg", type=pathlib.Path, default=DEFAULT_SVG)
    parser.add_argument("--full-png", type=pathlib.Path, default=DEFAULT_FULL_PNG)
    parser.add_argument("--full-pdf", type=pathlib.Path, default=DEFAULT_FULL_PDF)
    parser.add_argument("--full-svg", type=pathlib.Path, default=DEFAULT_FULL_SVG)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    lifespan = safe_read_tsv(args.lifespan)
    evidence = safe_read_tsv(args.evidence)
    table = build_plot_table(lifespan, evidence, args.score_variant)
    compact_table = select_compact_table(table)
    args.table.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.table, sep="\t", index=False)
    args.compact_table.parent.mkdir(parents=True, exist_ok=True)
    compact_table.to_csv(args.compact_table, sep="\t", index=False)
    plot_figure(
        compact_table,
        args.tree,
        args.png,
        args.pdf,
        args.svg,
        "Genome-panel phylogeny links lifespan residuals to sequence-audited repeat-control scores",
        "Compact manuscript view: all targeted validation birds plus anchor and clade-representative species. OpenTree topology shown as cladogram depth.",
        "A",
    )
    plot_figure(
        table,
        args.tree,
        args.full_png,
        args.full_pdf,
        args.full_svg,
        "Full 68-species genome-panel phylogeny with repeat-control score and validation evidence",
        "Supplementary full-panel view. OpenTree synthetic topology; branch lengths are shown as cladogram depth.",
        "S1",
    )
    write_report(table, compact_table, args.report, args.score_variant)
    print(
        f"Wrote {args.png}, {args.pdf}, {args.svg}, {args.full_png}, {args.full_pdf}, "
        f"{args.full_svg}, {args.table}, {args.compact_table}, and {args.report}"
    )


if __name__ == "__main__":
    main()
