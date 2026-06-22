"""Write Week 5 interpretation and manuscript-framing report."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DEFAULT_HUMAN = pathlib.Path("data/processed/human_mapping.tsv")
DEFAULT_TRANSLATION = pathlib.Path("data/processed/human_translation_priority.tsv")
DEFAULT_FIG_TABLE = pathlib.Path("data/processed/week5_translational_evidence_plot_table.tsv")
DEFAULT_MODULE_STATS = pathlib.Path("results/tables/week5_module_enrichment_tests.tsv")
DEFAULT_BACKGROUND_STATS = pathlib.Path("results/tables/week5_background_comparison_tests.tsv")
DEFAULT_PERM_STATS = pathlib.Path("results/tables/week5_permutation_tests.tsv")
DEFAULT_OUTPUT = pathlib.Path("results/reports/week5_interpretation_report.md")

FOCAL_GENES = ["TRIM28", "SETDB1", "MOV10", "PIWIL1", "PIWIL2"]
LEAD_GENES = ["TRIM28", "SETDB1", "MOV10"]
BACKGROUND_GENES = ["PIWIL1", "PIWIL2"]


def yes(value: object) -> bool:
    return str(value).lower() in {"yes", "1", "true", "keyword_supported"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--human", type=pathlib.Path, default=DEFAULT_HUMAN)
    parser.add_argument("--translation", type=pathlib.Path, default=DEFAULT_TRANSLATION)
    parser.add_argument("--fig-table", type=pathlib.Path, default=DEFAULT_FIG_TABLE)
    parser.add_argument("--module-stats", type=pathlib.Path, default=DEFAULT_MODULE_STATS)
    parser.add_argument("--background-stats", type=pathlib.Path, default=DEFAULT_BACKGROUND_STATS)
    parser.add_argument("--perm-stats", type=pathlib.Path, default=DEFAULT_PERM_STATS)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    human = pd.read_csv(args.human, sep="\t")
    translation = pd.read_csv(args.translation, sep="\t")
    fig = pd.read_csv(args.fig_table, sep="\t")
    merged = human.merge(
        translation[
            [
                "human_gene_symbol",
                "max_open_targets_score",
                "small_molecule_tractability",
                "top_disease_names",
                "week5_translation_priority",
            ]
        ],
        on="human_gene_symbol",
        how="left",
    ).merge(
        fig[
            [
                "human_gene_symbol",
                "week4_focal_transposon_gene",
                "repeat_control",
                "senescence",
                "small_molecule",
            ]
        ],
        on="human_gene_symbol",
        how="left",
    )

    focal = merged[merged["human_gene_symbol"].isin(FOCAL_GENES)].copy()
    lead = merged[merged["human_gene_symbol"].isin(LEAD_GENES)].copy()
    background = merged[merged["human_gene_symbol"].isin(BACKGROUND_GENES)].copy()
    module_stats = pd.read_csv(args.module_stats, sep="\t") if args.module_stats.exists() else pd.DataFrame()
    background_stats = pd.read_csv(args.background_stats, sep="\t") if args.background_stats.exists() else pd.DataFrame()
    perm_stats = pd.read_csv(args.perm_stats, sep="\t") if args.perm_stats.exists() else pd.DataFrame()

    high_priority = int((translation["week5_translation_priority"] == "high_translation_priority").sum())
    medium_priority = int((translation["week5_translation_priority"] == "medium_translation_priority").sum())
    curated_counts = {
        "GenAge human": int((human["genage_human_evidence"] == "yes").sum()),
        "LongevityMap": int((human["longevitymap_evidence"] == "yes").sum()),
        "CellAge": int((human["cellage_evidence"] == "yes").sum()),
    }
    focal_repeat = int((focal["repeat_control"] == 1).sum())
    focal_small = int((focal["small_molecule"] == 1).sum())

    lines = [
        "# Week 5 Interpretation Report",
        "",
        "## Position",
        "",
        "Week 5 connects the Week 4 comparative transposon-suppression signal to human ageing, disease, and tractability evidence. The result supports a translational framing, but not a direct claim that the focal transposon genes are proven human longevity genes.",
        "",
        "## Inputs",
        "",
        f"- Human ageing map: `{args.human}`",
        f"- Disease/tractability gene summary: `{args.translation}`",
        f"- Translational evidence plot table: `{args.fig_table}`",
        "",
        "## Summary Numbers",
        "",
        f"- Human maintenance genes mapped: {len(human)}",
        f"- High translation-priority genes: {high_priority}",
        f"- Medium translation-priority genes: {medium_priority}",
        f"- GenAge human hits: {curated_counts['GenAge human']}",
        f"- LongevityMap hits: {curated_counts['LongevityMap']}",
        f"- CellAge hits: {curated_counts['CellAge']}",
        f"- Focal transposon genes with repeat-control support: {focal_repeat}/5",
        f"- Focal transposon genes with small-molecule tractability: {focal_small}/5",
        "",
        "## Lead Translational Genes",
        "",
    ]
    for _, row in lead.sort_values("max_open_targets_score", ascending=False).iterrows():
        sm = row["small_molecule_tractability"] if isinstance(row["small_molecule_tractability"], str) and row["small_molecule_tractability"] else "none"
        lines.append(
            f"- `{row['human_gene_symbol']}`: priority={row['week5_translation_priority']}; "
            f"OpenTargetsMax={float(row['max_open_targets_score']):.3f}; "
            f"small_molecule={sm}; top_diseases={row['top_disease_names']}."
        )
    lines.extend(
        [
            "",
            "Interpretation: TRIM28, SETDB1, and MOV10 are the strongest translational bridge genes. They connect repeat/chromatin regulation to human disease spaces, and TRIM28/SETDB1 also carry small-molecule tractability annotations in Open Targets.",
            "",
            "## Mechanistic Background Genes",
            "",
        ]
    )
    for _, row in background.sort_values("max_open_targets_score", ascending=False).iterrows():
        lines.append(
            f"- `{row['human_gene_symbol']}`: priority={row['week5_translation_priority']}; "
            f"OpenTargetsMax={float(row['max_open_targets_score']):.3f}; "
            f"top_diseases={row['top_disease_names']}."
        )
    lines.extend(
        [
            "",
            "Interpretation: PIWIL1 and PIWIL2 remain important for the biological logic of piRNA/repeat suppression, but the current human disease and tractability layer is weaker. They should support the mechanism rather than carry the translational claim.",
            "",
            "## Claims Supported Now",
            "",
            "- The Week 4 comparative signal can be connected to human genes involved in repeat control, genome maintenance, senescence context, and disease associations.",
            "- The strongest translational hook is chromatin/repression and repeat-control biology, especially TRIM28, SETDB1, and MOV10.",
            "- The human layer supports relevance to disease spaces such as cancer, neurodegenerative disease, fertility/reproductive phenotypes, immune/inflammatory disease, and genome-instability syndromes.",
            "- Open Targets tractability supports prioritizing TRIM28 and SETDB1 for targetability discussion, while MOV10 is stronger as disease/mechanism context.",
            "",
            "## Statistical Support",
            "",
        ]
    )
    if module_stats.empty or background_stats.empty or perm_stats.empty:
        lines.append("- Statistical support tables have not yet been generated.")
    else:
        trans_repeat = module_stats[
            (module_stats["comparison_a"] == "transposon_suppression")
            & (module_stats["feature"] == "repeat_control")
        ]
        bg_trans_repeat = background_stats[
            (background_stats["comparison_a"] == "transposon_5")
            & (background_stats["feature"] == "repeat_control")
        ]
        lead_perm = perm_stats[
            (perm_stats["test"] == "lead_TRIM28_SETDB1_MOV10_vs_random3")
            & (perm_stats["universe"] == "maintenance_41")
        ]
        if not trans_repeat.empty:
            row = trans_repeat.iloc[0]
            lines.append(
                f"- Transposon_suppression is enriched for repeat-control evidence within the maintenance seed set: {int(row['a_hits'])}/{int(row['a_n'])} vs {int(row['b_hits'])}/{int(row['b_n'])}, p={float(row['p_value']):.3g}, BH={float(row['p_bh']):.3g}."
            )
        if not bg_trans_repeat.empty:
            row = bg_trans_repeat.iloc[0]
            lines.append(
                f"- The five transposon genes are enriched for repeat-control evidence against random HGNC protein-coding background genes: {int(row['a_hits'])}/{int(row['a_n'])} vs {int(row['b_hits'])}/{int(row['b_n'])}, p={float(row['p_value']):.3g}, BH={float(row['p_bh']):.3g}."
            )
        if not lead_perm.empty:
            row = lead_perm.iloc[0]
            lines.append(
                f"- The TRIM28/SETDB1/MOV10 lead trio is not statistically exceptional within the 41 maintenance genes by the current composite-score permutation test: observed={float(row['observed_score']):.3f}, p={float(row['p_value']):.3g}, BH={float(row['p_bh']):.3g}."
            )
    lines.extend(
        [
            "",
            "## Claims Not Supported Yet",
            "",
            "- Do not claim that the focal transposon genes directly explain human longevity.",
            "- Do not claim causal ageing effects from Open Targets disease association scores.",
            "- Do not claim all transposon-suppression genes are equally targetable; PIWIL1 and PIWIL2 currently lack small-molecule tractability in this pass.",
            "- Do not claim universal vertebrate conservation from Week 4, because the comparative signal remains bird-dependent and Tier1-sensitive.",
            "",
            "## Manuscript Framing",
            "",
            "A conservative manuscript sentence would be:",
            "",
            "> A sequence-supported transposon-suppression score is associated with longer-than-expected lifespan in a bird-heavy comparative panel, and the focal human orthologues connect this signal to repeat-control, chromatin repression, disease association, and selected tractability evidence, particularly through TRIM28, SETDB1, and MOV10.",
            "",
            "A stronger but still defensible framing for the Discussion would be:",
            "",
            "> These results suggest that long-lived flying lineages may have converged on enhanced somatic genome-maintenance programs, with transposon suppression emerging as a prioritized module that is mechanistically interpretable in human ageing-adjacent disease biology.",
            "",
            "## Week 6 Readiness",
            "",
            "Week 5 is ready to hand off to Week 6. The next step is not another broad data pass, but manuscript-level synthesis: final feasibility conclusion, figure panel selection, limitations, and a decision on whether to expand orthology validation or write a preprint-style report.",
        ]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
