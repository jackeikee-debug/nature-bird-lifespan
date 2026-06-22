"""Write Week 6 synthesis deliverables."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DEFAULT_OUT_REPORT = pathlib.Path("results/reports/week6_feasibility_synthesis_report.md")
DEFAULT_FIG_PANEL = pathlib.Path("results/tables/week6_figure_panel.tsv")
DEFAULT_EXPANDED_PLAN = pathlib.Path("docs/expanded_panel_v2_plan.md")


def fmt(x: float) -> str:
    return f"{x:.3g}"


def read_tsv(path: str | pathlib.Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_OUT_REPORT)
    parser.add_argument("--figure-panel", type=pathlib.Path, default=DEFAULT_FIG_PANEL)
    parser.add_argument("--expanded-plan", type=pathlib.Path, default=DEFAULT_EXPANDED_PLAN)
    args = parser.parse_args()

    week4 = read_tsv("results/tables/week4_full_sequence_validated_transposon_summary.tsv")
    sens = read_tsv("results/tables/week4_full_sequence_sensitivity_summary.tsv")
    module_stats = read_tsv("results/tables/week5_module_enrichment_tests.tsv")
    background_stats = read_tsv("results/tables/week5_background_comparison_tests.tsv")
    perm = read_tsv("results/tables/week5_permutation_tests.tsv")
    human_translation = read_tsv("data/processed/human_translation_priority.tsv")
    genome = read_tsv("data/processed/genome_availability_audit.tsv")

    strict = week4[
        (week4["score_variant"] == "transposon_sequence_strict")
        & (week4["model"] == "mass_clade_module")
    ].iloc[0]
    weak = week4[
        (week4["score_variant"] == "transposon_sequence_weak_inclusive")
        & (week4["model"] == "mass_clade_module")
    ].iloc[0]
    bird_loss = sens[
        (sens["test"] == "subset_sensitivity")
        & (sens["comparison"] == "leave_out_Aves")
    ].iloc[0]
    tier1_loss = sens[
        (sens["test"] == "subset_sensitivity")
        & (sens["comparison"] == "tier1_only")
    ].iloc[0]
    trans_repeat = module_stats[
        (module_stats["comparison_a"] == "transposon_suppression")
        & (module_stats["feature"] == "repeat_control")
    ].iloc[0]
    bg_trans_repeat = background_stats[
        (background_stats["comparison_a"] == "transposon_5")
        & (background_stats["feature"] == "repeat_control")
    ].iloc[0]
    lead_perm = perm[
        (perm["test"] == "lead_TRIM28_SETDB1_MOV10_vs_random3")
        & (perm["universe"] == "maintenance_41")
    ].iloc[0]
    genome_counts = genome["genome_analysis_tier"].value_counts().to_dict()
    missing_assembly = int((genome["genome_available"] == "no").sum())
    trans = human_translation[
        human_translation["maintenance_module"] == "transposon_suppression"
    ].copy()
    trans_order = {"TRIM28": 1, "SETDB1": 2, "MOV10": 3, "PIWIL1": 4, "PIWIL2": 5}
    trans["trans_order"] = trans["human_gene_symbol"].map(trans_order).fillna(99)
    trans = trans.sort_values(["trans_order", "human_gene_symbol"])

    figure_rows = [
        {
            "figure_id": "Fig1",
            "role": "main",
            "file_png": "results/figures/week4_lifespan_residual_vs_strict_transposon_score.png",
            "file_pdf": "results/figures/week4_lifespan_residual_vs_strict_transposon_score.pdf",
            "caption_short": "Strict sequence-supported transposon score versus lifespan residual.",
            "message": "Shows positive visual association; formal inference comes from PGLS.",
        },
        {
            "figure_id": "Fig2",
            "role": "main",
            "file_png": "results/figures/week4_transposon_main_effect_sensitivity_forest.png",
            "file_pdf": "results/figures/week4_transposon_main_effect_sensitivity_forest.pdf",
            "caption_short": "Main effect and sensitivity forest plot.",
            "message": "Shows the strict signal survives most sensitivity checks but fails leave-out-birds and Tier1-only.",
        },
        {
            "figure_id": "Fig3",
            "role": "main",
            "file_png": "results/figures/week4_sequence_validation_waterfall_counts.png",
            "file_pdf": "results/figures/week4_sequence_validation_waterfall_counts.pdf",
            "caption_short": "Week 4 sequence-validation counts.",
            "message": "Documents sequence support and unresolved rows for transposon orthology.",
        },
        {
            "figure_id": "Fig4",
            "role": "main_or_supplement",
            "file_png": "results/figures/week5_translational_evidence_map.png",
            "file_pdf": "results/figures/week5_translational_evidence_map.pdf",
            "caption_short": "Human translational evidence map.",
            "message": "Connects focal transposon genes to human repeat-control, disease, and tractability evidence.",
        },
        {
            "figure_id": "FigS1",
            "role": "supplement",
            "file_png": "results/figures/week4_clade_colored_transposon_score_distribution.png",
            "file_pdf": "results/figures/week4_clade_colored_transposon_score_distribution.pdf",
            "caption_short": "Strict transposon score distribution by clade.",
            "message": "Supports interpretation of bird-heavy and annotation-tier-sensitive score distribution.",
        },
    ]
    fig_panel = pd.DataFrame(figure_rows)
    args.figure_panel.parent.mkdir(parents=True, exist_ok=True)
    fig_panel.to_csv(args.figure_panel, sep="\t", index=False)

    report_lines = [
        "# Week 6 Feasibility Synthesis Report",
        "",
        "## Executive Decision",
        "",
        "**Decision: cautious-go.** The project has a coherent and statistically supported feasibility signal, but it is not yet ready to be framed as a full Nature-family comparative genomics paper without an expanded validation stage.",
        "",
        "The strongest current claim is not broad flight convergence across all vertebrates. The strongest claim is a bird-dependent, sequence-supported transposon-suppression association with longer-than-expected lifespan, with human repeat-control and disease-context support.",
        "",
        "## Core Results",
        "",
        f"- Species/lifespan base: 240 species curated, 237 model-ready species in Week 1.",
        f"- Genome mechanism panel: 68 primary species with Tier1/Tier2 usable genome annotation for the 41-gene maintenance seed panel.",
        f"- Genome availability limitation: {missing_assembly} species without assembly; tier counts: {genome_counts}.",
        f"- Strict transposon PGLS mass+clade: estimate={fmt(float(strict['module_estimate']))}, p={fmt(float(strict['module_p']))}, BH={fmt(float(strict['module_p_bh_by_variant_model']))}, n={int(strict['n'])}.",
        f"- Weak-inclusive transposon PGLS mass+clade: estimate={fmt(float(weak['module_estimate']))}, p={fmt(float(weak['module_p']))}, BH={fmt(float(weak['module_p_bh_by_variant_model']))}.",
        f"- Leave-out-birds sensitivity fails: estimate={fmt(float(bird_loss['estimate']))}, p={fmt(float(bird_loss['p']))}, BH={fmt(float(bird_loss['p_bh_by_variant_model']))}.",
        f"- Tier1-only sensitivity fails: estimate={fmt(float(tier1_loss['estimate']))}, p={fmt(float(tier1_loss['p']))}, BH={fmt(float(tier1_loss['p_bh_by_variant_model']))}.",
        "",
        "## Orthology and Sequence Evidence",
        "",
        "- Full sequence-validated transposon matrix: 297 sequence-supported rows, 3 weak-supported rows, and 40 not-supported rows across 340 transposon gene-species slots.",
        "- Direct NCBI candidates were mostly reliable: 240/244 reciprocal-supported in the direct cross-check.",
        "- Main residual risk remains in rescue/unresolved PIWIL2 and TRIM28 rows, especially PIWIL2 paralog ambiguity and low-coverage TRIM28-like hits.",
        "",
        "## Human Translation Layer",
        "",
        "- Human mapping covered 41/41 maintenance seed genes.",
        "- Open Targets disease/tractability mapping succeeded for 41/41 genes and produced 410 top disease-association rows.",
        "- Lead translational genes: TRIM28, SETDB1, MOV10.",
        "- Mechanistic background genes: PIWIL1, PIWIL2.",
    ]
    for _, row in trans.iterrows():
        sm = row["small_molecule_tractability"] if isinstance(row["small_molecule_tractability"], str) and row["small_molecule_tractability"] else "none"
        report_lines.append(
            f"- {row['human_gene_symbol']}: max Open Targets score={fmt(float(row['max_open_targets_score']))}; small_molecule={sm}; priority={row['week5_translation_priority']}."
        )
    report_lines.extend(
        [
            "",
            "## Week 5 Statistical Support",
            "",
            f"- Transposon_suppression is enriched for repeat-control evidence within the maintenance seed set: {int(trans_repeat['a_hits'])}/{int(trans_repeat['a_n'])} vs {int(trans_repeat['b_hits'])}/{int(trans_repeat['b_n'])}, p={fmt(float(trans_repeat['p_value']))}, BH={fmt(float(trans_repeat['p_bh']))}.",
            f"- Five transposon genes are enriched for repeat-control evidence against random HGNC protein-coding background genes: {int(bg_trans_repeat['a_hits'])}/{int(bg_trans_repeat['a_n'])} vs {int(bg_trans_repeat['b_hits'])}/{int(bg_trans_repeat['b_n'])}, p={fmt(float(bg_trans_repeat['p_value']))}, BH={fmt(float(bg_trans_repeat['p_bh']))}.",
            f"- TRIM28/SETDB1/MOV10 are not statistically exceptional within the 41-gene maintenance universe by the current composite-score permutation test: p={fmt(float(lead_perm['p_value']))}, BH={fmt(float(lead_perm['p_bh']))}.",
            "",
            "## Claim Boundary",
            "",
            "Supported:",
            "",
            "- A strict sequence-supported transposon-suppression score is positively associated with longer-than-expected lifespan in a 68-species primary genome panel after body-mass and clade adjustment.",
            "- The signal is not driven by a single transposon gene, human, bats, reptiles, non-flying mammals, or the five largest lifespan-residual outliers.",
            "- The transposon module is specifically enriched for repeat-control evidence in human mapping, both within the maintenance seed set and against random protein-coding genes.",
            "",
            "Not yet supported:",
            "",
            "- A universal vertebrate transposon-lifespan effect.",
            "- A bird-independent effect.",
            "- A Tier1-only annotation-robust effect.",
            "- A claim that TRIM28/SETDB1/MOV10 are statistically exceptional among all maintenance genes.",
            "- A claim that focal transposon genes directly cause human longevity.",
            "",
            "## Figure Panel",
            "",
            f"Figure panel table: `{args.figure_panel}`",
            "",
            "- Main Figure 1: lifespan residual vs strict transposon score.",
            "- Main Figure 2: main effect and sensitivity forest plot.",
            "- Main Figure 3: sequence-validation waterfall/count plot.",
            "- Main/Supplement Figure 4: human translational evidence map.",
            "- Supplement Figure S1: clade-colored transposon score distribution.",
            "",
            "## Go / No-Go",
            "",
            "The decision is **cautious-go**.",
            "",
            "Proceed to a second-stage project if the next phase focuses on expanded validation rather than overclaiming the current seed panel. The current evidence is strong enough to justify more compute, more orthology validation, and a preprint-style manuscript draft. It is not yet strong enough to claim a Nature-family-ready mechanistic discovery.",
            "",
            "## Required Next Validation",
            "",
            "- Expanded panel v2: 200-300 genes focused on repeat suppression, chromatin repression, DNA repair, autophagy/mitophagy, inflammation restraint, and cancer surveillance.",
            "- Better phylogeny: replace OpenTree+Grafen with a dated species tree or TimeTree-informed branch lengths where possible.",
            "- Orthology hardening: OMA/OrthoDB/Ensembl Compara or manual domain checks for the 39-row ambiguous queue.",
            "- Annotation-bias controls: explicitly model genome tier, annotation completeness, and module missingness.",
            "- Independent gene-set validation: test whether GO/Reactome transposon silencing and chromatin repression gene sets reproduce the current signal.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    expanded_lines = [
        "# Expanded Panel v2 Plan",
        "",
        "## Purpose",
        "",
        "Expanded panel v2 tests whether the Week 4 transposon-suppression signal survives beyond the 41-gene seed panel. The goal is not to add every ageing gene immediately, but to add enough curated genes to distinguish a specific repeat/chromatin signal from a broad genome-maintenance or annotation-quality artifact.",
        "",
        "## Target Size",
        "",
        "- Minimum: 150 genes",
        "- Preferred: 200-300 genes",
        "- Avoid immediate 2000+ gene expansion until the 200-300 gene panel is stable.",
        "",
        "## Proposed Modules",
        "",
        "1. Transposon and repeat suppression",
        "   - GO: piRNA metabolic process, transposon silencing, heterochromatin formation",
        "   - Include PIWI pathway, KRAB-ZNF cofactors, HUSH complex, LINE-1 restriction factors",
        "",
        "2. Chromatin repression and heterochromatin maintenance",
        "   - SETDB1/TRIM28 axis, HP1/CBX proteins, DNA methylation readers/writers, histone methyltransferases",
        "",
        "3. DNA repair and replication stress",
        "   - Reactome DNA repair, double-strand break repair, NHEJ, homologous recombination, Fanconi anemia pathway",
        "",
        "4. Proteostasis and autophagy/mitophagy",
        "   - chaperones, proteasome, ATG core, PINK1/PRKN mitophagy, mitochondrial dynamics",
        "",
        "5. Inflammation restraint and innate immune sensing",
        "   - NF-kB regulators, inflammasome, cGAS-STING, interferon pathway negative regulators",
        "",
        "6. Cancer surveillance and senescence control",
        "   - TP53/RB/CDKN pathways, apoptosis, senescence regulators, tumor suppressor modules",
        "",
        "## Data Sources",
        "",
        "- Reactome pathway gene sets",
        "- GO Biological Process gene sets",
        "- GenAge and CellAge",
        "- MSigDB Hallmark/CP where license permits",
        "- UniProt keywords for DNA repair, chromatin, transposon/repeat biology",
        "- HGNC families for PIWI, KRAB-ZNF, chromatin regulators",
        "",
        "## Analysis Plan",
        "",
        "- Build `maintenance_gene_sets_v2.tsv` with source and module labels.",
        "- Recompute ortholog matrix for the 68 primary species first.",
        "- Score modules using coverage-aware and confidence-weighted scores.",
        "- Re-run PGLS with body mass, clade, and annotation-tier sensitivity.",
        "- Run leave-one-submodule-out tests for transposon/repeat suppression.",
        "- Compare module ranks against random matched gene sets.",
        "",
        "## Success Criteria",
        "",
        "- Repeat/chromatin module remains positive in mass+clade PGLS.",
        "- Signal survives removal of obvious annotation-tier artifacts.",
        "- Signal is not erased by removing birds with weaker assemblies.",
        "- Repeat/chromatin modules rank above generic disease/annotation-rich modules in permutation tests.",
        "",
        "## Stop Criteria",
        "",
        "- Signal disappears after annotation-tier adjustment.",
        "- Many unrelated gene sets show equal or stronger associations, suggesting annotation completeness artifact.",
        "- Orthology ambiguity becomes too high for repeat/chromatin genes without a stronger domain-level pipeline.",
    ]
    args.expanded_plan.parent.mkdir(parents=True, exist_ok=True)
    args.expanded_plan.write_text("\n".join(expanded_lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.report}, {args.figure_panel}, and {args.expanded_plan}")


if __name__ == "__main__":
    main()
