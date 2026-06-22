"""Write manuscript outline v2 and executive summary."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DEFAULT_SYNTHESIS = pathlib.Path("results/reports/week6_feasibility_synthesis_report.md")
DEFAULT_FIG_PANEL = pathlib.Path("results/tables/week6_figure_panel.tsv")
DEFAULT_OUTLINE = pathlib.Path("docs/manuscript_outline_v2.md")
DEFAULT_EXEC = pathlib.Path("results/reports/week6_executive_summary.md")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--synthesis", type=pathlib.Path, default=DEFAULT_SYNTHESIS)
    parser.add_argument("--figure-panel", type=pathlib.Path, default=DEFAULT_FIG_PANEL)
    parser.add_argument("--outline", type=pathlib.Path, default=DEFAULT_OUTLINE)
    parser.add_argument("--executive-summary", type=pathlib.Path, default=DEFAULT_EXEC)
    args = parser.parse_args()

    if not args.synthesis.exists():
        raise FileNotFoundError(f"Missing synthesis report: {args.synthesis}")
    if not args.figure_panel.exists():
        raise FileNotFoundError(f"Missing figure panel: {args.figure_panel}")

    figures = pd.read_csv(args.figure_panel, sep="\t")
    _ = figures  # Kept for backward-compatible CLI inputs; the current outline uses a fixed figure plan.

    outline = [
        "# Manuscript Outline v2",
        "",
        "## Working Titles",
        "",
        "1. Bird-biased longevity residuals are associated with sequence-supported transposon-suppression gene retention",
        "2. Transposon suppression links avian lifespan residuals to human repeat-control disease biology",
        "3. A comparative genomic feasibility analysis prioritizes transposon suppression in long-lived birds",
        "",
        "Preferred current title:",
        "",
        "**A sequence-audited transposon/repeat suppression signal in avian-enriched longevity residuals**",
        "",
        "## Current Article Type",
        "",
        "This should be treated as a validation-focused comparative genomics manuscript or preprint-style paper, not yet as a full Nature-family discovery paper. The project now has a 200-gene primary matrix and a targeted sequence-evidence ladder, but high-coverage and formal bird-interaction analyses remain too weak for strong bird-specific or flight-convergence claims.",
        "",
        "## Abstract Draft",
        "",
        "Long-lived birds often exceed mammalian lifespan expectations for their body size, but the genome-maintenance mechanisms associated with this residual longevity remain difficult to test across species because lifespan, phylogeny, and genome annotation quality are strongly confounded. We assembled a comparative vertebrate dataset spanning birds, bats, non-flying mammals, reptiles, and anchor model species, and used body-mass- and clade-adjusted lifespan residuals rather than raw maximum lifespan as the primary phenotype. We then built a 68-species primary genome panel and scored a curated 200-gene genome-maintenance framework covering DNA repair and replication stress, proteostasis/autophagy/mitophagy, inflammation and innate immune restraint, cancer surveillance and senescence, chromatin repression, and transposon/repeat suppression. Across this annotation-aware module screen, transposon/repeat suppression remained a prioritized positive association with longevity residuals, but the signal was sensitive to species coverage and did not support a strong formal bird-specific interaction. To address the main inference risk, we performed targeted sequence validation for 140 high-priority bird species-gene rows from low-coverage transposon/repeat entries. Local assembly/GFF validation, CDS translation, reciprocal DIAMOND/BLAST validation, UniProt sensitivity support, and partial/paralog review placed these rows on an explicit evidence ladder: 62 rows had strict local sequence support, 65 had strict support including external sensitivity evidence, 16 were retained as partial or family-ambiguous but not scoreable, and 45 remained unresolved. Applying the combined validation layer preserved the positive residual-model association but did not make high-coverage or birds-only high-coverage analyses significant. These results support a cautious candidate model in which repeat-control and related chromatin-maintenance genes are prioritized in an avian-enriched longevity-residual framework, while highlighting annotation coverage and orthology validation as central constraints on comparative longevity genomics.",
        "",
        "## Central Claim",
        "",
        "A sequence-audited transposon/repeat suppression module is positively associated with longer-than-expected lifespan in an avian-enriched comparative genome panel, but the association remains coverage-sensitive and should be interpreted as a candidate signal rather than a proven bird-specific mechanism.",
        "",
        "## Claims To Avoid",
        "",
        "- Do not claim universal vertebrate convergence.",
        "- Do not claim a bird-independent effect.",
        "- Do not claim causal human longevity effects.",
        "- Do not claim TRIM28/SETDB1/MOV10 are statistically exceptional within all maintenance genes.",
        "- Do not claim annotation-quality independence until expanded panel and annotation-tier controls are complete.",
        "",
        "## Results Structure",
        "",
        "### Result 1: A curated cross-species lifespan panel identifies model-ready longevity residuals",
        "",
        "Use the Week 1 species and lifespan model results. Emphasize the shift from raw maximum lifespan to body-mass-adjusted residuals, and the need for phylogenetic correction.",
        "",
        "Key evidence:",
        "",
        "- 240 curated species",
        "- 237 model-ready species",
        "- birds, bats, non-flying mammals, reptiles, and anchor species",
        "",
        "### Result 2: Genome availability restricts the primary mechanism panel",
        "",
        "Explain why the primary genome mechanism panel is smaller than the lifespan panel. This is not a failure; it is an annotation-quality decision.",
        "",
        "Key evidence:",
        "",
        "- 68 primary Tier1/Tier2 species",
        "- 96 species without assemblies",
        "- Tier3 assembly-only species held out of primary claims",
        "",
        "### Result 3: A 200-gene genome-maintenance screen prioritizes transposon/repeat suppression as a candidate module",
        "",
        "The current core result is a module-ranking result, not a single seed-panel result. Emphasize that transposon/repeat suppression is prioritized within a broader maintenance landscape that also includes cancer/senescence, inflammation, chromatin repression, proteostasis, and DNA repair.",
        "",
        "Key evidence:",
        "",
        "- 200-gene current primary matrix across 68 genome-panel species",
        "- 236 validated candidate genes in the expanded v2 design",
        "- transposon/repeat module remains positive in all-primary residual and mass/clade models",
        "- candidate signal is not unique enough for causal language",
        "",
        "### Result 4: Coverage and bird-enriched sensitivity analyses define the boundary of the signal",
        "",
        "This section is crucial because it prevents overclaiming.",
        "",
        "Supported/cautionary patterns:",
        "",
        "- all-primary residual model remains positive and nominally supported",
        "- coverage >= 0.25 subset remains positive and nominally supported",
        "- high-coverage subsets at coverage >= 0.50 and >= 0.70 are not significant",
        "- birds-only coverage >= 0.50 is not significant",
        "- explicit bird interaction is not significant",
        "",
        "Interpretation:",
        "",
        "- acceptable wording: avian-enriched candidate signal",
        "- avoid: robust bird-specific effect",
        "- avoid: broad flight convergence",
        "",
        "### Result 5: Targeted sequence validation converts the main vulnerability into an evidence ladder",
        "",
        "This should become one of the most distinctive Results sections because it differentiates the study from ordinary cross-species association screens.",
        "",
        "Key evidence:",
        "",
        "- 140 high-priority bird transposon/repeat species-gene rows",
        "- 91 rows with assembly GFF annotation support",
        "- 58 local GFF protein strict sequence rows",
        "- 4 additional local CDS strict rows",
        "- 3 additional external UniProt sensitivity rows",
        "- 62 total local strict rows",
        "- 65 total strict rows including external sensitivity evidence",
        "- 16 partial/family not-scoreable rows",
        "- 45 unresolved rows",
        "",
        "### Result 6: Combined validation overlay preserves the association but does not establish high-coverage robustness",
        "",
        "Key evidence:",
        "",
        "- strict sequence rows applied to matrix: 7",
        "- mean transposon/repeat coverage: 0.776 to 0.779",
        "- all-module residual model baseline: rank 4, estimate = 0.3178, P = 0.00213",
        "- all-module residual model after validation overlay: rank 4, estimate = 0.3235, P = 0.00204",
        "- transposon coverage >= 0.50: P = 0.1931",
        "- birds-only transposon coverage >= 0.50: P = 0.4171",
        "",
        "Interpretation:",
        "",
        "- sequence validation hardens evidence and preserves the signal",
        "- it does not yet upgrade the claim to coverage-independent robustness",
        "- high-coverage weakness should be presented openly",
        "",
        "### Result 7: Human translation mapping remains supportive but secondary",
        "",
        "Use the existing Week 5 analyses as translational context, not as the main statistical pillar. Human mapping can support biological plausibility for repeat-control, chromatin repression, disease association, and tractability, but it does not prove human aging causality.",
        "",
        "Key evidence:",
        "",
        "- human mapping and disease/druggability files remain useful for discussion",
        "- module enrichment/background/permutation tests are supportive context",
        "- translational claims should remain secondary to the comparative sequence-audited result",
        "",
        "## Figure Plan",
        "",
    ]
    outline.extend(
        [
            "- Fig1 (main): Compact phylogenetic heatmap linking clade, lifespan residual, transposon/repeat score, coverage, and sequence-validation evidence. `results/figures/figure1_phylogenetic_heatmap_polished.png`",
            "- Fig2 (main): Maintenance-module ranking or lifespan-residual module association. `results/figures/phase2_final_module_ranking.png`",
            "- Fig3 (main): Targeted sequence evidence ladder. `results/figures/phase3_sequence_evidence_waterfall.png`",
            "- Fig4 (main or supplement): Combined validation-overlay model impact and coverage boundary. `results/figures/phase3_combined_rescue_model_impact.png`",
            "- FigS1 (supplement): Full 68-species phylogenetic heatmap. `results/figures/figureS1_phylogenetic_heatmap_full.png`",
            "- FigS2 (supplement): Gene-level sequence evidence classes. `results/figures/phase3_evidence_levels_by_gene.png`",
            "- FigS3 (supplement or discussion): Human translational evidence map. `results/figures/week5_translational_evidence_map.png`",
        ]
    )
    outline.extend(
        [
            "",
            "## Discussion Skeleton",
            "",
            "1. Transposon/repeat suppression emerges as a prioritized candidate module in an avian-enriched longevity-residual framework.",
            "2. The strongest contribution is the annotation-aware and sequence-audited workflow, not a causal proof of avian longevity.",
            "3. The evidence ladder reduces the risk that low-coverage bird rows are being interpreted blindly.",
            "4. The result remains biologically plausible because repeat suppression, chromatin repression, DNA damage, senescence, and inflammation are mechanistically connected.",
            "5. Human mapping supports translational relevance, but does not prove human longevity causality.",
            "6. The next validation stage must broaden orthology validation, improve high-coverage bird representation, and upgrade phylogenetic branch lengths.",
            "",
            "## Limitations",
            "",
            "- OpenTree synthetic topology with Grafen fallback branch lengths is acceptable for feasibility but not final inference.",
            "- The primary genome panel has only 68 species.",
            "- The current primary matrix has 200 genes, while the validated design contains 236 genes; final submission should explain this distinction.",
            "- High-coverage and birds-only high-coverage subsets are not significant.",
            "- Formal bird interaction models are not significant.",
            "- Orthology ambiguity remains for 45 unresolved high-priority rows plus 16 partial/family not-scoreable rows.",
            "- Current sequence validation deeply audits 10 priority transposon/repeat genes, not every gene in the 200-gene matrix.",
            "- Human translation evidence is disease/annotation evidence, not ageing causality.",
            "",
            "## Next-Stage Validation",
            "",
            "- Build expanded 200-300 gene panel.",
            "- Re-test repeat/chromatin modules using GO/Reactome gene sets.",
            "- Add annotation-tier and completeness covariates.",
            "- Use dated phylogeny or TimeTree-informed branch lengths.",
            "- Broaden sequence validation beyond the 10 priority transposon/repeat genes.",
            "- Improve high-coverage bird sampling and coverage-balanced bird-only models.",
            "- Use OMA, OrthoDB, Ensembl Compara, UniProt, and local assembly sequence as orthology support layers, with local sequence validation preferred for strict scoring.",
            "- Decide whether to post a conservative preprint before broader next-stage validation.",
            "",
            "## Current Submission Readiness",
            "",
            "- Internal feasibility report: ready.",
            "- Validation-focused preprint or specialized journal manuscript: plausible after polishing figures, Methods, and limitations.",
            "- Nature-family manuscript: not yet; requires high-coverage robustness, stronger phylogeny, broader orthology validation, and clearer independent support for the avian-enriched component.",
        ]
    )
    args.outline.parent.mkdir(parents=True, exist_ok=True)
    args.outline.write_text("\n".join(outline) + "\n", encoding="utf-8")

    executive = [
        "# Executive Summary",
        "",
        "## Decision",
        "",
        "**Cautious-go.** The project has moved from feasibility into validation-focused manuscript preparation, but the claim should remain candidate-level until high-coverage and broader orthology validation improve.",
        "",
        "## Best Current Claim",
        "",
        "A sequence-audited transposon/repeat suppression module is positively associated with longer-than-expected lifespan in an avian-enriched 68-species genome panel, but the association remains coverage-sensitive and should not yet be framed as a proven bird-specific or flight-convergence mechanism.",
        "",
        "## Why It Matters",
        "",
        "The strongest contribution is now the combination of lifespan residual modeling, annotation-aware module scoring, and a targeted sequence-evidence ladder for low-coverage bird transposon/repeat rows.",
        "",
        "## Main Weakness",
        "",
        "The signal remains weak in high-coverage and birds-only high-coverage subsets, and formal bird-interaction models are not significant.",
        "",
        "## Sequence-Validation Evidence",
        "",
        "- High-priority rows: 140",
        "- Strict local sequence support: 62",
        "- Strict support including external UniProt sensitivity: 65",
        "- Partial/family not-scoreable: 16",
        "- Unresolved: 45",
        "",
        "## Next Move",
        "",
        "Polish the validation-focused manuscript package, broaden sequence validation beyond the 10 priority genes, improve high-coverage bird representation, and upgrade the phylogeny before attempting stronger evolutionary claims.",
    ]
    args.executive_summary.parent.mkdir(parents=True, exist_ok=True)
    args.executive_summary.write_text("\n".join(executive) + "\n", encoding="utf-8")
    print(f"Wrote {args.outline} and {args.executive_summary}")


if __name__ == "__main__":
    main()
