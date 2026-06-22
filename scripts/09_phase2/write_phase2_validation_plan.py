"""Write Phase 2 expanded validation planning documents."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DEFAULT_PANEL = pathlib.Path("config/maintenance_gene_sets_v2_draft.tsv")
DEFAULT_DOCS = pathlib.Path("docs")
DEFAULT_REPORTS = pathlib.Path("results/reports")
DEFAULT_TABLES = pathlib.Path("results/tables")


PHASE2_TASKS = [
    {
        "stage": "P2.0",
        "name": "literature_novelty_surveillance",
        "objective": "Track directly overlapping avian longevity, bird-bat flight convergence, and transposon-ageing comparative studies before committing heavy compute.",
        "primary_output": "docs/phase2_novelty_scan.md",
        "decision_gate": "Project framing remains distinct from published bird/bat longevity studies and any new paper is triaged before analysis expansion.",
        "status": "initialized",
    },
    {
        "stage": "P2.1",
        "name": "gene_symbol_source_validation",
        "objective": "Validate HGNC/MyGene symbols and attach source evidence for the 236-gene draft panel.",
        "primary_output": "data/processed/maintenance_gene_sets_v2_validated.tsv",
        "decision_gate": "At least 90 percent of genes have validated current symbols and at least one source evidence tag.",
        "status": "not_started",
    },
    {
        "stage": "P2.2",
        "name": "orthology_feasibility_audit",
        "objective": "Estimate gene-by-species coverage and flag paralog/low-confidence risk before full scoring.",
        "primary_output": "data/processed/orthology_feasibility_v2.tsv",
        "decision_gate": "Each primary module has enough coverage for stable scoring; high-risk genes are moved to sensitivity tier.",
        "status": "not_started",
    },
    {
        "stage": "P2.3",
        "name": "annotation_bias_model_design",
        "objective": "Define annotation tier, completeness, and missingness covariates before rerunning lifespan models.",
        "primary_output": "results/reports/phase2_annotation_bias_design.md",
        "decision_gate": "Models explicitly separate biological score effects from genome annotation completeness.",
        "status": "not_started",
    },
    {
        "stage": "P2.4",
        "name": "expanded_module_scoring",
        "objective": "Build coverage-aware and confidence-weighted module scores for the validated v2 panel.",
        "primary_output": "data/processed/maintenance_scores_v2.tsv",
        "decision_gate": "Module scores have documented missingness and confidence weights.",
        "status": "not_started",
    },
    {
        "stage": "P2.5",
        "name": "expanded_pgls_and_sensitivity",
        "objective": "Rerun PGLS for v2 module scores with bird, tier, outlier, and leave-submodule-out sensitivity.",
        "primary_output": "results/tables/phase2_expanded_pgls_summary.tsv",
        "decision_gate": "Repeat/chromatin signal remains positive and ranks highly after bias controls.",
        "status": "not_started",
    },
    {
        "stage": "P2.6",
        "name": "matched_random_gene_set_tests",
        "objective": "Test whether repeat/chromatin modules outperform size- and annotation-matched random gene sets.",
        "primary_output": "results/tables/phase2_matched_random_set_tests.tsv",
        "decision_gate": "Observed repeat/chromatin score is stronger than matched random expectation.",
        "status": "not_started",
    },
    {
        "stage": "P2.7",
        "name": "phylogeny_upgrade_sensitivity",
        "objective": "Compare OpenTree+Grafen feasibility results with dated or TimeTree-informed branch lengths.",
        "primary_output": "results/tables/phase2_phylogeny_sensitivity.tsv",
        "decision_gate": "Main inference is not an artifact of feasibility-stage branch lengths.",
        "status": "not_started",
    },
    {
        "stage": "P2.8",
        "name": "manuscript_decision_gate",
        "objective": "Decide whether the project is ready for full manuscript drafting, preprint, or pivot.",
        "primary_output": "results/reports/phase2_final_decision_report.md",
        "decision_gate": "Go, cautious-go, pivot, or stop decision is made from predeclared evidence.",
        "status": "not_started",
    },
]


def module_counts(panel: pd.DataFrame) -> pd.DataFrame:
    return (
        panel.groupby("maintenance_module_v2", as_index=False)
        .agg(
            n_genes=("human_gene_symbol", "nunique"),
            n_seed=("seed_status", lambda x: int((x == "seed_v0").sum())),
            n_expanded=("seed_status", lambda x: int((x == "expanded_v2_candidate").sum())),
            high_orthology_priority=(
                "orthology_validation_priority",
                lambda x: int((x == "high").sum()),
            ),
        )
        .sort_values("maintenance_module_v2")
    )


def write_phase2_plan(panel: pd.DataFrame, counts: pd.DataFrame, docs: pathlib.Path) -> None:
    lines = [
        "# Phase 2 Expanded Validation Plan",
        "",
        "## Purpose",
        "",
        "Phase 2 converts the six-week feasibility result into a predeclared validation program. The goal is to test whether the Week 4 transposon/repeat-suppression signal survives a broader gene universe, stronger orthology checks, annotation-bias controls, and improved phylogenetic sensitivity.",
        "",
        "## Starting Point",
        "",
        "- Week 6 decision: cautious-go.",
        "- Current best signal: bird-dependent, sequence-supported transposon-suppression association with lifespan residuals.",
        "- Current weakness: 41-gene seed panel, Tier1-only sensitivity failure, bird dependence, approximate OpenTree+Grafen branch lengths.",
        f"- Expanded draft panel: {panel['human_gene_symbol'].nunique()} unique genes.",
        "",
        "## Module Structure",
        "",
    ]
    for _, row in counts.iterrows():
        lines.append(
            f"- {row['maintenance_module_v2']}: {row['n_genes']} genes "
            f"({row['n_seed']} seed, {row['n_expanded']} expanded, "
            f"{row['high_orthology_priority']} high-priority orthology checks)"
        )
    lines.extend(
        [
            "",
            "## Core Principle",
            "",
            "Phase 2 should behave like a pre-registration layer for the project. New results can change the interpretation, but the main tests, success criteria, and failure criteria should not be rewritten after seeing the answers.",
            "",
            "## Work Packages",
            "",
        ]
    )
    for task in PHASE2_TASKS:
        lines.extend(
            [
                f"### {task['stage']}: {task['name']}",
                "",
                f"Objective: {task['objective']}",
                "",
                f"Primary output: `{task['primary_output']}`",
                "",
                f"Decision gate: {task['decision_gate']}",
                "",
            ]
        )
    lines.extend(
        [
            "## Go Criteria",
            "",
            "- Repeat/chromatin modules remain positive after annotation-tier and missingness controls.",
            "- Repeat/chromatin modules rank above generic maintenance modules and matched random gene sets.",
            "- Bird-internal analyses support that the signal is not only a cross-clade annotation artifact.",
            "- Orthology ambiguity is documented and does not dominate the strict score.",
            "- Phylogeny sensitivity does not reverse the main conclusion.",
            "",
            "## Pivot Criteria",
            "",
            "- Multiple unrelated modules perform similarly, suggesting broad annotation/maintenance covariation rather than transposon specificity.",
            "- The effect becomes a bird-only descriptive pattern with insufficient orthology support for mechanism claims.",
            "- Human mapping remains useful only as prioritization, not as a manuscript-driving result.",
            "",
            "## Stop Criteria",
            "",
            "- The repeat/chromatin signal disappears after annotation controls.",
            "- Orthology coverage is too sparse or too ambiguous to score repeat/chromatin modules reliably.",
            "- The result depends on a small number of low-confidence species or genes.",
            "",
            "## Expected Manuscript Outcomes",
            "",
            "- Strong go: expanded repeat/chromatin signal survives controls; write full comparative genomics manuscript.",
            "- Cautious go: signal survives partially; write preprint/short report with clear limitations.",
            "- Pivot: broaden story to genome-maintenance/annotation-aware longevity patterns.",
            "- Stop: archive as feasibility result and avoid overclaiming.",
        ]
    )
    (docs / "phase2_expanded_validation_plan.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def write_decision_gates(docs: pathlib.Path) -> None:
    lines = [
        "# Phase 2 Decision Gates",
        "",
        "| Gate | Pass | Caution | Fail | Action |",
        "|---|---|---|---|---|",
        "| Literature novelty | No direct duplicate of residual-based repeat/chromatin validation | Adjacent bird/bat longevity papers overlap in framing | A paper already tests the same score/model/question | Reframe before compute expansion |",
        "| Gene symbol/source validation | >=90% validated symbols and source tags | 75-90% validated | <75% validated | Fix aliases or shrink panel before orthology |",
        "| Orthology feasibility | All primary modules scorable with documented confidence | One or two modules sparse | Repeat/chromatin sparse or paralog-dominated | Move weak genes to sensitivity or redesign panel |",
        "| Annotation bias | Repeat/chromatin effect survives tier/missingness controls | Effect weakens but remains directional | Effect disappears | Reframe as annotation artifact or stop |",
        "| Expanded PGLS | Repeat/chromatin positive and high-ranked | Positive but not top-ranked | Null or reversed | Pivot or stop transposon-specific claim |",
        "| Matched random sets | Observed module exceeds matched random expectation | Borderline percentile | Not above random expectation | Avoid specificity claim |",
        "| Phylogeny sensitivity | Same conclusion under upgraded branch lengths | Effect size shifts but direction stable | Conclusion reverses | Treat feasibility result as tree-sensitive |",
        "| Manuscript readiness | Claims survive all high-risk gates | One high-risk limitation remains | Multiple high-risk gates fail | Draft preprint, pivot, or archive |",
        "",
        "## Non-Negotiable Rule",
        "",
        "Do not upgrade the claim from feasibility to discovery unless the expanded panel, annotation-bias controls, and orthology checks all support the same qualitative conclusion.",
    ]
    (docs / "phase2_decision_gates.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def write_task_table(tables: pathlib.Path) -> None:
    tasks = pd.DataFrame(PHASE2_TASKS)
    tasks.to_csv(tables / "phase2_task_register.tsv", sep="\t", index=False)


def write_novelty_scan(docs: pathlib.Path) -> None:
    lines = [
        "# Phase 2 Novelty Scan",
        "",
        "Date: 2026-06-11",
        "",
        "## Bottom Line",
        "",
        "The project is not working in an empty field. There are already strong adjacent studies on avian longevity, bird-bat flight/longevity convergence, mammalian longevity genomics, and transposable elements in ageing. However, the current project remains distinct if Phase 2 focuses on residual-based lifespan modeling, sequence-supported repeat/chromatin module scores, annotation-bias controls, and predeclared expanded validation.",
        "",
        "## Closest Adjacent Studies",
        "",
        "1. Martini et al. 2025, Aging Cell: Avian Lifespan Network Reveals Shared Mechanisms and New Key Players in Animal Longevity. This study analyzed 141 bird genomes and built an avian lifespan network, emphasizing metabolism and cell-cycle control.",
        "   - DOI: https://doi.org/10.1111/acel.70156",
        "   - PubMed: https://pubmed.ncbi.nlm.nih.gov/40580161/",
        "",
        "2. Matsuda and Makino 2024, Proceedings B: Comparative genomics reveals convergent signals associated with the high metabolism and longevity in birds and bats. This is the nearest conceptual competitor for the flight/longevity framing; it tested evolutionary rate shifts and convergent amino-acid substitutions in flying species.",
        "   - DOI: https://doi.org/10.1098/rspb.2024.1068",
        "   - PubMed: https://pubmed.ncbi.nlm.nih.gov/39191281/",
        "",
        "3. Li et al. 2021, Molecular Biology and Evolution: Comparative Analysis of Mammal Genomes Unveils Key Genomic Variability for Human Life Span. This is a mammal-focused comparative longevity genomics paper identifying many longevity-associated amino-acid changes and pathways.",
        "   - DOI: https://doi.org/10.1093/molbev/msab219",
        "",
        "4. Guio/Vieira review 2025, Genome Biology and Evolution: Exploring the Relationship of Transposable Elements and Ageing. This supports the biological plausibility of TE/ageing links but is a review, not an avian/bird-bat residual comparative score study.",
        "   - DOI: https://doi.org/10.1093/gbe/evaf088",
        "",
        "## Current Differentiation",
        "",
        "- Uses lifespan residuals rather than raw maximum lifespan or only long/short categories.",
        "- Focuses on sequence-supported transposon/repeat suppression and chromatin repression as testable modules.",
        "- Explicitly models genome availability and annotation-tier sensitivity as risks rather than hiding them.",
        "- Builds a predeclared 236-gene Phase 2 validation panel rather than interpreting only a 41-gene seed result.",
        "- Connects comparative signal to human ageing/disease/tractability evidence, while avoiding causal human longevity claims.",
        "",
        "## Main Risk",
        "",
        "The 2024 bird-bat paper is close enough that a generic 'flight explains longevity through genome maintenance' manuscript would look derivative. The project must therefore avoid a broad flight-convergence claim unless Phase 2 produces new evidence. The stronger niche is repeat/chromatin maintenance as an annotation-aware, sequence-supported module associated with bird-biased lifespan residuals.",
        "",
        "## Phase 2 Rule",
        "",
        "Before heavy orthology and PGLS expansion, rerun this novelty scan with the exact keywords: avian longevity network, bird bat longevity comparative genomics, transposon suppression longevity birds, repeat silencing lifespan residuals, and chromatin repression avian lifespan.",
    ]
    (docs / "phase2_novelty_scan.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def write_status_report(
    panel: pd.DataFrame, counts: pd.DataFrame, reports: pathlib.Path
) -> None:
    high_priority = int((panel["orthology_validation_priority"] == "high").sum())
    lines = [
        "# Phase 2 Status Report",
        "",
        "## Current Status",
        "",
        "Phase 2 has been initialized as a structured expanded-validation program. No validation results have been run yet; this report records the starting state and prevents the next phase from drifting into post-hoc interpretation.",
        "",
        "## Starting Assets",
        "",
        f"- Expanded panel draft genes: {panel['human_gene_symbol'].nunique()}",
        f"- Modules: {panel['maintenance_module_v2'].nunique()}",
        f"- High-priority orthology genes: {high_priority}",
        "- Main current claim: bird-dependent transposon/repeat-suppression signal in lifespan residuals.",
        "- Main current limitation: 41-gene seed panel and annotation-tier sensitivity.",
        "- Literature status: adjacent bird/bat longevity and avian lifespan-network studies exist; no direct duplicate has been identified for the Phase 2 residual-based repeat/chromatin validation framing.",
        "",
        "## Module Counts",
        "",
    ]
    for _, row in counts.iterrows():
        lines.append(f"- {row['maintenance_module_v2']}: {row['n_genes']} genes")
    lines.extend(
        [
            "",
            "## Next Immediate Step",
            "",
            "Run P2.1 gene symbol/source validation for the 236-gene draft panel, producing `data/processed/maintenance_gene_sets_v2_validated.tsv`.",
        ]
    )
    (reports / "phase2_status_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", type=pathlib.Path, default=DEFAULT_PANEL)
    parser.add_argument("--docs", type=pathlib.Path, default=DEFAULT_DOCS)
    parser.add_argument("--reports", type=pathlib.Path, default=DEFAULT_REPORTS)
    parser.add_argument("--tables", type=pathlib.Path, default=DEFAULT_TABLES)
    args = parser.parse_args()

    if not args.panel.exists():
        raise FileNotFoundError(f"Missing expanded panel draft: {args.panel}")

    args.docs.mkdir(parents=True, exist_ok=True)
    args.reports.mkdir(parents=True, exist_ok=True)
    args.tables.mkdir(parents=True, exist_ok=True)

    panel = pd.read_csv(args.panel, sep="\t")
    counts = module_counts(panel)
    write_phase2_plan(panel, counts, args.docs)
    write_decision_gates(args.docs)
    write_novelty_scan(args.docs)
    write_task_table(args.tables)
    write_status_report(panel, counts, args.reports)

    print(
        "Wrote Phase 2 validation plan for "
        f"{panel['human_gene_symbol'].nunique()} genes"
    )


if __name__ == "__main__":
    main()
