"""Write Week 6 phase 3 manuscript package and expanded-panel v2 design."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DEFAULT_OUTDIR_DOCS = pathlib.Path("docs")
DEFAULT_OUTDIR_RESULTS = pathlib.Path("results/reports")
DEFAULT_PANEL = pathlib.Path("config/maintenance_gene_sets_v2_draft.tsv")
DEFAULT_COUNTS = pathlib.Path("results/tables/week6_expanded_panel_v2_module_counts.tsv")


MODULES = {
    "transposon_repeat_suppression": {
        "order": 1,
        "submodules": {
            "piRNA_germline_repeat_control": [
                "PIWIL1",
                "PIWIL2",
                "PIWIL3",
                "PIWIL4",
                "TDRD1",
                "TDRD5",
                "TDRD6",
                "TDRD7",
                "TDRD9",
                "TDRD12",
                "TDRKH",
                "DDX4",
                "MAEL",
                "MOV10L1",
                "HENMT1",
                "PLD6",
                "ASZ1",
                "GTSF1",
                "RNF17",
                "FKBP6",
            ],
            "somatic_retroelement_restriction": [
                "MOV10",
                "APOBEC3A",
                "APOBEC3B",
                "APOBEC3C",
                "APOBEC3F",
                "APOBEC3G",
                "APOBEC3H",
                "SAMHD1",
                "TREX1",
                "ZCCHC3",
                "MORC3",
                "MORC4",
            ],
            "repeat_chromatin_repression": [
                "TRIM28",
                "SETDB1",
                "SETDB2",
                "DNMT1",
                "DNMT3A",
                "DNMT3B",
                "DNMT3L",
                "UHRF1",
                "HELLS",
                "MBD1",
                "MBD2",
                "MBD3",
            ],
        },
    },
    "chromatin_repression_heterochromatin": {
        "order": 2,
        "submodules": {
            "hp1_polycomb_heterochromatin": [
                "CBX1",
                "CBX3",
                "CBX5",
                "CBX2",
                "CBX4",
                "CBX6",
                "CBX7",
                "CBX8",
                "SUV39H1",
                "SUV39H2",
                "EHMT1",
                "EHMT2",
                "EZH1",
                "EZH2",
                "EED",
                "SUZ12",
                "RING1",
                "RNF2",
            ],
            "chromatin_remodeling_and_silencing": [
                "CHD4",
                "CHD3",
                "MTA1",
                "MTA2",
                "HDAC1",
                "HDAC2",
                "SIN3A",
                "SIN3B",
                "SMARCA4",
                "SMARCA5",
                "ATRX",
                "DAXX",
            ],
        },
    },
    "DNA_repair_replication_stress": {
        "order": 3,
        "submodules": {
            "damage_signaling_checkpoint": [
                "ATM",
                "ATR",
                "CHEK1",
                "CHEK2",
                "TP53BP1",
                "MDC1",
                "H2AX",
                "NBN",
                "MRE11",
                "RAD50",
            ],
            "homologous_recombination_fanconi": [
                "BRCA1",
                "BRCA2",
                "PALB2",
                "RAD51",
                "RAD51C",
                "RAD51D",
                "XRCC2",
                "XRCC3",
                "FANCA",
                "FANCC",
                "FANCD2",
                "FANCI",
                "FANCM",
                "BRIP1",
                "BLM",
                "WRN",
            ],
            "end_joining_base_excision_nucleotide_excision": [
                "XRCC5",
                "XRCC6",
                "PRKDC",
                "LIG4",
                "XRCC4",
                "NHEJ1",
                "PARP1",
                "PARP2",
                "XRCC1",
                "LIG3",
                "OGG1",
                "MUTYH",
                "APEX1",
                "ERCC1",
                "ERCC2",
                "ERCC3",
                "ERCC4",
                "ERCC5",
                "XPA",
                "XPC",
            ],
        },
    },
    "proteostasis_autophagy_mitophagy": {
        "order": 4,
        "submodules": {
            "chaperone_heat_shock": [
                "HSPA8",
                "HSPA1A",
                "HSPA1B",
                "HSP90AA1",
                "HSP90AB1",
                "HSPB1",
                "HSF1",
                "DNAJA1",
                "DNAJB1",
                "BAG1",
                "BAG3",
            ],
            "ubiquitin_proteasome": [
                "UBB",
                "UBC",
                "PSMA1",
                "PSMA3",
                "PSMB5",
                "PSMB7",
                "PSMC1",
                "PSMD1",
                "PSMD4",
                "VCP",
                "UBE2D1",
                "UBE2N",
            ],
            "autophagy_mitophagy_mito_dynamics": [
                "BECN1",
                "ATG3",
                "ATG4B",
                "ATG5",
                "ATG7",
                "ATG10",
                "ATG12",
                "ATG16L1",
                "MAP1LC3B",
                "GABARAPL1",
                "SQSTM1",
                "OPTN",
                "TBK1",
                "PINK1",
                "PRKN",
                "MFN1",
                "MFN2",
                "OPA1",
                "DNM1L",
                "BNIP3",
                "BNIP3L",
                "FUNDC1",
                "PARK7",
                "SOD2",
                "GPX1",
            ],
        },
    },
    "inflammation_innate_immune_restraint": {
        "order": 5,
        "submodules": {
            "nfkb_inflammasome_restraint": [
                "NFKB1",
                "NFKBIA",
                "RELA",
                "IKBKB",
                "TNFAIP3",
                "CYLD",
                "NLRP3",
                "PYCARD",
                "CASP1",
                "IL1B",
                "IL6",
                "TNF",
                "SIRT1",
                "SIRT6",
            ],
            "cytosolic_dna_rna_sensing": [
                "MB21D1",
                "STING1",
                "TBK1",
                "IRF3",
                "IFIH1",
                "DDX58",
                "MAVS",
                "ADAR",
                "ADARB1",
                "RNASEH2A",
                "RNASEH2B",
                "RNASEH2C",
            ],
            "interferon_negative_regulation": [
                "SOCS1",
                "SOCS3",
                "USP18",
                "ISG15",
                "TREX1",
                "SAMHD1",
                "NLRC3",
                "NLRX1",
                "PIN1",
                "PPARG",
            ],
        },
    },
    "cancer_surveillance_senescence": {
        "order": 6,
        "submodules": {
            "tumor_suppressor_cell_cycle": [
                "TP53",
                "RB1",
                "PTEN",
                "CDKN2A",
                "CDKN1A",
                "CDKN1B",
                "MDM2",
                "CDKN2B",
                "ATM",
                "CHEK2",
                "GADD45A",
                "FOXO3",
            ],
            "apoptosis_senescence_sasp": [
                "BAX",
                "BAK1",
                "BCL2",
                "BCL2L1",
                "CASP3",
                "CASP8",
                "CASP9",
                "FAS",
                "FASLG",
                "LMNB1",
                "GLB1",
                "SERPINE1",
            ],
            "telomere_and_genome_stability": [
                "TERT",
                "TERC",
                "DKC1",
                "TINF2",
                "POT1",
                "ACD",
                "TERF1",
                "TERF2",
                "RTEL1",
                "RECQL4",
                "STAG2",
                "SMC1A",
                "SMC3",
            ],
        },
    },
}


SEED_GENES = {
    "ATM",
    "ATR",
    "BRCA1",
    "BRCA2",
    "RAD51",
    "XRCC5",
    "XRCC6",
    "PARP1",
    "HSPA8",
    "HSP90AA1",
    "HSF1",
    "UBB",
    "PSMA1",
    "PSMB5",
    "BECN1",
    "ATG5",
    "ATG7",
    "MAP1LC3B",
    "SQSTM1",
    "PINK1",
    "PRKN",
    "MFN1",
    "MFN2",
    "OPA1",
    "BNIP3",
    "NFKB1",
    "RELA",
    "IL6",
    "TNF",
    "NLRP3",
    "SIRT1",
    "PIWIL1",
    "PIWIL2",
    "MOV10",
    "TRIM28",
    "SETDB1",
    "TP53",
    "RB1",
    "PTEN",
    "CDKN2A",
    "CHEK2",
}


def build_panel() -> pd.DataFrame:
    rows = []
    seen: dict[str, str] = {}
    for module, module_data in MODULES.items():
        module_order = module_data["order"]
        gene_order = 0
        for submodule, genes in module_data["submodules"].items():
            for gene in genes:
                if gene in seen:
                    continue
                seen[gene] = module
                gene_order += 1
                rows.append(
                    {
                        "human_gene_symbol": gene,
                        "maintenance_module_v2": module,
                        "submodule_v2": submodule,
                        "module_order": module_order,
                        "gene_order_within_module": gene_order,
                        "seed_status": "seed_v0" if gene in SEED_GENES else "expanded_v2_candidate",
                        "inclusion_tier": "core" if gene in SEED_GENES else "candidate",
                        "source_plan": "curated_draft_to_validate_against_hgnc_reactome_go_uniprot",
                        "orthology_validation_priority": (
                            "high"
                            if module
                            in {
                                "transposon_repeat_suppression",
                                "chromatin_repression_heterochromatin",
                            }
                            else "standard"
                        ),
                    }
                )
    panel = pd.DataFrame(rows)
    if not 200 <= len(panel) <= 300:
        raise ValueError(f"Expanded panel draft has {len(panel)} genes, expected 200-300")
    return panel


def write_docs(panel: pd.DataFrame, args: argparse.Namespace) -> None:
    counts = (
        panel.groupby("maintenance_module_v2", as_index=False)
        .agg(
            n_genes=("human_gene_symbol", "nunique"),
            n_seed=("seed_status", lambda x: int((x == "seed_v0").sum())),
            n_expanded=("seed_status", lambda x: int((x == "expanded_v2_candidate").sum())),
        )
        .sort_values("maintenance_module_v2")
    )
    args.counts.parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(args.counts, sep="\t", index=False)

    narrative = [
        "# Week 6 Phase 3 Manuscript Package",
        "",
        "## Purpose",
        "",
        "This package turns the Week 6 feasibility result into manuscript-facing language while keeping the current claim bounded. It should be used as the bridge between the feasibility sprint and an expanded validation manuscript.",
        "",
        "## Results Narrative Draft",
        "",
        "### Result 1: Comparative lifespan residuals define the phenotype",
        "",
        "We first assembled a 240-species comparative life-history panel spanning birds, bats, non-flying mammals, reptiles, and anchor model species. Because raw maximum lifespan is strongly confounded by body size and clade, the primary phenotype was a body-mass-adjusted lifespan residual rather than lifespan itself. This produced 237 model-ready species for the initial lifespan layer.",
        "",
        "### Result 2: Genome availability narrows the mechanism panel",
        "",
        "The mechanism layer was intentionally restricted to species with usable genome annotation. The primary genome-maintenance panel contains 68 Tier1/Tier2 species, while species without assemblies or with assembly-only records were held out of primary mechanism claims. This reduces sample size but protects the analysis from treating missing annotation as biology.",
        "",
        "### Result 3: Transposon suppression is the strongest current module",
        "",
        "Across the 41-gene seed panel, the strict sequence-supported transposon-suppression score was positively associated with lifespan residuals after body-mass and clade adjustment. The mass+clade PGLS estimate was 0.392 with p = 0.00887 and BH = 0.0345. A weak-inclusive score remained positive and BH-significant.",
        "",
        "### Result 4: Sensitivity analyses define a bird-dependent signal",
        "",
        "The signal survived exclusion of bats, reptiles, non-flying mammals, human, top lifespan-residual outliers, and each individual transposon gene. It did not survive removal of birds, and it was not recovered in the Tier1-only annotation subset. The correct interpretation is therefore bird-dependent and annotation-tier-sensitive, not universal vertebrate convergence.",
        "",
        "### Result 5: Sequence validation substantially hardens the orthology layer",
        "",
        "The final strict transposon score used sequence-supported rows. Most direct NCBI candidates were reciprocal-supported, while unresolved PIWIL2 and low-coverage TRIM28-like cases remain a defined validation queue rather than hidden uncertainty.",
        "",
        "### Result 6: Human mapping supports relevance but not causality",
        "",
        "Human mapping connected focal repeat/chromatin genes to repeat-control evidence, disease context, and selected tractability. TRIM28, SETDB1, and MOV10 are the strongest translation leads, while PIWIL1 and PIWIL2 remain important mechanistic background genes. This layer supports prioritization, not causal human longevity claims.",
        "",
        "## Methods Skeleton",
        "",
        "1. Species curation and life-history harmonization.",
        "2. Body-mass-adjusted lifespan residual modeling.",
        "3. OpenTree synthetic topology and feasibility-stage Grafen branch lengths.",
        "4. Genome availability tiering and primary/sensitivity panel definition.",
        "5. Maintenance gene-set construction and ortholog matrix scoring.",
        "6. Sequence validation for transposon candidates using reciprocal protein evidence.",
        "7. PGLS and sensitivity models for module scores.",
        "8. Human ageing, disease, and tractability mapping.",
        "9. Enrichment, background comparison, and permutation tests.",
        "10. Expanded panel v2 validation plan.",
        "",
        "## Statistical Language To Use",
        "",
        "- Use 'associated with', not 'drives' or 'causes'.",
        "- Use 'longer-than-expected lifespan' for residual-based analyses.",
        "- Use 'bird-dependent' instead of 'avian-universal'.",
        "- Use 'sequence-supported score' rather than 'complete ortholog set'.",
        "- Use 'human disease-context support' rather than 'human ageing mechanism'.",
        "",
        "## Draft Figure Captions",
        "",
        "**Figure 1.** Strict sequence-supported transposon-suppression score versus body-mass-adjusted lifespan residuals across the 68-species primary genome panel. Formal inference is based on PGLS rather than the visual regression alone.",
        "",
        "**Figure 2.** Main-effect and sensitivity forest plot for the strict transposon-suppression score. The association survives most subset and leave-one-gene tests but is lost when birds are removed and in the Tier1-only subset.",
        "",
        "**Figure 3.** Sequence-validation waterfall for transposon-suppression gene-species rows. The plot separates sequence-supported, weak-supported, and unsupported rows and makes remaining orthology ambiguity explicit.",
        "",
        "**Figure 4.** Human translational evidence map for maintenance seed genes. Repeat/chromatin genes show interpretable repeat-control, disease-association, and tractability context, but this evidence is not treated as causal ageing proof.",
        "",
        "## Immediate Manuscript Decision",
        "",
        "A short preprint-style manuscript is possible now, but the stronger route is to complete expanded panel v2 before aiming at a high-profile comparative genomics venue.",
    ]
    phase3 = args.docs / "week6_phase3_manuscript_package.md"
    phase3.write_text("\n".join(narrative) + "\n", encoding="utf-8")

    risk = [
        "# Week 6 Reviewer Risk Table",
        "",
        "| Risk | Why It Matters | Current Answer | Next Validation | Severity |",
        "|---|---|---|---|---|",
        "| 41-gene seed panel may be cherry-picked | A focused panel can inflate module-specific signals | The current result is explicitly framed as feasibility | Build and test 200-300 gene expanded panel v2 | High |",
        "| Signal disappears without birds | Weakens flight-convergence framing | Current claim is bird-dependent, not universal vertebrate convergence | Test birds internally and compare bats only as convergent contrast | High |",
        "| Tier1-only subset loses signal | Annotation quality may affect score estimates | Tier1-only loss is disclosed as a limitation | Add annotation tier/completeness covariates and missingness controls | High |",
        "| OpenTree+Grafen branch lengths are approximate | PGLS inference can depend on branch lengths | Acceptable for feasibility, not final inference | Use dated tree or TimeTree-informed branch lengths | High |",
        "| Orthology ambiguity remains | Repeat/chromatin families have paralog and low-complexity pitfalls | Ambiguous rows are isolated rather than hidden | Validate with OMA/OrthoDB/Compara/domain checks | Medium-High |",
        "| Human disease mapping is not ageing causality | Reviewers may reject translational overreach | Current language says prioritization only | Add independent ageing/disease enrichment controls | Medium |",
        "| Many maintenance modules may correlate with annotation completeness | Could make transposon signal non-specific | Week 5 background tests help but are not enough | Random matched gene sets and module-rank permutation | High |",
        "| Maximum lifespan data are noisy | Captive records and uneven sampling affect residuals | Outlier audits and sensitivity tests were run | Add sample-size/certainty weights if available | Medium |",
    ]
    (args.docs / "week6_reviewer_risk_table.md").write_text(
        "\n".join(risk) + "\n", encoding="utf-8"
    )

    design = [
        "# Expanded Panel v2 Design",
        "",
        "## Decision",
        "",
        "Yes: the 200-300 gene expansion should be designed now and treated as the main validation gate for any manuscript stronger than a feasibility/preprint report.",
        "",
        "## Rationale",
        "",
        "The current 41-gene panel is useful because it identified a coherent transposon-suppression signal, but it is too small and too curated to rule out hand-selection, annotation density, or generic maintenance-gene effects. Expanded panel v2 tests whether repeat/chromatin biology remains prioritized when placed inside a broader, predeclared maintenance universe.",
        "",
        "## Draft Panel",
        "",
        f"Draft file: `{args.panel}`",
        "",
        f"Total draft genes: {panel['human_gene_symbol'].nunique()}",
        "",
        "Module counts:",
        "",
    ]
    for _, row in counts.iterrows():
        design.append(
            f"- {row['maintenance_module_v2']}: {row['n_genes']} genes "
            f"({row['n_seed']} seed, {row['n_expanded']} expanded)"
        )
    design.extend(
        [
            "",
            "## Validation Gates",
            "",
            "1. Symbol validation: confirm HGNC/MyGene-valid symbols and resolve aliases.",
            "2. Source validation: attach GO/Reactome/UniProt/GenAge/CellAge/source-family evidence.",
            "3. Orthology feasibility: estimate species coverage before full scoring.",
            "4. Annotation-bias control: model genome tier, annotation completeness, and per-module missingness.",
            "5. Biological specificity: compare repeat/chromatin rank against DNA repair, proteostasis, mitophagy, inflammation, and cancer/senescence modules.",
            "6. Matched random sets: compare against size- and annotation-richness-matched random gene sets.",
            "",
            "## Primary Test",
            "",
            "Fit PGLS models for each module score against lifespan residuals with body mass and clade covariates. The primary success criterion is that repeat/chromatin modules remain positive and rank near the top after annotation-tier and missingness sensitivity analyses.",
            "",
            "## Failure Interpretation",
            "",
            "If the signal disappears after annotation controls, the correct conclusion is that the 41-gene result was likely an annotation-quality artifact. If many unrelated modules perform equally well, the project should pivot from transposon-specific claims to a broader annotation/maintenance covariation analysis.",
        ]
    )
    (args.docs / "expanded_panel_v2_design.md").write_text(
        "\n".join(design) + "\n", encoding="utf-8"
    )

    report = [
        "# Week 6 Phase 3 Report",
        "",
        "## Outputs",
        "",
        f"- Manuscript package: `{phase3}`",
        f"- Reviewer risk table: `{args.docs / 'week6_reviewer_risk_table.md'}`",
        f"- Expanded panel v2 design: `{args.docs / 'expanded_panel_v2_design.md'}`",
        f"- Expanded panel v2 draft genes: `{args.panel}`",
        f"- Module counts: `{args.counts}`",
        "",
        "## Status",
        "",
        "Week 6 phase 3 is complete as a manuscript-planning and validation-design step. The expanded panel is a curated draft, not yet a validated HGNC/Reactome/GO-backed final panel.",
    ]
    args.results.mkdir(parents=True, exist_ok=True)
    (args.results / "week6_phase3_report.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--docs", type=pathlib.Path, default=DEFAULT_OUTDIR_DOCS)
    parser.add_argument("--results", type=pathlib.Path, default=DEFAULT_OUTDIR_RESULTS)
    parser.add_argument("--panel", type=pathlib.Path, default=DEFAULT_PANEL)
    parser.add_argument("--counts", type=pathlib.Path, default=DEFAULT_COUNTS)
    args = parser.parse_args()

    panel = build_panel()
    args.panel.parent.mkdir(parents=True, exist_ok=True)
    panel.to_csv(args.panel, sep="\t", index=False)
    args.docs.mkdir(parents=True, exist_ok=True)
    write_docs(panel, args)
    print(f"Wrote Week 6 phase 3 package with {len(panel)} draft genes")


if __name__ == "__main__":
    main()
