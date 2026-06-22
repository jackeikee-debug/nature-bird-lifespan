# Journal of Molecular Evolution reproducibility guide

## Scope

This guide reproduces the final analysis tables, sensitivity analyses, main figures, and supplementary figures from the processed inputs included in the linked data archive. Database downloads and assembly rescue are upstream acquisition stages and are not repeated by the final analysis rebuild.

## Environment

- Python dependencies: `env/environment.yml` or `env/requirements.txt`
- R library: `env/R_library` when present, otherwise the active R library
- Workflow definition: `workflow/Snakefile`
- Expected R version: 4.5.2

Run commands from the project root.

## Final analysis rebuild

```powershell
python scripts/12_revision/build_jme_revision3_sensitivities.py
Rscript scripts/12_revision/run_module_weight_sensitivity_pgls.R
Rscript scripts/12_revision/run_samhd1_alignment_sensitivity_pgls.R
python scripts/12_revision/build_jme_figures_and_tables.py
```

These commands perform the following stages in order:

1. Audits matching quality for the 5,000 retained random gene sets.
2. Recalculates module scores under three evidence-confidence weight schemes.
3. Fits DateLife-tree PGLS weight-sensitivity models in R.
4. Audits the SAMHD1 alignment by species and human-reference position.
5. Repeats SAMHD1 PGLS after removal of gap-heavy alignment columns.
6. Rebuilds Figures 1-4, Supplementary Figure 7, Table 1, and figure-source mappings.

## Principal inputs

- `data/processed/ortholog_matrix_primary_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv`
- `data/processed/maintenance_lifespan_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv`
- `data/processed/phase2_strict_v2_scoring_eligibility_sequence_updated.tsv`
- `data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre`
- `data/interim/protein_conservation/SAMHD1.aligned.faa`
- `results/tables/phase2_W3_full_background_matched_random_set_null.tsv`

## Principal outputs

- `results/tables/jme_matched_random_gene_set_audit.tsv`
- `results/tables/module_weight_sensitivity_pgls.tsv`
- `results/tables/samhd1_alignment_position_qc.tsv`
- `results/tables/samhd1_alignment_species_qc.tsv`
- `results/tables/samhd1_alignment_sensitivity_pgls.tsv`
- `results/figures/jme_figure1_study_design.*`
- `results/figures/jme_figure2_module_forest.*`
- `results/figures/jme_figure3_orthology_ladder.*`
- `results/figures/jme_figure4_samhd1_domain_evolution.*`
- `results/figures/jme_figureS7_samhd1_clade_heatmap.*`

## Random-set matching definition

For each target module, candidate sets are drawn without replacement within a set from the 139 final maintenance genes outside the chromatin and transposon/repeat target modules. The closest 5,000 of 150,000 candidate draws are retained by absolute difference in mean gene observability. Matching does not include gene length or sequence conservation. The retained null evaluates non-phylogenetic residual correlation and residual slope, not a PGLS null.

## Manual database stages

NCBI, UniProt, OMA, OrthoDB, Ensembl, GFF/CDS rescue, DIAMOND/BLAST, InterProScan, MAFFT, and IQ-TREE steps require the database versions and local assets recorded in the row-level evidence tables and reports. Their cached outputs are retained so the final analysis rebuild does not depend on live network services.
