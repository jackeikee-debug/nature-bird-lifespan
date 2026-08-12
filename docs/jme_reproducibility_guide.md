# Journal of Molecular Evolution reproducibility guide

## Scope

This guide rebuilds the analyses added for the second major revision, followed by the retained sequence and sensitivity analyses. Database downloads and assembly rescue are upstream acquisition stages and are not repeated here.

Run all commands from the project root. The matching versioned Zenodo archive contains the processed inputs required below.

## Environment

- Python dependencies: `env/environment.yml` or `env/requirements.txt`
- R packages: `env/R_packages.tsv`
- Workflow definition: `workflow/Snakefile`
- R version used for the submitted analyses: 4.5.2

## Major revision analyses

```powershell
python scripts/13_major_revision2/build_study_effort_covariates.py
python scripts/13_major_revision2/build_revision2_analysis_table.py
Rscript scripts/13_major_revision2/run_effort_busco_pgls.R
python scripts/13_major_revision2/build_degradation_inputs.py
Rscript scripts/13_major_revision2/run_annotation_degradation_simulation.R
python scripts/13_major_revision2/plot_revision2_core_results.py
python scripts/13_major_revision2/build_revision2_main_figures.py
Rscript scripts/13_major_revision2/plot_circular_phylogenetic_heatmap.R
python scripts/13_major_revision2/plot_figure5_degradation_sequence.py
```

These commands build the study-effort covariates, merge BUSCO and assembly data, fit the revised PGLS models, run 500-replicate matrix-degradation experiments, and rebuild Figures 1-5. The degradation analysis is an intervention on the observed binary matrix, not physical assembly fragmentation and reannotation.

## Retained sequence and sensitivity analyses

```powershell
python scripts/12_revision/build_jme_revision3_sensitivities.py
Rscript scripts/12_revision/run_module_weight_sensitivity_pgls.R
Rscript scripts/12_revision/run_samhd1_alignment_sensitivity_pgls.R
python scripts/12_revision/build_jme_figures_and_tables.py
```

These commands reproduce the matched random gene sets, evidence weighting sensitivity, branch-length sensitivity, sequence classification, protein-family trees, and supplementary protein and domain analyses.

## Principal inputs

- `data/processed/ortholog_matrix_primary_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv`
- `data/processed/maintenance_lifespan_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv`
- `data/processed/phase2_strict_v2_scoring_eligibility_sequence_updated.tsv`
- `data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre`
- `data/interim/protein_conservation/SAMHD1.aligned.faa`
- `results/tables/phase2_W3_full_background_matched_random_set_null.tsv`

## Principal revision outputs

- `data/processed/jme_revision2_analysis_table.tsv`
- `results/tables/jme_revision2_effort_busco_pgls.tsv`
- `results/tables/jme_revision2_module_correlations.tsv`
- `results/tables/jme_revision2_degradation_replicates.tsv`
- `results/tables/jme_revision2_degradation_summary.tsv`
- `results/figures/jme_revision2_figure1_circular_phylogeny.*`
- `results/figures/jme_revision2_figure1_study_design.*` (submitted as Figure 2)
- `results/figures/jme_revision2_figure2_shared_module_axis.*` (submitted as Figure 3)
- `results/figures/jme_revision2_effort_busco_forest.*` (submitted as Figure 4)
- `results/figures/jme_revision2_figure5_degradation_sequence.*`

## Interpretation boundaries

The six module scores are strongly correlated and are summarized by a shared first principal component. Their concordance is not interpreted as six independent biological confirmations. BUSCO complete-case analyses are sensitivity analyses because BUSCO values were available for 49 of 68 species and were calculated with heterogeneous vertebrate lineage datasets. Study-effort proxies are observational controls, not direct measures of lifespan sampling intensity.
