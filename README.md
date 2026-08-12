# Gene observability in vertebrate lifespan comparative genomics

Code for a methodologically focused comparative genomics study of gene observability, genome maintenance genes, and vertebrate lifespan evolution. The workflow combines phylogenetic lifespan models, observation-bias simulations, genome-completeness and study-effort controls, sequence-based orthology checks, protein-family trees, and Pfam-domain conservation analyses.

The repository name is historical. The associated manuscript concerns gene observability and comparative genomic inference, not a bird-specific lifespan mechanism.

## Study design

- 68-species primary vertebrate genome panel.
- 200 scoreable genome maintenance genes in six prespecified modules.
- PGLS models adjusted for body mass, clade, genome completeness, and phenotype study effort.
- Matrix-degradation experiments with random and informative observation loss across 500 replicates per condition.
- Six correlated module scores summarized as one shared principal-component axis.
- Sequence classification of 140 low-observability avian species-gene rows.
- Targeted protein and domain analyses, with SAMHD1 retained only as an exploratory supplementary hypothesis.

The analysis does not claim flight convergence, positive selection, altered biochemical activity, a uniquely supported repeat-control mechanism, or six independent pathway confirmations.

## Repository layout

```text
config/       Study definitions, species groups, pathways, and source metadata
docs/         Reproducibility guide and methodological documentation
env/          Lightweight Python and R environment specifications
scripts/      Python and R analysis stages
workflow/     Snakemake workflow
```

Large database downloads, local software installations, processed data, generated results, and journal submission files are excluded from GitHub. Third-party assemblies, annotations, and protein records are referenced by database accession and version.

## Environment

```bash
conda env create -f env/environment.yml
conda activate bird-lifespan
```

Alternatively install the Python packages in `env/requirements.txt` and the R packages listed in `env/R_packages.tsv`. External tools used in sequence follow-up include BLAST+, DIAMOND, MAFFT, IQ-TREE, InterProScan/Pfam, HMMER, samtools, seqkit, and bedtools.

## Reproduction

Run commands from the repository root. The staged commands for the empirical PGLS analyses, degradation experiments, figures, and retained sequence analyses are in [`docs/jme_reproducibility_guide.md`](docs/jme_reproducibility_guide.md).

## Data archive

Processed inputs, figure source data, supplementary tables, alignments, trees, accession manifests, and checksums are archived in Zenodo under concept DOI https://doi.org/10.5281/zenodo.20798436. Use the versioned archive identified by the matching GitHub release.

## Associated publication

This repository contains the analysis code for:

*Gene observability limits comparative genomic inference of vertebrate lifespan evolution.* *Journal of Molecular Evolution*. DOI: pending.

## Maintainer

For questions about the code, please contact the corresponding author listed in the paper.

## Citation

Please cite the associated paper and the versioned archived data DOI.
