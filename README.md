# Gene observability in vertebrate lifespan comparative genomics

Code for a methodologically focused comparative genomics study of gene observability, genome maintenance genes, and vertebrate lifespan evolution. The workflow combines phylogenetic lifespan models, gene-observability controls, orthology evidence auditing, targeted protein-family trees, and Pfam-domain conservation analyses.

The repository name is historical. The associated manuscript is currently framed around gene observability and comparative genomic inference, not a claim for a bird-specific lifespan mechanism.

## Study design

- 68-species primary vertebrate genome panel.
- 200 scoreable genome maintenance genes in six prespecified modules.
- Mass- and clade-adjusted PGLS association and sensitivity models.
- Annotation-tier, coverage, matched-random-set, and branch-length controls.
- Sequence auditing of 140 low-observability avian species-gene rows.
- Targeted protein and domain analyses, including an exploratory SAMHD1 result.

The analysis does not claim flight convergence, positive selection, altered biochemical activity, or a uniquely supported repeat-control mechanism.

## Repository layout

```text
config/       Study definitions, species groups, pathways, and source metadata
docs/         Analysis plans, variable definitions, and methodological reports
env/          Lightweight Python/R environment specifications
scripts/      Python and R analysis stages
workflow/     Snakemake workflow
```

Large database downloads, local software installations, processed data, generated results, and journal submission files are intentionally excluded from GitHub. Analysis-ready data and frozen outputs will be deposited in a DOI-bearing archive.

## Environment

```bash
conda env create -f env/environment.yml
conda activate bird-lifespan
```

Alternatively install the Python packages in `env/requirements.txt` and the R packages listed in `env/R_packages.tsv`. External tools used in sequence-validation stages include BLAST+, DIAMOND, MAFFT, IQ-TREE, InterProScan/Pfam, HMMER, samtools, seqkit, and bedtools.

## Workflow

Run commands from the repository root. The workflow entry point is:

```bash
snakemake --snakefile workflow/Snakefile --cores 4
```

Individual stages can also be run from the numbered directories under `scripts/`. Upstream database acquisition requires the source releases and accession manifests described in `config/` and the project documentation.

## Reproducibility and data availability

The public release is split into two parts:

1. This GitHub repository contains code, configuration, documentation, and environment specifications.
2. The linked Zenodo archive contains processed inputs, figure source data, supplementary tables, alignments, trees, accession manifests, and checksums.

Third-party assemblies, annotations, and protein records are referenced by database accession and version rather than redistributed when they remain publicly retrievable.

## Data archive

Processed data and frozen outputs are archived in Zenodo: https://doi.org/10.5281/zenodo.21063039

## Associated publication

This repository contains the analysis code for:

*Gene observability limits comparative genomic inference of vertebrate lifespan evolution.* *Journal of Molecular Evolution*. DOI: pending.

## Maintainer

For questions about the code, please contact the corresponding author listed in the paper.

## Citation

Please cite the associated paper and the archived data DOI.
