# Phylogeny and PGLS Prep

This stage prepares the 237 strict model-ready species for tree matching and
PGLS. It does not assume that a single external tree source will cover all
species cleanly.

## Outputs

- `data/processed/pgls_species.tsv`: model-ready phenotype and residual table
  with tree labels.
- `data/processed/tree_label_audit.tsv`: taxonomy and tree-match review table.
- `data/processed/phylogeny_inputs/opentree_tnrs_names.txt`: all species names
  for OpenTree TNRS matching.
- `data/processed/phylogeny_inputs/birdtree_species.txt`: bird-only list.
- `data/processed/phylogeny_inputs/vertlife_mammal_species.txt`: mammal and bat
  list.
- `data/processed/phylogeny_inputs/reptile_species.txt`: reptile list.
- `results/reports/phylogeny_prep_report.md`: summary counts and risk flags.
- Optional after network matching: `data/processed/opentree_tnrs_matches.tsv`
  and `data/processed/tree_label_audit_opentree.tsv`.
- `data/processed/pgls_tree_ready.tsv`: phenotype, residual, OpenTree ID, and
  final tree-label table.
- `data/processed/pgls_tree_manual_review.tsv`: rows that need manual review
  before final tree pruning.

## Tree Source Strategy

Use multiple sources rather than forcing all species into one tree immediately:

- Birds: BirdTree or OpenTree.
- Mammals and bats: VertLife or OpenTree.
- Reptiles: OpenTree or a published reptile supertree.

After matching, prune to the intersection of available tree tips and
`pgls_species.tsv`, then fit OLS and PGLS on the same species subset.

## Manual Review Rules

Review rows flagged as:

- `trinomial_or_subspecies`
- `anage_alias_or_subspecies`
- `anage_quality_low`
- `anage_quality_questionable`
- `sample_size_tiny`
- `sample_size_small`
- `missing_sexual_maturity`

The tree-match audit table contains empty columns for `ott_id`,
`matched_tree_tip`, and `manual_tree_label`. These should be filled only after
checking the actual tree or TNRS result.

## OpenTree TNRS

Run this after `prepare_tree_labels.py` when network access is available:

```powershell
python scripts\02_phylogeny\match_opentree_tnrs.py
```

The script writes OpenTree taxonomy identifiers and an updated tree-label audit
table, but PGLS should still be run only after checking that the final tree tips
match these labels.

Then join the TNRS results to the phenotype table:

```powershell
python scripts\02_phylogeny\build_pgls_tree_ready.py
```

Fetch an induced OpenTree subtree and build the matching trait table:

```powershell
python scripts\02_phylogeny\fetch_opentree_subtree.py
```

This writes:

- `data/processed/phylogeny_inputs/opentree_induced_subtree.tre`
- `data/processed/pgls_trait_table.tsv`
- `data/processed/opentree_subtree_audit.tsv`

Current OpenTree result:

- 237 OTT IDs requested.
- 233 tips recovered into the species-level trait table.
- 4 taxa were reported in the OpenTree `broken` field and should be excluded
  from species-level PGLS unless manually resolved:
  - `Apteryx australis`
  - `Tyto alba`
  - `Miniopterus schreibersii`
  - `Myotis lucifugus`

## First OLS Baseline

The script below fits OLS on the exact 233-tip subset that will enter PGLS:

```powershell
python scripts\03_lifespan_models\fit_ols_tree_subset.py
```

Outputs:

- `results/tables/tree_subset_ols_models.tsv`
- `data/processed/tree_subset_ols_residuals.tsv`
- `results/reports/tree_subset_ols_report.md`

## PGLS Template

`scripts/03_lifespan_models/fit_pgls_first_pass.R` is a first-pass R template.
It requires `ape`, `caper`, and `nlme`. The OpenTree synthetic tree is
topology-only, so the template applies Grafen branch lengths with
`ape::compute.brlen`. This is acceptable for feasibility screening but should
be replaced by better branch lengths before manuscript-scale inference.

Current first-pass outputs:

- `results/tables/pgls_first_pass_models.tsv`
- `results/models/pgls_first_pass_summary.txt`
- `data/processed/pgls_first_pass_residuals.tsv`

Current model set:

```text
Model A: log10(max_lifespan_years) ~ log10(body_mass_g)
Model B: log10(max_lifespan_years) ~ log10(body_mass_g) + flight_status
Model C: log10(max_lifespan_years) ~ log10(body_mass_g) + clade
```

Do not fit `flight_status + clade` together at this stage: powered flight and
bat clade membership are structurally confounded in the current panel, causing
a singular design matrix.
