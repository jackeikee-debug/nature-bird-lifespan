# Data Gap Resolution Notes

Generated during Week 1 audit after building:

- `data/processed/taxonomy_audit.tsv`
- `data/processed/outlier_audit.tsv`
- `data/processed/replacement_candidates.tsv`

## Current Missing Model Fields

Three species remain incomplete for the lifespan residual model:

| species | current issue | recommended default |
| --- | --- | --- |
| `Aptenodytes forsteri` | AnAge has body mass and maturity but no maximum longevity | keep out of model unless manual supplement is approved |
| `Gekko japonicus` | no same-species AnAge record | replace with AnAge-backed gecko if reptile sample needs this slot |
| `Python bivittatus` | no same-species AnAge record | replace with AnAge-backed pythonid; do not alias by default |

## Supplement Policy

Manual supplements live in `config/manual_life_history_supplements.tsv`.

A supplement enters `species_master.tsv` only when:

1. `use_in_model` is set to `yes`.
2. The filled value is in the same unit as the target field.
3. A source name or URL is recorded.
4. The note explains why the value is comparable to the AnAge field.

Secondary species accounts are allowed as candidates but should not be promoted
to model values until manually reviewed.

## Candidate Supplement

`Aptenodytes forsteri` has a candidate maximum lifespan value from Animal
Diversity Web:

- Candidate maximum lifespan: 50 years
- Source: https://animaldiversity.org/accounts/Aptenodytes_forsteri/
- Status: staged in `manual_life_history_supplements.tsv` with `use_in_model=no`

This value is not currently used in the model.

## Replacement Options

If a species is retained for genome-level anchoring but lacks lifespan data, do
not impute lifespan from related species. Instead, replace the species in the
life-history panel while keeping the original species in a genome-anchor list if
needed.

Top replacement candidates from AnAge:

| target | preferred candidate | reason |
| --- | --- | --- |
| `Aptenodytes forsteri` | `Aptenodytes patagonicus` | same genus, model-ready AnAge record |
| `Gekko japonicus` | `Gekko gecko` | same genus, model-ready AnAge record |
| `Python bivittatus` | `Python molurus` | same genus and historical taxonomy candidate, but not a default alias |

See `data/processed/replacement_candidates.tsv` for additional same-family
options.

## Recommendation

For the next PGLS-ready panel, use strict species-level AnAge values by default:

1. Keep the current 237 model-ready species for residual modeling.
2. Keep `Aptenodytes forsteri`, `Gekko japonicus`, and `Python bivittatus` in
   the project as genome or biological-interest anchors, but exclude them from
   lifespan residual modeling unless reviewed supplements are approved.
3. If exact group balance is required, replace the missing species using
   `replacement_candidates.tsv` rather than imputing values.

