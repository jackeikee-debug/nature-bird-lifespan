# Variable Dictionary

| variable | description |
| --- | --- |
| `scientific_name` | Binomial or trinomial species name used for joining data sources. |
| `anage_matched_name` | Name used to join AnAge after applying taxonomy aliases. |
| `common_name` | Common name from AnAge where available. |
| `clade` | Project-level clade label. |
| `flight_status` | `flighted`, `flightless`, `powered_flight`, `non_flying`, or `gliding`. |
| `body_mass_g` | Adult body mass in grams, normalized from AnAge. |
| `max_lifespan_years` | Maximum reported lifespan in years from AnAge. |
| `sexual_maturity_years` | Sexual maturity age in years where recoverable from AnAge. |
| `habitat` | Coarse ecological habitat label. |
| `diet` | Coarse diet label. |
| `genome_available` | Whether a genome assembly is expected to be available. |
| `annotation_quality` | Working genome annotation confidence: `high`, `medium`, or `low`. |
| `lifespan_residual_log10` | Observed minus predicted log10 lifespan after body-mass correction. |
| `lifespan_residual_ratio` | Fold change implied by the log10 residual. |
| `life_history_data_source` | Source layer used for mass, lifespan, and maturity values. |
| `manual_supplement_source` | URL or citation note for manually supplemented values. |
| `taxonomy_issue` | Audit classification for missing aliases, missing AnAge records, and incomplete fields. |
| `outlier_flags` | Semicolon-separated residual, lifespan, and mass audit flags. |
| `audit_priority` | Whether a row should receive manual review before downstream modeling. |
| `recommended_use` | How a replacement candidate should be interpreted. |
| `tree_label` | Underscore-separated species label for tree tip matching. |
| `tree_search_name` | Scientific name submitted to tree or taxonomy matching tools. |
| `source_hint` | Suggested phylogeny source family for the species. |
| `tree_match_status` | Manual or automated status after matching against an external tree. |
| `ott_id` | OpenTree taxonomy identifier, when available. |
