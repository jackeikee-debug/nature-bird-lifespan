#!/usr/bin/env python3
"""Build prespecified scoring strategies for partial/family-ambiguous rows."""

from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path("data/processed/maintenance_lifespan_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv")
EVIDENCE = Path("data/processed/phase3_evidence_ladder.tsv")
OUTPUT = Path("results/tables/ambiguous_row_scoring_strategies_species.tsv")
REPORT = Path("results/reports/ambiguous_row_scoring_strategies.md")

VARIANT = "phase2_W3_full_background_sensitivity"
SCORE = "transposon_repeat_suppression_score"
COVERAGE = "transposon_repeat_suppression_coverage"
GENE_COUNT = 31


def main() -> None:
    base = pd.read_csv(BASE, sep="\t")
    base = base.loc[base["score_variant"].eq(VARIANT)].drop_duplicates("scientific_name").copy()
    evidence = pd.read_csv(EVIDENCE, sep="\t", dtype=str)
    ambiguous = evidence.loc[evidence["phase3_scoring_status"].eq("not_scoreable_not_absence")].copy()
    counts = ambiguous.groupby("scientific_name").size().rename("ambiguous_rows_n")
    base = base.merge(counts, on="scientific_name", how="left")
    base["ambiguous_rows_n"] = base["ambiguous_rows_n"].fillna(0).astype(int)
    base["base_weighted_numerator"] = pd.to_numeric(base[SCORE], errors="coerce") * GENE_COUNT
    base["base_found_n"] = np.rint(pd.to_numeric(base[COVERAGE], errors="coerce") * GENE_COUNT).astype(int)

    strategies = {
        "ambiguous_as_missing": {
            "description": "Ambiguous rows contribute zero while the fixed 31-gene denominator is retained.",
            "weight": 0.0,
            "exclude_from_denominator": False,
        },
        "ambiguous_present_like_0.5": {
            "description": "Ambiguous rows contribute a conservative 0.5 confidence weight and count as observed.",
            "weight": 0.5,
            "exclude_from_denominator": False,
        },
        "ambiguous_excluded": {
            "description": "Ambiguous rows contribute zero and are removed from the species-specific denominator.",
            "weight": 0.0,
            "exclude_from_denominator": True,
        },
    }

    rows = []
    for strategy, spec in strategies.items():
        x = base.copy()
        x["ambiguity_strategy"] = strategy
        x["strategy_description"] = spec["description"]
        x["strategy_gene_denominator"] = GENE_COUNT
        if spec["exclude_from_denominator"]:
            x["strategy_gene_denominator"] = GENE_COUNT - x["ambiguous_rows_n"]
        x["strategy_weighted_numerator"] = (
            x["base_weighted_numerator"] + spec["weight"] * x["ambiguous_rows_n"]
        )
        x["strategy_found_n"] = x["base_found_n"]
        if spec["weight"] > 0:
            x["strategy_found_n"] = x["base_found_n"] + x["ambiguous_rows_n"]
        x["strategy_score"] = x["strategy_weighted_numerator"] / x["strategy_gene_denominator"]
        x["strategy_coverage"] = x["strategy_found_n"] / x["strategy_gene_denominator"]
        rows.append(x)

    output = pd.concat(rows, ignore_index=True)
    ordered = [
        "ambiguity_strategy",
        "strategy_description",
        "scientific_name",
        "opentree_tip_label",
        "clade",
        "flight_status",
        "genome_analysis_tier",
        "body_mass_g",
        "log10_body_mass_g",
        "max_lifespan_years",
        "pgls_model_c_mass_clade_residual",
        "ambiguous_rows_n",
        "base_found_n",
        "base_weighted_numerator",
        "strategy_found_n",
        "strategy_weighted_numerator",
        "strategy_gene_denominator",
        "strategy_coverage",
        "strategy_score",
    ]
    output = output[ordered].sort_values(["ambiguity_strategy", "scientific_name"])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT, sep="\t", index=False)

    affected = output.loc[output["ambiguous_rows_n"].gt(0)]
    summary = (
        affected.groupby("ambiguity_strategy")
        .agg(
            affected_species=("scientific_name", "nunique"),
            ambiguous_rows=("ambiguous_rows_n", lambda x: int(x.groupby(affected.loc[x.index, "scientific_name"]).first().sum())),
            mean_score=("strategy_score", "mean"),
            mean_coverage=("strategy_coverage", "mean"),
        )
        .reset_index()
    )
    lines = [
        "# Ambiguous-Row Scoring Strategies",
        "",
        f"The evidence ladder contains {len(ambiguous)} partial/family-ambiguous rows across {ambiguous['scientific_name'].nunique()} bird species. These rows are not treated as confirmed absence. Three prespecified score encodings test the numerical consequence of that decision.",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- {row['ambiguity_strategy']}: affected species = {int(row['affected_species'])}; "
            f"mean affected-species score = {row['mean_score']:.3f}; mean affected-species coverage = {row['mean_coverage']:.3f}."
        )
    lines.extend(
        [
            "",
            "The present-like strategy uses a conservative 0.5 confidence weight rather than treating partial fragments as strict presence. The excluded strategy changes only the denominator. None of these encodings asserts gene loss or functional activity.",
        ]
    )
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT} and {REPORT}")


if __name__ == "__main__":
    main()
