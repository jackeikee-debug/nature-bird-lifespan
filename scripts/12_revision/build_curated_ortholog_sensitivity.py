#!/usr/bin/env python3
"""Audit strict ortholog rules and build diagnostic module-level scores."""

from pathlib import Path

import pandas as pd


MATRIX = Path("data/processed/ortholog_matrix_primary_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv")
TRAITS = Path("data/processed/maintenance_lifespan_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv")
ROW_OUTPUT = Path("results/tables/curated_ortholog_rule_rows.tsv")
SCORE_OUTPUT = Path("results/tables/curated_ortholog_module_scores.tsv")
SUMMARY_OUTPUT = Path("results/tables/curated_ortholog_rule_clade_summary.tsv")
ELIGIBILITY_OUTPUT = Path("results/tables/curated_ortholog_rule_eligibility.tsv")
REPORT = Path("results/reports/curated_ortholog_sensitivity_audit.md")

VARIANT = "phase2_W3_full_background_sensitivity"
STRICT_PHASE3 = {
    "phase3_gff_sequence_supported",
    "phase3_cds_translation_sequence_supported",
    "phase3_uniprot_sequence_supported",
}


def main() -> None:
    matrix = pd.read_csv(MATRIX, sep="\t", dtype=str)
    traits = pd.read_csv(TRAITS, sep="\t")
    traits = traits.loc[traits["score_variant"].eq(VARIANT)].drop_duplicates("scientific_name")
    matrix["rule_exact_symbol_only"] = matrix["ortholog_status"].eq("exact_symbol")
    matrix["rule_exact_plus_strict_rescue"] = matrix["ortholog_status"].isin({"exact_symbol"} | STRICT_PHASE3)
    matrix["rule_high_confidence_only"] = matrix["week4_candidate_confidence"].isin({"high", "medium_high_external"})
    matrix["rule_reciprocal_sequence_only"] = (
        matrix["week4_candidate_status"].eq("week4_sequence_supported_candidate")
        | matrix["ortholog_status"].isin(STRICT_PHASE3)
    )
    rule_columns = [column for column in matrix.columns if column.startswith("rule_")]

    row_columns = [
        "scientific_name", "clade", "flight_status", "maintenance_module", "human_gene_symbol",
        "ortholog_status", "week4_candidate_status", "week4_candidate_confidence",
    ] + rule_columns
    ROW_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    matrix[row_columns].to_csv(ROW_OUTPUT, sep="\t", index=False)

    score_rows = []
    clade_rows = []
    eligibility_rows = []
    for rule in rule_columns:
        for (species, module), group in matrix.groupby(["scientific_name", "maintenance_module"], sort=False):
            first = group.iloc[0]
            score_rows.append(
                {
                    "ortholog_rule": rule.removeprefix("rule_"),
                    "scientific_name": species,
                    "clade": first["clade"],
                    "flight_status": first["flight_status"],
                    "maintenance_module": module,
                    "genes_total": len(group),
                    "curated_genes_n": int(group[rule].sum()),
                    "curated_fraction": group[rule].mean(),
                }
            )
        for clade, group in matrix.groupby("clade"):
            clade_rows.append(
                {
                    "ortholog_rule": rule.removeprefix("rule_"),
                    "clade": clade,
                    "curated_rows_n": int(group[rule].sum()),
                    "rows_total": len(group),
                    "curated_fraction": group[rule].mean(),
                }
            )
        clade_fractions = matrix.groupby("clade")[rule].mean()
        species_fractions = matrix.groupby("scientific_name")[rule].mean()
        min_clade = clade_fractions.min()
        max_clade = clade_fractions.max()
        overall = matrix[rule].mean()
        adequate = min_clade >= 0.5
        balanced = (max_clade - min_clade) <= 0.2
        nondegenerate = species_fractions.std() >= 0.05
        eligibility_rows.append(
            {
                "ortholog_rule": rule.removeprefix("rule_"),
                "overall_curated_fraction": overall,
                "min_clade_fraction": min_clade,
                "max_clade_fraction": max_clade,
                "max_minus_min_clade_fraction": max_clade - min_clade,
                "species_fraction_sd": species_fractions.std(),
                "adequate_min_clade_coverage": adequate,
                "balanced_across_clades": balanced,
                "nondegenerate_species_variation": nondegenerate,
                "eligible_as_primary_full_matrix_control": adequate and balanced and nondegenerate,
            }
        )

    scores = pd.DataFrame(score_rows).merge(
        traits[
            [
                "scientific_name", "opentree_tip_label", "genome_analysis_tier", "body_mass_g",
                "log10_body_mass_g", "max_lifespan_years", "pgls_model_c_mass_clade_residual",
            ]
        ],
        on="scientific_name",
        how="left",
    )
    summaries = pd.DataFrame(clade_rows)
    eligibility = pd.DataFrame(eligibility_rows)
    scores.to_csv(SCORE_OUTPUT, sep="\t", index=False)
    summaries.to_csv(SUMMARY_OUTPUT, sep="\t", index=False)
    eligibility.to_csv(ELIGIBILITY_OUTPUT, sep="\t", index=False)

    lines = [
        "# Curated Ortholog-Only Full-Matrix Audit",
        "",
        "Four stricter ortholog rules were evaluated before model interpretation. A rule was eligible for primary full-matrix sensitivity only if every clade retained at least 50% of rows, the maximum clade coverage difference was at most 0.20, and species-level fractions retained non-degenerate variation (SD at least 0.05).",
        "",
    ]
    for _, row in eligibility.iterrows():
        lines.append(
            f"- {row['ortholog_rule']}: overall = {row['overall_curated_fraction']:.3f}; "
            f"minimum clade = {row['min_clade_fraction']:.3f}; max-min clade difference = "
            f"{row['max_minus_min_clade_fraction']:.3f}; species SD = {row['species_fraction_sd']:.3f}; "
            f"primary eligible = {bool(row['eligible_as_primary_full_matrix_control'])}."
        )
    eligible_n = int(eligibility["eligible_as_primary_full_matrix_control"].sum())
    lines.extend(
        [
            "",
            f"Primary-eligible rules: {eligible_n}/{len(eligibility)}.",
            "",
            "Exact-symbol and NCBI-centered rules retain more rows but are strongly clade-imbalanced. High-confidence and reciprocal-sequence-only rules are more balanced but cover only a few percent of the 200-gene matrix. Their PGLS fits are reported as diagnostics, not as a replacement primary analysis. The targeted 10-gene sequence and Pfam-domain layers remain the defensible orthology-aware validation route.",
        ]
    )
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {ROW_OUTPUT}, {SCORE_OUTPUT}, {SUMMARY_OUTPUT}, {ELIGIBILITY_OUTPUT}, and {REPORT}")


if __name__ == "__main__":
    main()
