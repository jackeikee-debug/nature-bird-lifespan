"""Create the observed 68 x 200 weighted matrix for degradation experiments."""

from __future__ import annotations

import argparse
import pathlib

import numpy as np
import pandas as pd


FOUND_STATUSES = {
    "ncbi_gene_candidate",
    "gff_rescue_candidate",
    "diamond_validated_protein_candidate",
    "week4_sequence_supported_candidate",
}
WEIGHTS = {"high": 1.0, "medium": 0.8, "low": 0.5, "": 0.0}
BASE_GROUPS = {
    "strict_ready",
    "strict_sequence_supported",
    "domain_supported_paralog_guard",
    "domain_supported_manual_upgrade_candidate",
    "crossdb_confirm",
}
W3_MODULES = {
    "DNA_repair_replication_stress",
    "proteostasis_autophagy_mitophagy",
    "cancer_surveillance_senescence",
    "inflammation_innate_immune_restraint",
}


def first_text(row: pd.Series, columns: list[str]) -> str:
    for column in columns:
        value = row.get(column, "")
        if isinstance(value, str) and value and value != "nan":
            return value
    return ""


def gene_weight(row: pd.Series) -> float:
    status = first_text(row, ["week4_candidate_status", "final_candidate_status", "combined_candidate_status"])
    if status not in FOUND_STATUSES:
        return 0.0
    confidence = first_text(row, ["week4_candidate_confidence", "final_candidate_confidence", "ortholog_confidence"])
    return WEIGHTS.get(confidence, 0.0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, default=pathlib.Path("data/processed/ortholog_matrix_primary_phase2_W3_full_background_expanded.tsv"))
    parser.add_argument("--eligibility", type=pathlib.Path, default=pathlib.Path("data/processed/phase2_strict_v2_scoring_eligibility_sequence_updated.tsv"))
    parser.add_argument("--analysis", type=pathlib.Path, default=pathlib.Path("data/processed/jme_revision2_analysis_table.tsv"))
    parser.add_argument("--output", type=pathlib.Path, default=pathlib.Path("data/processed/jme_revision2_degradation_matrix.tsv"))
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t", low_memory=False)
    eligibility = pd.read_csv(args.eligibility, sep="\t")
    analysis = pd.read_csv(args.analysis, sep="\t")
    panel_mask = eligibility["v2_scoring_group"].isin(BASE_GROUPS) | (
        eligibility["maintenance_module_v2"].isin(W3_MODULES)
        & eligibility["v2_scoring_group"].eq("standard_mapping_pending")
    )
    panel = eligibility.loc[panel_mask, ["human_gene_symbol", "maintenance_module_v2"]].drop_duplicates()
    if len(panel) != 200:
        raise RuntimeError(f"Expected 200 panel genes, observed {len(panel)}")
    species = analysis["scientific_name"].tolist()
    selected = matrix[matrix["human_gene_symbol"].isin(panel["human_gene_symbol"])].copy()
    selected = selected[selected["scientific_name"].isin(species)]
    selected["observed_weight"] = selected.apply(gene_weight, axis=1)
    selected = selected.drop(columns=["maintenance_module_v2"], errors="ignore")
    selected = selected.merge(panel, on="human_gene_symbol", how="inner")
    selected = selected[["scientific_name", "human_gene_symbol", "maintenance_module_v2", "observed_weight"]]
    selected = selected.drop_duplicates(["scientific_name", "human_gene_symbol"])
    complete_index = pd.MultiIndex.from_product(
        [species, panel["human_gene_symbol"].tolist()], names=["scientific_name", "human_gene_symbol"]
    )
    selected = selected.set_index(["scientific_name", "human_gene_symbol"]).reindex(complete_index).reset_index()
    selected = selected.merge(panel, on="human_gene_symbol", how="left", suffixes=("", "_fill"))
    selected["maintenance_module_v2"] = selected["maintenance_module_v2"].fillna(selected["maintenance_module_v2_fill"])
    selected["observed_weight"] = selected["observed_weight"].fillna(0.0)
    selected = selected.drop(columns=["maintenance_module_v2_fill"])
    if len(selected) != 68 * 200:
        raise RuntimeError(f"Expected 13,600 rows, observed {len(selected)}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {args.output}: {len(selected)} rows, {(selected.observed_weight > 0).sum()} observed gene rows")


if __name__ == "__main__":
    main()
