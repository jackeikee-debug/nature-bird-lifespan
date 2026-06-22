"""Build Phase 2 orthology feasibility audit for expanded panel v2."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


HIGH_PARALOG_PATTERNS = (
    "PIWIL",
    "TDRD",
    "APOBEC",
    "MORC",
    "TRIM",
    "CBX",
    "HSPA",
    "HSP",
    "PSM",
    "FANC",
    "RAD51",
    "RNASEH",
    "CASP",
)


def gene_family_risk(symbol: str, module: str) -> str:
    if any(symbol.startswith(pattern) for pattern in HIGH_PARALOG_PATTERNS):
        return "high_paralog_family"
    if module in {
        "transposon_repeat_suppression",
        "chromatin_repression_heterochromatin",
    }:
        return "module_high_priority"
    return "standard"


def summarize_seed_coverage(matrix: pd.DataFrame) -> pd.DataFrame:
    tested = matrix.copy()
    supported_status = {
        "sequence_supported",
        "sequence_supported_weak",
        "not_tested_non_transposon",
    }
    tested["has_candidate"] = tested["week4_candidate_status"].fillna("").ne(
        "candidate_not_found"
    )
    tested["sequence_or_candidate_supported"] = (
        tested["week4_sequence_status"].fillna("").isin(supported_status)
        | tested["has_candidate"]
    )
    return (
        tested.groupby("human_gene_symbol", as_index=False)
        .agg(
            species_rows=("scientific_name", "nunique"),
            candidate_supported_species=("sequence_or_candidate_supported", "sum"),
            tier1_species=(
                "genome_analysis_tier",
                lambda x: int((x == "tier1_refseq_annotated_chromosome").sum()),
            ),
            sequence_supported_species=(
                "week4_sequence_status",
                lambda x: int((x == "sequence_supported").sum()),
            ),
            sequence_weak_species=(
                "week4_sequence_status",
                lambda x: int((x == "sequence_supported_weak").sum()),
            ),
            sequence_not_supported_species=(
                "week4_sequence_status",
                lambda x: int((x == "sequence_not_supported").sum()),
            ),
        )
        .assign(
            candidate_supported_fraction=lambda d: d["candidate_supported_species"]
            / d["species_rows"]
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validated-panel", type=pathlib.Path, required=True)
    parser.add_argument("--matrix", type=pathlib.Path, required=True)
    parser.add_argument("--species", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--module-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    panel = pd.read_csv(args.validated_panel, sep="\t")
    matrix = pd.read_csv(args.matrix, sep="\t")
    species = pd.read_csv(args.species, sep="\t")
    n_species = species["scientific_name"].nunique()
    seed_cov = summarize_seed_coverage(matrix)

    out = panel.merge(seed_cov, on="human_gene_symbol", how="left")
    for col in [
        "species_rows",
        "candidate_supported_species",
        "tier1_species",
        "sequence_supported_species",
        "sequence_weak_species",
        "sequence_not_supported_species",
    ]:
        out[col] = out[col].fillna(0).astype(int)
    out["candidate_supported_fraction"] = out["candidate_supported_fraction"].fillna(0.0)
    out["expected_species_panel_size"] = n_species
    out["gene_family_risk"] = [
        gene_family_risk(sym, mod)
        for sym, mod in zip(out["human_gene_symbol"], out["maintenance_module_v2"])
    ]
    out["orthology_data_status"] = out["seed_status"].map(
        {
            "seed_v0": "existing_seed_matrix_observed",
            "expanded_v2_candidate": "pending_v2_mapping",
        }
    )
    out["orthology_feasibility_class"] = "standard_mapping_candidate"
    out.loc[
        (out["seed_status"] == "seed_v0")
        & (out["candidate_supported_fraction"] >= 0.80),
        "orthology_feasibility_class",
    ] = "observed_high_coverage_seed"
    out.loc[
        (out["seed_status"] == "seed_v0")
        & (out["candidate_supported_fraction"] < 0.80),
        "orthology_feasibility_class",
    ] = "observed_seed_needs_review"
    out.loc[
        (out["seed_status"] != "seed_v0")
        & (
            (out["orthology_validation_priority"] == "high")
            | (out["gene_family_risk"] == "high_paralog_family")
        ),
        "orthology_feasibility_class",
    ] = "new_high_priority_validation_required"

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    module = (
        out.groupby("maintenance_module_v2", as_index=False)
        .agg(
            n_genes=("human_gene_symbol", "nunique"),
            n_seed=("seed_status", lambda x: int((x == "seed_v0").sum())),
            n_pending=("orthology_data_status", lambda x: int((x == "pending_v2_mapping").sum())),
            n_high_priority=(
                "orthology_feasibility_class",
                lambda x: int((x == "new_high_priority_validation_required").sum()),
            ),
            mean_observed_seed_coverage=(
                "candidate_supported_fraction",
                lambda x: float(x[x > 0].mean()) if (x > 0).any() else 0.0,
            ),
        )
        .sort_values("maintenance_module_v2")
    )
    args.module_output.parent.mkdir(parents=True, exist_ok=True)
    module.to_csv(args.module_output, sep="\t", index=False)

    seed = out[out["seed_status"] == "seed_v0"]
    high_priority = int(
        (out["orthology_feasibility_class"] == "new_high_priority_validation_required").sum()
    )
    seed_high = int((seed["candidate_supported_fraction"] >= 0.80).sum())
    gate = "caution"
    if seed_high / max(len(seed), 1) >= 0.90 and high_priority <= 90:
        gate = "pass_for_mapping_design"
    if high_priority > 120:
        gate = "fail_panel_too_ambiguous"

    lines = [
        "# Phase 2 P2.2 Orthology Feasibility Audit Report",
        "",
        f"Primary species panel size: {n_species}",
        f"Expanded panel genes: {len(out)}",
        f"Seed genes with existing orthology observations: {len(seed)}",
        f"Seed genes with >=80% candidate-supported species: {seed_high}/{len(seed)}",
        f"New high-priority validation-required genes: {high_priority}",
        "",
        f"Decision gate: **{gate}**",
        "",
        "Interpretation: this gate does not claim that all 236 genes are orthology-validated. It determines whether the panel is feasible enough to proceed to targeted v2 mapping and which genes need strict/sensitivity handling.",
    ]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} with gate={gate}")


if __name__ == "__main__":
    main()
