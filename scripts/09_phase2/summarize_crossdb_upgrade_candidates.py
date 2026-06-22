"""Combine NCBI Gene and UniProt pilot evidence into strict-upgrade decisions."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


PRIORITY1_GENES = {
    "ASZ1",
    "DDX4",
    "FKBP6",
    "GTSF1",
    "HENMT1",
    "MAEL",
    "MOV10L1",
    "PLD6",
    "RNF17",
    "TDRD1",
    "TDRD12",
    "TDRD5",
    "TDRD6",
    "TDRD7",
    "TDRD9",
    "TDRKH",
}


def decision(row: pd.Series) -> str:
    ncbi_fraction = float(row.get("ncbi_candidate_fraction", 0))
    uniprot_fraction = float(row.get("uniprot_symbol_like_fraction", 0))
    any_domain_route = "domain_check" in str(row.get("crossdb_route_examples", ""))
    if ncbi_fraction >= 0.8 and uniprot_fraction >= 0.8 and not any_domain_route:
        return "strict_upgrade_candidate_after_matrix_merge"
    if ncbi_fraction >= 0.8 and any_domain_route:
        return "sequence_or_domain_confirmation_required_high_paralog"
    if ncbi_fraction >= 0.8 and uniprot_fraction < 0.8:
        return "sequence_confirmation_required_uniprot_sparse"
    if ncbi_fraction >= 0.6:
        return "expand_crossdb_species_then_sequence_check"
    return "hold_manual_review_not_absence"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ncbi-summary", type=pathlib.Path, required=True)
    parser.add_argument("--uniprot-summary", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    ncbi = pd.read_csv(args.ncbi_summary, sep="\t")
    uni = pd.read_csv(args.uniprot_summary, sep="\t")

    ncbi = ncbi[ncbi["human_gene_symbol"].isin(PRIORITY1_GENES)].copy()
    ncbi = ncbi.rename(
        columns={
            "candidate_fraction": "ncbi_candidate_fraction",
            "candidate_found": "ncbi_candidate_found_species",
            "no_candidate": "ncbi_no_candidate_species",
            "candidate_exact_symbol": "ncbi_exact_symbol_species",
        }
    )
    keep_ncbi = [
        "human_gene_symbol",
        "ncbi_candidate_fraction",
        "ncbi_candidate_found_species",
        "ncbi_no_candidate_species",
        "ncbi_exact_symbol_species",
    ]
    for col in keep_ncbi:
        if col not in ncbi.columns:
            ncbi[col] = pd.NA

    merged = uni.merge(ncbi[keep_ncbi], on="human_gene_symbol", how="left")
    merged["ncbi_candidate_fraction"] = merged["ncbi_candidate_fraction"].fillna(0)
    merged["combined_strict_upgrade_decision"] = merged.apply(decision, axis=1)
    merged["recommended_next_step"] = merged["combined_strict_upgrade_decision"].map(
        {
            "strict_upgrade_candidate_after_matrix_merge": "merge_into_confirmed_crossdb_matrix_then_recompute_strict_panel",
            "sequence_or_domain_confirmation_required_high_paralog": "run_domain_architecture_or_reciprocal_sequence_check_before_strict_upgrade",
            "sequence_confirmation_required_uniprot_sparse": "fetch_candidate_proteins_and_run_sequence_confirmation_before_strict_upgrade",
            "expand_crossdb_species_then_sequence_check": "add_more_species_or_databases_before_sequence_confirmation",
            "hold_manual_review_not_absence": "hold_out_of_strict_and_do_not_score_absence",
        }
    )

    sort_cols = [
        "crossdb_priority",
        "combined_strict_upgrade_decision",
        "human_gene_symbol",
    ]
    merged = merged.sort_values(sort_cols)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)

    counts = merged["combined_strict_upgrade_decision"].value_counts().sort_index()
    lines = [
        "# Phase 2 Cross-Database Upgrade Candidate Report",
        "",
        "## Summary",
        "",
        f"Priority 1 genes assessed: {merged['human_gene_symbol'].nunique()}",
        "",
        "## Combined Decisions",
        "",
    ]
    for key, count in counts.items():
        lines.append(f"- {key}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "NCBI Gene support is strong for the priority 1 piRNA/repeat genes, but UniProt is sparse across the five pilot birds. Therefore these genes should not be upgraded directly into strict v2 from UniProt alone. The next evidence layer should be targeted protein sequence and domain confirmation, with extra caution for TDRD-family paralogs.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
