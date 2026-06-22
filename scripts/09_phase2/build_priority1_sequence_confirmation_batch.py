"""Build the priority 1 protein sequence/domain confirmation batch."""

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


DOMAIN_RULES = {
    "ASZ1": "ankyrin_repeat_plus_SAM_basic_region_check",
    "DDX4": "DEAD_box_helicase_domain_check",
    "FKBP6": "FKBP_peptidyl_prolyl_isomerase_domain_check",
    "GTSF1": "CHHC_zinc_finger_domain_check",
    "HENMT1": "RNA_methyltransferase_domain_check",
    "MAEL": "MAEL_domain_or_HMG_like_domain_check",
    "MOV10L1": "RNA_helicase_domain_check_and_MOV10_family_disambiguation",
    "PLD6": "phospholipase_D_nuclease_domain_check",
    "RNF17": "RING_finger_plus_Tudor_related_region_check",
    "TDRKH": "Tudor_domain_related_region_check",
}


def domain_rule(symbol: str) -> str:
    if symbol.startswith("TDRD"):
        return "Tudor_domain_architecture_check_and_TDRD_paralog_disambiguation"
    return DOMAIN_RULES.get(symbol, "reciprocal_sequence_similarity_check")


def confirmation_route(symbol: str) -> str:
    if symbol.startswith("TDRD"):
        return "domain_architecture_plus_reciprocal_sequence"
    return "candidate_protein_sequence_plus_reciprocal_similarity"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ncbi-pilot-results", type=pathlib.Path, required=True)
    parser.add_argument("--upgrade-candidates", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    ncbi = pd.read_csv(args.ncbi_pilot_results, sep="\t")
    upgrades = pd.read_csv(args.upgrade_candidates, sep="\t")
    decisions = upgrades[
        [
            "human_gene_symbol",
            "combined_strict_upgrade_decision",
            "recommended_next_step",
        ]
    ]

    batch = ncbi[
        ncbi["human_gene_symbol"].isin(PRIORITY1_GENES)
        & (ncbi["ncbi_pilot_status"] == "candidate_found")
    ].copy()
    batch = batch.merge(decisions, on="human_gene_symbol", how="left")
    batch["ncbi_gene_id"] = batch["ncbi_gene_id"].astype("Int64").astype(str)
    batch["entrezgene"] = batch["entrezgene"].astype("Int64").astype(str)
    batch["sequence_confirmation_route"] = batch["human_gene_symbol"].map(confirmation_route)
    batch["domain_or_sequence_rule"] = batch["human_gene_symbol"].map(domain_rule)
    batch["protein_fetch_query"] = batch.apply(
        lambda row: f"{row['ncbi_gene_id']}[Gene ID] AND txid{row['species_taxid']}[Organism]",
        axis=1,
    )
    batch["strict_upgrade_allowed_after_this_step"] = False
    batch["absence_claim_allowed"] = False

    keep = [
        "human_gene_symbol",
        "maintenance_module_v2",
        "submodule_v2",
        "scientific_name",
        "species_taxid",
        "best_assembly_accession",
        "ncbi_gene_id",
        "entrezgene",
        "ncbi_gene_symbol",
        "ncbi_gene_description",
        "ncbi_gene_url",
        "ncbi_symbol_confidence",
        "combined_strict_upgrade_decision",
        "recommended_next_step",
        "sequence_confirmation_route",
        "domain_or_sequence_rule",
        "protein_fetch_query",
        "strict_upgrade_allowed_after_this_step",
        "absence_claim_allowed",
    ]
    batch = batch[keep].sort_values(["human_gene_symbol", "scientific_name"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    batch.to_csv(args.output, sep="\t", index=False)

    gene_counts = batch.groupby("human_gene_symbol")["scientific_name"].nunique()
    lines = [
        "# Phase 2 Priority 1 Sequence Confirmation Batch Report",
        "",
        "## Summary",
        "",
        f"Rows: {len(batch)}",
        f"Genes: {batch['human_gene_symbol'].nunique()}",
        f"Species: {batch['scientific_name'].nunique()}",
        "",
        "## Gene Coverage",
        "",
    ]
    for gene, count in gene_counts.items():
        lines.append(f"- {gene}: {count} species")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This is the protein/domain validation batch for priority 1 repeat/piRNA genes. Rows are NCBI Gene-supported candidates, but strict upgrade remains false until protein sequence/domain checks are completed.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
