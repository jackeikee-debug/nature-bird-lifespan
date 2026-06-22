"""Build targeted Phase 2 orthology queues for high-priority v2 genes."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def validation_route(row: pd.Series) -> str:
    module = row["maintenance_module_v2"]
    submodule = row["submodule_v2"]
    symbol = row["human_gene_symbol"]
    if module == "transposon_repeat_suppression":
        if submodule == "piRNA_germline_repeat_control":
            return "domain_plus_reciprocal_sequence_piwi_tudor"
        if submodule == "somatic_retroelement_restriction":
            return "reciprocal_sequence_plus_ncbi_symbol"
        return "chromatin_repeat_repressor_crossdb_plus_sequence"
    if module == "chromatin_repression_heterochromatin":
        if symbol.startswith("CBX"):
            return "chromodomain_family_crossdb_then_sequence"
        if symbol in {"SUV39H1", "SUV39H2", "EHMT1", "EHMT2", "EZH1", "EZH2", "SETDB2"}:
            return "set_domain_crossdb_then_sequence"
        return "chromatin_complex_crossdb_then_symbol"
    if row["gene_family_risk"] == "high_paralog_family":
        return "paralog_family_crossdb_then_sequence"
    return "ncbi_symbol_then_crossdb"


def gene_priority(row: pd.Series) -> int:
    module = row["maintenance_module_v2"]
    submodule = row["submodule_v2"]
    symbol = row["human_gene_symbol"]
    if module == "transposon_repeat_suppression" and submodule == "piRNA_germline_repeat_control":
        return 1
    if module == "transposon_repeat_suppression":
        return 2
    if module == "chromatin_repression_heterochromatin":
        return 3
    if row["gene_family_risk"] == "high_paralog_family":
        return 4
    if symbol in {"RAD51C", "RAD51D", "XRCC2", "XRCC3", "FANCA", "FANCC", "FANCD2", "FANCI", "FANCM"}:
        return 5
    return 6


def strictness_tier(row: pd.Series) -> str:
    if row["gene_priority_rank"] <= 3:
        return "strict_required_before_claim"
    if row["gene_family_risk"] == "high_paralog_family":
        return "strict_or_sensitivity_required"
    return "symbol_candidate_allowed_for_feasibility"


def species_priority(row: pd.Series) -> int:
    if row["clade"] == "Aves" and row["genome_analysis_tier"] == "tier1_refseq_annotated_chromosome":
        return 1
    if row["clade"] == "Aves":
        return 2
    if row["flight_status"] == "flighted":
        return 3
    if row["genome_analysis_tier"] == "tier1_refseq_annotated_chromosome":
        return 4
    return 5


def batch_label(gene_rank: int, species_rank: int) -> str:
    if gene_rank <= 2 and species_rank <= 2:
        return "batch1_repeat_genes_birds"
    if gene_rank <= 3 and species_rank <= 2:
        return "batch2_chromatin_genes_birds"
    if gene_rank <= 3 and species_rank <= 4:
        return "batch3_repeat_chromatin_nonbird_controls"
    if gene_rank <= 5 and species_rank <= 2:
        return "batch4_other_paralog_risk_birds"
    return "batch5_sensitivity_remaining"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthology-audit", type=pathlib.Path, required=True)
    parser.add_argument("--species", type=pathlib.Path, required=True)
    parser.add_argument("--gene-queue-output", type=pathlib.Path, required=True)
    parser.add_argument("--query-plan-output", type=pathlib.Path, required=True)
    parser.add_argument("--batch-summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--pilot-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    audit = pd.read_csv(args.orthology_audit, sep="\t")
    species = pd.read_csv(args.species, sep="\t")

    genes = audit[
        audit["orthology_feasibility_class"] == "new_high_priority_validation_required"
    ].copy()
    genes["gene_priority_rank"] = genes.apply(gene_priority, axis=1)
    genes["validation_route"] = genes.apply(validation_route, axis=1)
    genes["claim_strictness_tier"] = genes.apply(strictness_tier, axis=1)
    genes["human_reference_gene_url"] = genes["entrezgene"].apply(
        lambda value: f"https://www.ncbi.nlm.nih.gov/gene/{int(value)}"
        if pd.notna(value) and str(value) not in {"", "nan"}
        else ""
    )
    gene_cols = [
        "gene_priority_rank",
        "human_gene_symbol",
        "validated_symbol",
        "entrezgene",
        "gene_name",
        "maintenance_module_v2",
        "submodule_v2",
        "gene_family_risk",
        "orthology_validation_priority",
        "validation_route",
        "claim_strictness_tier",
        "source_evidence_tags",
        "human_reference_gene_url",
    ]
    genes = genes[gene_cols].sort_values(
        ["gene_priority_rank", "maintenance_module_v2", "submodule_v2", "human_gene_symbol"]
    )

    species_keep = species.copy()
    species_keep["species_priority_rank"] = species_keep.apply(species_priority, axis=1)
    species_cols = [
        "species_priority_rank",
        "scientific_name",
        "clade",
        "flight_status",
        "species_taxid",
        "best_assembly_accession",
        "genome_analysis_tier",
        "genome_quality_risk",
        "has_annotation_report",
        "busco_complete",
        "busco_missing",
        "ftp_path_refseq",
        "ftp_path_genbank",
    ]
    species_keep = species_keep[species_cols]

    query = genes.merge(species_keep, how="cross")
    query["validation_batch"] = [
        batch_label(g, s)
        for g, s in zip(query["gene_priority_rank"], query["species_priority_rank"])
    ]
    query["ncbi_gene_query_term"] = (
        '("' + query["human_gene_symbol"] + '"[Gene Name] OR "'
        + query["human_gene_symbol"]
        + '"[Preferred Symbol]) AND txid'
        + query["species_taxid"].astype(str)
        + "[Organism]"
    )
    query["query_status"] = "not_started"
    query["external_orthology_sources"] = "NCBI_Gene;OMA;OrthoDB;Ensembl_Compara"
    query = query.sort_values(
        [
            "gene_priority_rank",
            "species_priority_rank",
            "maintenance_module_v2",
            "human_gene_symbol",
            "scientific_name",
        ]
    )

    pilot = query[
        query["validation_batch"].isin(
            {"batch1_repeat_genes_birds", "batch2_chromatin_genes_birds"}
        )
    ].copy()
    pilot = (
        pilot.sort_values(
            [
                "gene_priority_rank",
                "species_priority_rank",
                "human_gene_symbol",
                "scientific_name",
            ]
        )
        .groupby("human_gene_symbol", as_index=False)
        .head(5)
    )
    pilot["pilot_reason"] = "top_5_bird_species_per_priority_gene"

    batch_summary = (
        query.groupby("validation_batch", as_index=False)
        .agg(
            rows=("human_gene_symbol", "size"),
            genes=("human_gene_symbol", "nunique"),
            species=("scientific_name", "nunique"),
            aves_rows=("clade", lambda x: int((x == "Aves").sum())),
            tier1_rows=(
                "genome_analysis_tier",
                lambda x: int((x == "tier1_refseq_annotated_chromosome").sum()),
            ),
        )
        .sort_values("validation_batch")
    )

    args.gene_queue_output.parent.mkdir(parents=True, exist_ok=True)
    genes.to_csv(args.gene_queue_output, sep="\t", index=False)
    args.query_plan_output.parent.mkdir(parents=True, exist_ok=True)
    query.to_csv(args.query_plan_output, sep="\t", index=False)
    args.batch_summary_output.parent.mkdir(parents=True, exist_ok=True)
    batch_summary.to_csv(args.batch_summary_output, sep="\t", index=False)
    args.pilot_output.parent.mkdir(parents=True, exist_ok=True)
    pilot.to_csv(args.pilot_output, sep="\t", index=False)

    module_counts = genes.groupby("maintenance_module_v2")["human_gene_symbol"].nunique()
    lines = [
        "# Phase 2 P2.2 Targeted Orthology Queue Report",
        "",
        f"High-priority genes queued: {genes['human_gene_symbol'].nunique()}",
        f"Primary species: {species_keep['scientific_name'].nunique()}",
        f"Full gene-species query rows: {len(query)}",
        f"Pilot query rows: {len(pilot)}",
        "",
        "## High-Priority Genes by Module",
        "",
    ]
    for module, count in module_counts.items():
        lines.append(f"- {module}: {count}")
    lines.extend(["", "## Query Rows by Batch", ""])
    for _, row in batch_summary.iterrows():
        lines.append(
            f"- {row['validation_batch']}: {row['rows']} rows, "
            f"{row['genes']} genes, {row['species']} species"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This queue is a targeted validation plan, not an orthology result. Batch 1 and batch 2 should be executed first because they protect the main repeat/chromatin claim. Downstream scoring should keep strict-required genes separate from feasibility-only symbol candidates.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.query_plan_output} with {len(query)} rows")


if __name__ == "__main__":
    main()
