"""Build cross-database confirmation queue for v2 crossdb_confirm genes."""

from __future__ import annotations

import argparse
import pathlib
import urllib.parse

import pandas as pd


def quote(value: str) -> str:
    return urllib.parse.quote(str(value))


def route(row: pd.Series) -> str:
    module = row["maintenance_module_v2"]
    gene = row["human_gene_symbol"]
    if module == "chromatin_repression_heterochromatin":
        return "OMA_OrthoDB_Ensembl_symbol_then_domain_family_check"
    if gene.startswith("TDRD"):
        return "OMA_OrthoDB_symbol_plus_Tudor_domain_check"
    if gene in {"DNMT1", "DNMT3A", "DNMT3B", "SETDB2", "HELLS", "UHRF1"}:
        return "OMA_OrthoDB_symbol_plus_chromatin_domain_check"
    return "OMA_OrthoDB_symbol_then_sequence_confirmation"


def priority(row: pd.Series) -> int:
    module = row["maintenance_module_v2"]
    submodule = row["submodule_v2"]
    if module == "transposon_repeat_suppression" and submodule == "piRNA_germline_repeat_control":
        return 1
    if module == "transposon_repeat_suppression":
        return 2
    if module == "chromatin_repression_heterochromatin":
        return 3
    return 4


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eligibility", type=pathlib.Path, required=True)
    parser.add_argument("--species", type=pathlib.Path, required=True)
    parser.add_argument("--gene-output", type=pathlib.Path, required=True)
    parser.add_argument("--species-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    eligibility = pd.read_csv(args.eligibility, sep="\t")
    species = pd.read_csv(args.species, sep="\t")
    genes = eligibility[eligibility["v2_scoring_group"] == "crossdb_confirm"].copy()
    genes["crossdb_priority"] = genes.apply(priority, axis=1)
    genes["crossdb_route"] = genes.apply(route, axis=1)
    genes["oma_symbol_url"] = genes["human_gene_symbol"].apply(
        lambda g: f"https://omabrowser.org/oma/search/?type=Gene&q={quote(g)}"
    )
    genes["orthodb_symbol_url"] = genes["human_gene_symbol"].apply(
        lambda g: f"https://www.orthodb.org/?query={quote(g)}"
    )
    genes["uniprot_symbol_url"] = genes["human_gene_symbol"].apply(
        lambda g: f"https://www.uniprot.org/uniprotkb?query={quote(g + ' Aves OR bird')}"
    )
    genes["ncbi_gene_human_url"] = genes["human_gene_symbol"].apply(
        lambda g: f"https://www.ncbi.nlm.nih.gov/gene/?term={quote(g + '[sym] AND Homo sapiens[orgn]')}"
    )

    gene_cols = [
        "crossdb_priority",
        "human_gene_symbol",
        "maintenance_module_v2",
        "submodule_v2",
        "v2_scoring_group",
        "candidate_fraction",
        "triage_class",
        "gene_family_risk",
        "crossdb_route",
        "oma_symbol_url",
        "orthodb_symbol_url",
        "uniprot_symbol_url",
        "ncbi_gene_human_url",
        "next_validation_step",
        "claim_use",
    ]
    gene_queue = genes[gene_cols].sort_values(
        ["crossdb_priority", "maintenance_module_v2", "submodule_v2", "human_gene_symbol"]
    )
    args.gene_output.parent.mkdir(parents=True, exist_ok=True)
    gene_queue.to_csv(args.gene_output, sep="\t", index=False)

    pilot_species = species[
        (species["clade"] == "Aves")
        & (species["genome_analysis_tier"] == "tier1_refseq_annotated_chromosome")
    ].copy()
    pilot_species = pilot_species.sort_values(
        ["genome_quality_risk", "busco_missing", "scientific_name"]
    ).head(5)
    rows = []
    for _, gene in gene_queue.iterrows():
        for _, sp in pilot_species.iterrows():
            rows.append(
                {
                    "crossdb_priority": gene["crossdb_priority"],
                    "human_gene_symbol": gene["human_gene_symbol"],
                    "maintenance_module_v2": gene["maintenance_module_v2"],
                    "submodule_v2": gene["submodule_v2"],
                    "scientific_name": sp["scientific_name"],
                    "species_taxid": sp["species_taxid"],
                    "best_assembly_accession": sp["best_assembly_accession"],
                    "crossdb_route": gene["crossdb_route"],
                    "oma_species_symbol_url": f"https://omabrowser.org/oma/search/?type=Gene&q={quote(gene['human_gene_symbol'] + ' ' + sp['scientific_name'])}",
                    "orthodb_species_symbol_url": f"https://www.orthodb.org/?query={quote(gene['human_gene_symbol'] + ' ' + sp['scientific_name'])}",
                    "uniprot_species_symbol_url": f"https://www.uniprot.org/uniprotkb?query={quote(gene['human_gene_symbol'] + ' ' + sp['scientific_name'])}",
                    "manual_crossdb_decision": "",
                    "manual_crossdb_notes": "",
                }
            )
    species_queue = pd.DataFrame(rows).sort_values(
        ["crossdb_priority", "human_gene_symbol", "scientific_name"]
    )
    args.species_output.parent.mkdir(parents=True, exist_ok=True)
    species_queue.to_csv(args.species_output, sep="\t", index=False)

    module_counts = gene_queue.groupby("maintenance_module_v2")["human_gene_symbol"].nunique()
    lines = [
        "# Phase 2 Cross-Database Confirmation Queue Report",
        "",
        f"crossdb_confirm genes: {len(gene_queue)}",
        f"pilot species rows: {len(species_queue)}",
        "",
        "## Genes by Module",
        "",
    ]
    for module, count in module_counts.items():
        lines.append(f"- {module}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "These genes have enough pilot support to be promising, but they are not strict-ready. Cross-database confirmation should upgrade only well-supported rows into the strict v2 score.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.gene_output}")


if __name__ == "__main__":
    main()
