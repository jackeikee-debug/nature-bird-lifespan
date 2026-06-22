"""Build external validation plan for Phase 2 P2.2 pilot review genes."""

from __future__ import annotations

import argparse
import pathlib
import urllib.parse

import pandas as pd


DOMAIN_RULES = {
    "APOBEC3A": "APOBEC-like cytidine deaminase domain; lineage-specific APOBEC3 expansion/loss check",
    "APOBEC3B": "APOBEC-like cytidine deaminase domain; lineage-specific APOBEC3 expansion/loss check",
    "APOBEC3C": "APOBEC-like cytidine deaminase domain; lineage-specific APOBEC3 expansion/loss check",
    "APOBEC3F": "APOBEC-like cytidine deaminase domain; lineage-specific APOBEC3 expansion/loss check",
    "APOBEC3G": "APOBEC-like cytidine deaminase domain; lineage-specific APOBEC3 expansion/loss check",
    "APOBEC3H": "APOBEC-like cytidine deaminase domain; lineage-specific APOBEC3 expansion/loss check",
    "PIWIL3": "PAZ+PIWI domain architecture; distinguish PIWIL1/2/4 paralogs",
    "PIWIL4": "PAZ+PIWI domain architecture; distinguish PIWIL1/2/3 paralogs",
    "DNMT3L": "DNMT3-like regulatory protein; distinguish true absence from DNMT3A/B annotation collapse",
    "TREX1": "3-prime repair exonuclease domain; check symbol synonym and gene model fragmentation",
    "ZCCHC3": "zinc finger CCHC domain; check gene model fragmentation and naming",
    "MBD1": "methyl-CpG binding domain; distinguish MBD family paralogs",
    "MORC4": "MORC ATPase/CW-type zinc finger architecture; distinguish MORC family paralogs",
    "RING1": "polycomb RING finger; distinguish RING1/RNF2 paralogs",
    "EHMT2": "SET domain histone methyltransferase; distinguish EHMT1/EHMT2",
    "SUV39H1": "SET domain histone methyltransferase; distinguish SUV39H1/SUV39H2",
}


def quote(value: str) -> str:
    return urllib.parse.quote(str(value))


def external_priority(row: pd.Series) -> int:
    gene = row["human_gene_symbol"]
    if gene.startswith("APOBEC3"):
        return 1
    if gene in {"PIWIL3", "PIWIL4"}:
        return 1
    if gene in {"DNMT3L", "TREX1", "ZCCHC3"}:
        return 2
    if row["candidate_fraction"] == 0:
        return 2
    return 3


def decision_rule(row: pd.Series) -> str:
    if row["candidate_fraction"] == 0:
        return "exclude_from_strict_absence_scoring_until_external_or_domain_support"
    return "hold_for_crossdb_resolution_before_strict_scoring"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--external-review", type=pathlib.Path, required=True)
    parser.add_argument("--pilot-results", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--species-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    review = pd.read_csv(args.external_review, sep="\t")
    pilot = pd.read_csv(args.pilot_results, sep="\t")

    review = review.copy()
    review["external_review_priority"] = review.apply(external_priority, axis=1)
    review["domain_validation_rule"] = review["human_gene_symbol"].map(DOMAIN_RULES).fillna(
        "family/domain architecture check"
    )
    review["strict_scoring_rule"] = review.apply(decision_rule, axis=1)
    review["oma_symbol_search_url"] = review["human_gene_symbol"].apply(
        lambda g: f"https://omabrowser.org/oma/search/?type=Gene&q={quote(g)}"
    )
    review["orthodb_symbol_search_url"] = review["human_gene_symbol"].apply(
        lambda g: f"https://www.orthodb.org/?query={quote(g)}"
    )
    review["uniprot_symbol_search_url"] = review["human_gene_symbol"].apply(
        lambda g: f"https://www.uniprot.org/uniprotkb?query={quote(g + ' Aves OR bird')}"
    )
    review["ncbi_gene_human_url"] = review["human_gene_symbol"].apply(
        lambda g: f"https://www.ncbi.nlm.nih.gov/gene/?term={quote(g + '[sym] AND Homo sapiens[orgn]')}"
    )

    species_rows = []
    for _, row in review.iterrows():
        gene = row["human_gene_symbol"]
        gene_pilot = pilot[pilot["human_gene_symbol"] == gene]
        for _, prow in gene_pilot.iterrows():
            species = prow["scientific_name"]
            species_query = species.replace(" ", "+")
            species_rows.append(
                {
                    "external_review_priority": row["external_review_priority"],
                    "human_gene_symbol": gene,
                    "maintenance_module_v2": row["maintenance_module_v2"],
                    "submodule_v2": row["submodule_v2"],
                    "scientific_name": species,
                    "species_taxid": prow["species_taxid"],
                    "best_assembly_accession": prow["best_assembly_accession"],
                    "ncbi_pilot_status": prow["ncbi_pilot_status"],
                    "ncbi_gene_id": prow.get("ncbi_gene_id", ""),
                    "ncbi_gene_symbol": prow.get("ncbi_gene_symbol", ""),
                    "domain_validation_rule": row["domain_validation_rule"],
                    "oma_species_symbol_url": f"https://omabrowser.org/oma/search/?type=Gene&q={quote(gene + ' ' + species)}",
                    "orthodb_species_symbol_url": f"https://www.orthodb.org/?query={quote(gene + ' ' + species)}",
                    "uniprot_species_symbol_url": f"https://www.uniprot.org/uniprotkb?query={quote(gene + ' ' + species)}",
                    "ncbi_gene_species_url": f"https://www.ncbi.nlm.nih.gov/gene/?term={quote(gene + '[sym] AND ' + species + '[orgn]')}",
                    "ncbi_protein_species_url": f"https://www.ncbi.nlm.nih.gov/protein/?term={quote(gene + ' ' + species_query)}",
                    "manual_external_decision": "",
                    "manual_external_notes": "",
                }
            )
    species_plan = pd.DataFrame(species_rows).sort_values(
        ["external_review_priority", "human_gene_symbol", "scientific_name"]
    )

    review = review.sort_values(
        ["external_review_priority", "candidate_fraction", "human_gene_symbol"]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    review.to_csv(args.output, sep="\t", index=False)
    args.species_output.parent.mkdir(parents=True, exist_ok=True)
    species_plan.to_csv(args.species_output, sep="\t", index=False)

    priority_counts = review["external_review_priority"].value_counts().sort_index()
    lines = [
        "# Phase 2 P2.2 External Review Plan Report",
        "",
        f"Review genes: {len(review)}",
        f"Gene-species review rows: {len(species_plan)}",
        "",
        "## Priority Counts",
        "",
    ]
    for priority, count in priority_counts.items():
        lines.append(f"- priority {priority}: {count}")
    lines.extend(
        [
            "",
            "## Key Rule",
            "",
            "Genes with zero bird NCBI symbol hits are not counted as absences. They remain excluded from strict absence scoring until external orthology or domain-level evidence is available.",
            "",
            "## Immediate Focus",
            "",
            "Priority 1 genes are APOBEC3 family members and PIWIL3/PIWIL4. These are the highest collision-risk genes for the repeat-suppression story because they may reflect lineage-specific biology, naming divergence, or annotation gaps.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} and {args.species_output}")


if __name__ == "__main__":
    main()
