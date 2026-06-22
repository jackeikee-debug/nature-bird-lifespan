"""Audit P2.4 scoring readiness and build full-species expansion queues."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def presence_status(row: pd.Series) -> str:
    for col in ["week4_candidate_status", "final_candidate_status", "combined_candidate_status"]:
        value = str(row.get(col, ""))
        if value:
            return value
    return "missing_matrix_row"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eligibility", type=pathlib.Path, required=True)
    parser.add_argument("--ortholog-matrix", type=pathlib.Path, required=True)
    parser.add_argument("--species", type=pathlib.Path, required=True)
    parser.add_argument("--gene-audit-output", type=pathlib.Path, required=True)
    parser.add_argument("--expansion-queue-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    eligibility = pd.read_csv(args.eligibility, sep="\t")
    matrix = pd.read_csv(args.ortholog_matrix, sep="\t")
    species = pd.read_csv(args.species, sep="\t")
    expected_species = species["scientific_name"].nunique()

    matrix["presence_status_for_audit"] = matrix.apply(presence_status, axis=1)
    found_statuses = {
        "ncbi_gene_candidate",
        "gff_rescue_candidate",
        "diamond_validated_protein_candidate",
        "week4_sequence_supported_candidate",
    }
    matrix["candidate_found_for_audit"] = matrix["presence_status_for_audit"].isin(found_statuses)

    gene_rows = []
    for _, gene in eligibility.iterrows():
        symbol = gene["human_gene_symbol"]
        sub = matrix[matrix["human_gene_symbol"] == symbol]
        species_rows = sub["scientific_name"].nunique()
        found_species = sub.loc[sub["candidate_found_for_audit"], "scientific_name"].nunique()
        if species_rows == expected_species:
            readiness = "species_matrix_complete"
        elif species_rows > 0:
            readiness = "species_matrix_partial"
        else:
            readiness = "species_matrix_missing"
        if bool(gene["strict_score_allowed"]) and readiness != "species_matrix_complete":
            next_step = "expand_to_full_species_matrix_before_pgls"
            pgls_ready = False
        elif bool(gene["strict_score_allowed"]):
            next_step = "ready_for_strict_species_level_scoring"
            pgls_ready = True
        elif bool(gene["sensitivity_score_allowed"]) and readiness == "species_matrix_complete":
            next_step = "ready_for_sensitivity_species_level_scoring"
            pgls_ready = True
        elif bool(gene["sensitivity_score_allowed"]):
            next_step = "expand_if_included_in_sensitivity_pgls"
            pgls_ready = False
        else:
            next_step = "not_ready_not_in_scoring_panel"
            pgls_ready = False
        gene_rows.append(
            {
                "human_gene_symbol": symbol,
                "maintenance_module_v2": gene["maintenance_module_v2"],
                "submodule_v2": gene["submodule_v2"],
                "v2_scoring_group": gene["v2_scoring_group"],
                "strict_score_allowed": gene["strict_score_allowed"],
                "sensitivity_score_allowed": gene["sensitivity_score_allowed"],
                "absence_scoring_allowed": gene["absence_scoring_allowed"],
                "claim_use": gene["claim_use"],
                "expected_species": expected_species,
                "species_rows_in_current_matrix": species_rows,
                "candidate_found_species_in_current_matrix": found_species,
                "matrix_species_coverage_fraction": species_rows / expected_species,
                "candidate_found_species_fraction": found_species / expected_species,
                "p2_4_readiness": readiness,
                "p2_4_pgls_ready_now": pgls_ready,
                "p2_4_next_step": next_step,
            }
        )
    gene_audit = pd.DataFrame(gene_rows)

    priority_groups = {
        "strict_sequence_supported",
        "domain_supported_paralog_guard",
        "domain_supported_manual_upgrade_candidate",
    }
    priority_genes = eligibility[eligibility["v2_scoring_group"].isin(priority_groups)].copy()
    queue_rows = []
    existing_pairs = set(zip(matrix["human_gene_symbol"], matrix["scientific_name"]))
    for _, gene in priority_genes.iterrows():
        for _, sp in species.iterrows():
            pair = (gene["human_gene_symbol"], sp["scientific_name"])
            queue_rows.append(
                {
                    "human_gene_symbol": gene["human_gene_symbol"],
                    "maintenance_module_v2": gene["maintenance_module_v2"],
                    "submodule_v2": gene["submodule_v2"],
                    "v2_scoring_group": gene["v2_scoring_group"],
                    "claim_use": gene["claim_use"],
                    "strict_score_allowed": gene["strict_score_allowed"],
                    "sensitivity_score_allowed": gene["sensitivity_score_allowed"],
                    "scientific_name": sp["scientific_name"],
                    "clade": sp["clade"],
                    "flight_status": sp["flight_status"],
                    "species_taxid": sp["species_taxid"],
                    "best_assembly_accession": sp["best_assembly_accession"],
                    "genome_analysis_tier": sp["genome_analysis_tier"],
                    "genome_quality_risk": sp["genome_quality_risk"],
                    "has_annotation_report": sp["has_annotation_report"],
                    "busco_complete": sp["busco_complete"],
                    "already_in_current_matrix": pair in existing_pairs,
                    "ncbi_gene_query_term": f'("{gene["human_gene_symbol"]}"[Gene Name] OR "{gene["human_gene_symbol"]}"[Preferred Symbol]) AND txid{sp["species_taxid"]}[Organism]',
                    "recommended_route": "NCBI_Gene_then_candidate_protein_sequence_confirmation",
                }
            )
    queue = pd.DataFrame(queue_rows)

    args.gene_audit_output.parent.mkdir(parents=True, exist_ok=True)
    args.expansion_queue_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    gene_audit.to_csv(args.gene_audit_output, sep="\t", index=False)
    queue.to_csv(args.expansion_queue_output, sep="\t", index=False)

    readiness_counts = gene_audit["p2_4_readiness"].value_counts().sort_index()
    strict = gene_audit[gene_audit["strict_score_allowed"].astype(bool)]
    sensitivity = gene_audit[gene_audit["sensitivity_score_allowed"].astype(bool)]
    queue_new = queue[~queue["already_in_current_matrix"].astype(bool)]
    lines = [
        "# Phase 2 P2.4 Scoring Readiness Audit",
        "",
        "## Summary",
        "",
        f"Expected species in primary mechanism panel: {expected_species}",
        f"Eligibility genes assessed: {gene_audit['human_gene_symbol'].nunique()}",
        f"Strict genes: {strict['human_gene_symbol'].nunique()}",
        f"Strict genes PGLS-ready now: {int(strict['p2_4_pgls_ready_now'].sum())}",
        f"Sensitivity genes: {sensitivity['human_gene_symbol'].nunique()}",
        f"Sensitivity genes PGLS-ready now: {int(sensitivity['p2_4_pgls_ready_now'].sum())}",
        f"Priority-1 expansion queue rows: {len(queue)}",
        f"Priority-1 expansion rows not already in matrix: {len(queue_new)}",
        "",
        "## Readiness Classes",
        "",
    ]
    for readiness, count in readiness_counts.items():
        lines.append(f"- {readiness}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The sequence-updated strict panel is evidence-ready at the gene level, but the 7 newly upgraded strict genes still require full-species matrix expansion before formal PGLS. Until that expansion is complete, PGLS should use the previous 41-gene species-complete matrix or be labelled as pilot-only.",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.gene_audit_output}")


if __name__ == "__main__":
    main()
