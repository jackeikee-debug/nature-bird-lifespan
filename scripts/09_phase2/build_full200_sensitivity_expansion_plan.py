"""Build a staged expansion plan for the full 200-gene sensitivity panel."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DONE_GROUPS = {
    "strict_ready",
    "strict_sequence_supported",
    "domain_supported_paralog_guard",
    "domain_supported_manual_upgrade_candidate",
}


def route_for_group(group: str) -> tuple[str, str, str, str]:
    if group == "strict_ready":
        return (
            "W0_species_matrix_complete",
            "complete",
            "already_species_complete_in_week4_matrix",
            "eligible_for_main_strict_or_baseline_sensitivity",
        )
    if group == "strict_sequence_supported":
        return (
            "W1_priority1_sequence_expanded",
            "complete",
            "priority1_NCBI_Gene_expansion_done_after_DIAMOND_BLAST_support",
            "eligible_for_strict48_presence_score_no_absence_claim",
        )
    if group in {"domain_supported_paralog_guard", "domain_supported_manual_upgrade_candidate"}:
        return (
            "W1_priority1_domain_expanded",
            "complete",
            "priority1_NCBI_Gene_expansion_done_after_InterProScan_domain_support",
            "sensitivity_only_no_absence_claim",
        )
    if group == "crossdb_confirm":
        return (
            "W2_crossdb_confirm_full_species",
            "next",
            "NCBI_Gene_full_species_then_NCBI_Protein_or_external_orthology_for_no_hits",
            "sensitivity_after_mapping_no_absence_until_sequence_or_external_support",
        )
    if group == "standard_mapping_pending":
        return (
            "W3_standard_mapping_full_species",
            "next_after_W2",
            "NCBI_Gene_full_species_symbol_taxid_mapping_then_low_cost_QC",
            "sensitivity_after_mapping_no_absence_for_low_confidence_or_no_hit_rows",
        )
    return (
        "W9_not_in_full200_sensitivity",
        "hold",
        "excluded_from_current_full200_sensitivity_expansion",
        "do_not_score_without_new_evidence",
    )


def add_batch_ids(gene_plan: pd.DataFrame, genes_per_batch: int = 5) -> pd.DataFrame:
    gene_plan = gene_plan.copy()
    gene_plan["gene_batch_rank"] = ""
    gene_plan["batch_id"] = ""
    gene_plan["batch_priority"] = ""
    todo = gene_plan[gene_plan["species_rows_needed"] > 0].copy()
    for (wave, module), sub in todo.groupby(["expansion_wave", "maintenance_module_v2"], sort=True):
        sub = sub.sort_values(["gene_family_risk", "human_gene_symbol"]).copy()
        for idx, original_index in enumerate(sub.index):
            batch_no = idx // genes_per_batch + 1
            rank = idx + 1
            short_module = module.replace("_", "-")
            batch_id = f"{wave}__{short_module}__batch{batch_no:02d}"
            gene_plan.loc[original_index, "gene_batch_rank"] = str(rank)
            gene_plan.loc[original_index, "batch_id"] = batch_id
            gene_plan.loc[original_index, "batch_priority"] = (
                "run_next" if wave.startswith("W2_") else "run_after_W2"
            )
    gene_plan.loc[gene_plan["species_rows_needed"] == 0, "batch_priority"] = "complete"
    return gene_plan


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eligibility", type=pathlib.Path, required=True)
    parser.add_argument("--species", type=pathlib.Path, required=True)
    parser.add_argument("--priority1-expansion", type=pathlib.Path, required=True)
    parser.add_argument("--gene-plan-output", type=pathlib.Path, required=True)
    parser.add_argument("--queue-output", type=pathlib.Path, required=True)
    parser.add_argument("--batch-summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--docs-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    eligibility = pd.read_csv(args.eligibility, sep="\t")
    species = pd.read_csv(args.species, sep="\t")
    priority1 = pd.read_csv(args.priority1_expansion, sep="\t")

    sensitivity = eligibility[eligibility["sensitivity_score_allowed"].astype(bool)].copy()
    priority1_genes = set(priority1["human_gene_symbol"].dropna())
    species_count = species["scientific_name"].nunique()

    plan_rows = []
    for _, row in sensitivity.iterrows():
        wave, status, route, claim_policy = route_for_group(str(row["v2_scoring_group"]))
        already_expanded = row["v2_scoring_group"] == "strict_ready" or row["human_gene_symbol"] in priority1_genes
        if already_expanded:
            rows_needed = 0
        else:
            rows_needed = species_count
        plan_rows.append(
            {
                "human_gene_symbol": row["human_gene_symbol"],
                "maintenance_module_v2": row["maintenance_module_v2"],
                "submodule_v2": row["submodule_v2"],
                "v2_scoring_group": row["v2_scoring_group"],
                "claim_use": row["claim_use"],
                "gene_family_risk": row["gene_family_risk"],
                "orthology_feasibility_class": row["orthology_feasibility_class"],
                "expansion_wave": wave,
                "expansion_status": status,
                "recommended_route": route,
                "claim_policy_after_mapping": claim_policy,
                "already_species_expanded": already_expanded,
                "species_rows_needed": rows_needed,
            }
        )
    gene_plan = add_batch_ids(pd.DataFrame(plan_rows))

    queue_genes = gene_plan[gene_plan["species_rows_needed"] > 0].copy()
    queue_rows = []
    for _, gene in queue_genes.iterrows():
        for _, sp in species.iterrows():
            queue_rows.append(
                {
                    "expansion_wave": gene["expansion_wave"],
                    "batch_id": gene["batch_id"],
                    "batch_priority": gene["batch_priority"],
                    "gene_batch_rank": gene["gene_batch_rank"],
                    "human_gene_symbol": gene["human_gene_symbol"],
                    "maintenance_module_v2": gene["maintenance_module_v2"],
                    "submodule_v2": gene["submodule_v2"],
                    "v2_scoring_group": gene["v2_scoring_group"],
                    "gene_family_risk": gene["gene_family_risk"],
                    "scientific_name": sp["scientific_name"],
                    "clade": sp["clade"],
                    "flight_status": sp["flight_status"],
                    "species_taxid": sp["species_taxid"],
                    "best_assembly_accession": sp["best_assembly_accession"],
                    "genome_analysis_tier": sp["genome_analysis_tier"],
                    "genome_quality_risk": sp["genome_quality_risk"],
                    "busco_complete": sp["busco_complete"],
                    "query_route": gene["recommended_route"],
                    "claim_policy_after_mapping": gene["claim_policy_after_mapping"],
                    "ncbi_gene_query_term": f'("{gene["human_gene_symbol"]}"[Gene Name] OR "{gene["human_gene_symbol"]}"[Preferred Symbol]) AND txid{sp["species_taxid"]}[Organism]',
                }
            )
    queue = pd.DataFrame(queue_rows)

    batch_summary = (
        gene_plan.groupby(
            ["expansion_wave", "expansion_status", "maintenance_module_v2", "batch_id", "batch_priority"],
            as_index=False,
        )
        .agg(
            genes=("human_gene_symbol", "nunique"),
            species_rows_needed=("species_rows_needed", "sum"),
        )
        .sort_values(["expansion_wave", "maintenance_module_v2", "batch_id"])
    )

    args.gene_plan_output.parent.mkdir(parents=True, exist_ok=True)
    args.queue_output.parent.mkdir(parents=True, exist_ok=True)
    args.batch_summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.docs_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    gene_plan.to_csv(args.gene_plan_output, sep="\t", index=False)
    queue.to_csv(args.queue_output, sep="\t", index=False)
    batch_summary.to_csv(args.batch_summary_output, sep="\t", index=False)

    wave_summary = (
        gene_plan.groupby(["expansion_wave", "expansion_status"], as_index=False)
        .agg(genes=("human_gene_symbol", "nunique"), species_rows_needed=("species_rows_needed", "sum"))
        .sort_values("expansion_wave")
    )
    module_summary = (
        gene_plan.groupby("maintenance_module_v2", as_index=False)
        .agg(
            sensitivity_genes=("human_gene_symbol", "nunique"),
            remaining_rows=("species_rows_needed", "sum"),
        )
        .sort_values("maintenance_module_v2")
    )

    doc_lines = [
        "# Full 200-Gene Sensitivity Expansion Plan",
        "",
        "## Purpose",
        "",
        "Expand the sequence-updated 200-gene sensitivity panel into species-level orthology matrices without blurring strict, sensitivity, and absence-claim boundaries.",
        "",
        "## Current State",
        "",
        f"- Sensitivity genes: {gene_plan['human_gene_symbol'].nunique()}",
        f"- Primary mechanism-panel species: {species_count}",
        f"- Genes already species-expanded: {int(gene_plan['already_species_expanded'].sum())}",
        f"- Genes still needing full-species expansion: {int((gene_plan['species_rows_needed'] > 0).sum())}",
        f"- Species-query rows still needed: {int(gene_plan['species_rows_needed'].sum())}",
        "",
        "## Expansion Waves",
        "",
    ]
    for _, row in wave_summary.iterrows():
        doc_lines.append(
            f"- {row['expansion_wave']} ({row['expansion_status']}): "
            f"{row['genes']} genes, {row['species_rows_needed']} species-query rows"
        )
    doc_lines.extend(["", "## Module Burden", ""])
    for _, row in module_summary.iterrows():
        doc_lines.append(
            f"- {row['maintenance_module_v2']}: {row['sensitivity_genes']} sensitivity genes, "
            f"{row['remaining_rows']} remaining species-query rows"
        )
    doc_lines.extend(["", "## Executable Batches", ""])
    executable_batches = batch_summary[batch_summary["species_rows_needed"] > 0]
    for _, row in executable_batches.iterrows():
        doc_lines.append(
            f"- {row['batch_id']} ({row['batch_priority']}): "
            f"{row['genes']} genes, {row['species_rows_needed']} rows"
        )
    doc_lines.extend(
        [
            "",
            "## Execution Rules",
            "",
            "- W2 crossdb_confirm genes are the next batch. Run full-species NCBI Gene mapping first, then protein/external orthology only for no-hit or non-exact rows.",
            "- W3 standard_mapping_pending genes should be run after W2, preferably module-balanced rather than all at once.",
            "- TDRD/domain-supported genes remain sensitivity-only unless tree/HMM paralog discrimination is added.",
            "- No newly expanded row is allowed to support a biological absence claim without independent sequence/external orthology evidence.",
            "- PGLS using partial expansion must be labelled interim; full200 sensitivity PGLS starts only after W2+W3 species matrices are complete and missingness covariates are attached.",
            "",
            "## Stop/Pause Rules",
            "",
            "- Pause expansion if any wave has repeated NCBI query errors above 5%.",
            "- Pause interpretation if a module has mean species coverage below 50% after W2+W3 mapping.",
            "- Do not escalate to manuscript-level claims if repeat/chromatin effects depend only on low-confidence or no-hit-imputed rows.",
        ]
    )
    args.docs_output.write_text("\n".join(doc_lines) + "\n", encoding="utf-8")

    report_lines = [
        "# Phase 2 Full200 Sensitivity Expansion Plan Report",
        "",
        f"Sensitivity genes: {gene_plan['human_gene_symbol'].nunique()}",
        f"Already species-expanded genes: {int(gene_plan['already_species_expanded'].sum())}",
        f"Genes requiring expansion: {int((gene_plan['species_rows_needed'] > 0).sum())}",
        f"Species-query rows requiring expansion: {int(gene_plan['species_rows_needed'].sum())}",
        "",
        "## Wave Summary",
        "",
    ]
    for _, row in wave_summary.iterrows():
        report_lines.append(f"- {row['expansion_wave']}: {row['genes']} genes, {row['species_rows_needed']} rows")
    args.report.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.gene_plan_output}")


if __name__ == "__main__":
    main()
