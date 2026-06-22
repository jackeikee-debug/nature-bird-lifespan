"""Build a Week 4 queue for cross-database orthology validation."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DEFAULT_MATRIX = pathlib.Path("data/processed/ortholog_matrix_primary_diamond_validated_candidates.tsv")
DEFAULT_SUMMARY = pathlib.Path("results/tables/week4_score_variant_pgls_summary.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/week4_orthology_validation_queue.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week4_orthology_validation_queue_report.md")

MODULE_PRIORITY = {
    "transposon_suppression": 1,
    "DNA_repair": 2,
    "mitochondrial_quality_control": 2,
    "proteostasis": 2,
    "cancer_surveillance": 2,
    "autophagy": 3,
    "inflammation_control": 4,
}


def source_class(status: str) -> str:
    if status == "ncbi_gene_candidate":
        return "ncbi_gene"
    if status == "gff_rescue_candidate":
        return "gff_rescue"
    if status == "diamond_validated_protein_candidate":
        return "diamond_protein"
    return "unresolved"


def validation_reason(row: pd.Series) -> str:
    status = row["final_candidate_status"]
    confidence = row.get("final_candidate_confidence", "")
    module = row["maintenance_module"]
    reasons = []
    if module == "transposon_suppression":
        reasons.append("top_week4_module")
    if status not in {"ncbi_gene_candidate", "gff_rescue_candidate", "diamond_validated_protein_candidate"}:
        reasons.append("unresolved_candidate")
    if status in {"gff_rescue_candidate", "diamond_validated_protein_candidate"}:
        reasons.append("rescue_candidate")
    if confidence in {"", "low", "medium"}:
        reasons.append(f"{confidence or 'unlabeled'}_confidence")
    return ";".join(reasons) if reasons else "source_control_check"


def suggested_validation(row: pd.Series) -> str:
    status = row["final_candidate_status"]
    if status == "diamond_validated_protein_candidate":
        return "reciprocal_BLAST_or_DIAMOND_plus_OMA_OrthoDB_check"
    if status == "gff_rescue_candidate":
        return "protein_sequence_reciprocal_search_plus_external_orthology_check"
    if status == "ncbi_gene_candidate":
        return "cross_check_OMA_OrthoDB_Ensembl_Compara"
    return "search_OMA_OrthoDB_Ensembl_then_targeted_protein_similarity"


def priority_rank(row: pd.Series) -> int:
    module_rank = MODULE_PRIORITY.get(row["maintenance_module"], 5)
    status = row["final_candidate_status"]
    if row["maintenance_module"] == "transposon_suppression" and status != "ncbi_gene_candidate":
        return 1
    if row["maintenance_module"] == "transposon_suppression":
        return 2
    if module_rank == 2 and status != "ncbi_gene_candidate":
        return 3
    if module_rank == 2:
        return 4
    if status != "ncbi_gene_candidate":
        return 5
    return 6


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, default=DEFAULT_MATRIX)
    parser.add_argument("--summary", type=pathlib.Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t")
    summary = pd.read_csv(args.summary, sep="\t")
    priority_by_module = dict(zip(summary["maintenance_module"], summary["week4_priority"], strict=False))

    queue = matrix.copy()
    queue["week4_module_priority"] = queue["maintenance_module"].map(priority_by_module).fillna("not_ranked")
    queue["source_class"] = queue["final_candidate_status"].map(source_class)
    queue["validation_reason"] = queue.apply(validation_reason, axis=1)
    queue["suggested_validation"] = queue.apply(suggested_validation, axis=1)
    queue["priority_rank"] = queue.apply(priority_rank, axis=1)
    queue["validation_batch"] = queue["priority_rank"].map(
        {
            1: "transposon_rescue_or_unresolved",
            2: "transposon_ncbi_crosscheck",
            3: "secondary_rescue_or_unresolved",
            4: "secondary_ncbi_crosscheck",
            5: "exploratory_rescue_or_unresolved",
            6: "exploratory_ncbi_crosscheck",
        }
    )

    keep = [
        "priority_rank",
        "validation_batch",
        "week4_module_priority",
        "maintenance_module",
        "human_gene_symbol",
        "scientific_name",
        "clade",
        "flight_status",
        "genome_analysis_tier",
        "species_taxid",
        "best_assembly_accession",
        "final_candidate_status",
        "final_candidate_source",
        "final_candidate_confidence",
        "ortholog_gene_id",
        "ortholog_gene_symbol",
        "diamond_validation_status",
        "diamond_pident",
        "diamond_query_coverage",
        "diamond_subject_coverage",
        "validation_reason",
        "suggested_validation",
    ]
    queue = queue[keep].sort_values(
        ["priority_rank", "maintenance_module", "human_gene_symbol", "scientific_name"]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    queue.to_csv(args.output, sep="\t", index=False)

    batch_counts = queue["validation_batch"].value_counts().sort_index()
    module_counts = queue.groupby(["priority_rank", "maintenance_module"]).size().reset_index(name="rows")
    lines = [
        "# Week 4 Orthology Validation Queue Report",
        "",
        f"Rows queued: {len(queue)}",
        f"Species: {queue['scientific_name'].nunique()}",
        f"Genes: {queue['human_gene_symbol'].nunique()}",
        "",
        "## Rows by Batch",
        "",
    ]
    for batch, count in batch_counts.items():
        lines.append(f"- {batch}: {count}")
    lines.extend(["", "## Rows by Module Priority", ""])
    for _, rec in module_counts.iterrows():
        lines.append(f"- rank {rec['priority_rank']} / {rec['maintenance_module']}: {rec['rows']}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The first validation batch should start with transposon-suppression rescue or unresolved rows, then cross-check direct NCBI Gene transposon candidates. Secondary batches cover DNA repair, mitochondrial quality control, proteostasis, and cancer surveillance.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
