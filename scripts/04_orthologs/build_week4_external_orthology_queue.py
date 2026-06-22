"""Build an external orthology queue for ambiguous Week 4 transposon rows."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DEFAULT_VALIDATION = pathlib.Path("data/processed/week4_transposon_reciprocal_validation.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/week4_external_orthology_queue.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week4_external_orthology_queue_report.md")

QUEUE_STATUSES = {
    "weak_forward_support",
    "not_reciprocal",
    "reciprocal_weak",
    "not_validated",
}


def priority(row: pd.Series) -> int:
    gene = row["human_gene_symbol"]
    status = row["week4_diamond_status"]
    if gene in {"PIWIL2", "TRIM28"} and status in {"not_reciprocal", "weak_forward_support"}:
        return 1
    if status in {"not_reciprocal", "not_validated"}:
        return 2
    if status == "reciprocal_weak":
        return 3
    return 4


def suggested_sources(row: pd.Series) -> str:
    gene = row["human_gene_symbol"]
    if gene.startswith("PIWIL"):
        return "OMA;OrthoDB;Ensembl_Compara;manual_PIWI_family_tree"
    if gene == "TRIM28":
        return "OMA;OrthoDB;Ensembl_Compara;domain_architecture_check"
    return "OMA;OrthoDB;Ensembl_Compara"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validation", type=pathlib.Path, default=DEFAULT_VALIDATION)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    rows = pd.read_csv(args.validation, sep="\t")
    queue = rows[rows["week4_diamond_status"].isin(QUEUE_STATUSES)].copy()
    queue["external_validation_priority"] = queue.apply(priority, axis=1)
    queue["suggested_external_sources"] = queue.apply(suggested_sources, axis=1)
    queue["external_validation_question"] = queue.apply(
        lambda row: (
            f"Does {row['scientific_name']} {row['human_gene_symbol']} candidate "
            f"{row.get('forward_target_protein_id', '')} resolve to the expected ortholog?"
        ),
        axis=1,
    )
    keep = [
        "external_validation_priority",
        "human_gene_symbol",
        "scientific_name",
        "clade",
        "flight_status",
        "genome_analysis_tier",
        "species_taxid",
        "best_assembly_accession",
        "final_candidate_status",
        "week4_diamond_status",
        "week4_diamond_reason",
        "forward_target_protein_id",
        "forward_pident",
        "forward_qcovhsp",
        "forward_scovhsp",
        "forward_bitscore",
        "reciprocal_gene",
        "reciprocal_top_reference",
        "suggested_external_sources",
        "external_validation_question",
    ]
    queue = queue[keep].sort_values(
        ["external_validation_priority", "human_gene_symbol", "scientific_name"]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    queue.to_csv(args.output, sep="\t", index=False)

    status_counts = queue["week4_diamond_status"].value_counts().sort_index()
    gene_counts = queue["human_gene_symbol"].value_counts().sort_index()
    lines = [
        "# Week 4 External Orthology Queue Report",
        "",
        f"Rows queued: {len(queue)}",
        f"Species: {queue['scientific_name'].nunique()}",
        "",
        "## Rows by DIAMOND Status",
        "",
    ]
    for status, count in status_counts.items():
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Rows by Gene", ""])
    for gene, count in gene_counts.items():
        lines.append(f"- {gene}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "These rows should be checked against OMA, OrthoDB, Ensembl Compara, or manual family/domain evidence before they are counted as positive transposon-suppression orthologs.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
