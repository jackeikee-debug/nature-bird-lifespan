"""Build a protein FASTA batch for Phase 2 priority-1 domain confirmation."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DOMAIN_RULES = {
    "TDRD1": "tudor_domain_required",
    "TDRD5": "tudor_domain_required",
    "TDRD6": "tudor_domain_required",
    "TDRD7": "tudor_domain_required",
    "TDRD9": "tudor_plus_helicase_expected",
    "TDRD12": "tudor_plus_helicase_expected",
    "DDX4": "dead_box_helicase_required",
    "GTSF1": "zinc_knuckle_required",
    "TDRKH": "tudor_or_kh_domain_required",
}


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-proteins", type=pathlib.Path, required=True)
    parser.add_argument("--domain-queue", type=pathlib.Path, required=True)
    parser.add_argument("--manual-queue", type=pathlib.Path, required=True)
    parser.add_argument("--fasta-output", type=pathlib.Path, required=True)
    parser.add_argument("--manifest-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    proteins = pd.read_csv(args.candidate_proteins, sep="\t")
    domain = pd.read_csv(args.domain_queue, sep="\t")
    manual = pd.read_csv(args.manual_queue, sep="\t")

    target_genes = sorted(
        set(domain["human_gene_symbol"].dropna()) | set(manual["human_gene_symbol"].dropna())
    )
    batch = proteins[
        proteins["human_gene_symbol"].isin(target_genes)
        & (proteins["protein_fetch_status"] == "protein_sequence_found")
    ].copy()
    batch["domain_rule"] = batch["human_gene_symbol"].map(DOMAIN_RULES).fillna("manual_domain_review")
    batch["domain_batch_record_id"] = (
        batch["human_gene_symbol"].astype(str)
        + "|"
        + batch["scientific_name"].astype(str).str.replace(" ", "_", regex=False)
        + "|"
        + batch["protein_accession"].astype(str)
        + "|rank:"
        + batch["protein_rank_for_gene_species"].astype(str)
    )
    batch = batch.sort_values(
        ["human_gene_symbol", "scientific_name", "protein_rank_for_gene_species", "protein_accession"]
    )

    args.fasta_output.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    with args.fasta_output.open("w", encoding="utf-8", newline="\n") as handle:
        for _, row in batch.iterrows():
            header = row["domain_batch_record_id"]
            desc = str(row.get("fasta_header", "")).replace("\n", " ")
            handle.write(f">{header} {desc}\n{wrap(str(row['sequence']))}\n")

    manifest_cols = [
        "domain_batch_record_id",
        "human_gene_symbol",
        "scientific_name",
        "species_taxid",
        "best_assembly_accession",
        "ncbi_gene_id",
        "protein_accession",
        "protein_rank_for_gene_species",
        "protein_length",
        "domain_rule",
        "combined_strict_upgrade_decision",
        "recommended_next_step",
        "domain_or_sequence_rule",
        "fasta_header",
    ]
    batch[manifest_cols].to_csv(args.manifest_output, sep="\t", index=False)

    counts = batch.groupby(["human_gene_symbol", "domain_rule"]).size().reset_index(name="protein_records")
    lines = [
        "# Phase 2 Priority 1 Domain Confirmation Batch Report",
        "",
        "## Summary",
        "",
        f"Target genes: {len(target_genes)}",
        f"Protein records written: {len(batch)}",
        f"Species represented: {batch['scientific_name'].nunique()}",
        "",
        "## Gene Counts",
        "",
    ]
    for _, row in counts.iterrows():
        lines.append(f"- {row['human_gene_symbol']} ({row['domain_rule']}): {row['protein_records']}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This batch is designed for InterProScan/Pfam domain validation of high-paralog or manual-review priority-1 genes. Domain support should be interpreted as protein-family evidence, not by itself as exact orthology for TDRD paralogs.",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.fasta_output}")


if __name__ == "__main__":
    main()
