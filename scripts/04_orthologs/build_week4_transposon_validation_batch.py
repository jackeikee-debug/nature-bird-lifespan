"""Build the first Week 4 transposon orthology-validation batch."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DEFAULT_QUEUE = pathlib.Path("data/processed/week4_orthology_validation_queue.tsv")
DEFAULT_GENOME_AUDIT = pathlib.Path("data/processed/genome_availability_audit.tsv")
DEFAULT_RESCUE_MANIFEST = pathlib.Path("data/processed/annotation_rescue_manifest.tsv")
DEFAULT_BATCH = pathlib.Path("data/processed/week4_transposon_validation_batch.tsv")
DEFAULT_FASTA_MANIFEST = pathlib.Path("data/processed/week4_transposon_protein_fasta_manifest.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week4_transposon_validation_batch_report.md")


def https_url(url: str) -> str:
    if not isinstance(url, str) or not url:
        return ""
    return url.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")


def protein_url_from_ftp(ftp_path: str) -> str:
    ftp_path = https_url(ftp_path).rstrip("/")
    if not ftp_path:
        return ""
    basename = ftp_path.split("/")[-1]
    return f"{ftp_path}/{basename}_protein.faa.gz"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=pathlib.Path, default=DEFAULT_QUEUE)
    parser.add_argument("--genome-audit", type=pathlib.Path, default=DEFAULT_GENOME_AUDIT)
    parser.add_argument("--rescue-manifest", type=pathlib.Path, default=DEFAULT_RESCUE_MANIFEST)
    parser.add_argument("--batch-output", type=pathlib.Path, default=DEFAULT_BATCH)
    parser.add_argument("--fasta-manifest-output", type=pathlib.Path, default=DEFAULT_FASTA_MANIFEST)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    parser.add_argument("--validation-batch", default="transposon_rescue_or_unresolved")
    args = parser.parse_args()

    queue = pd.read_csv(args.queue, sep="\t")
    genome = pd.read_csv(args.genome_audit, sep="\t")
    rescue = pd.read_csv(args.rescue_manifest, sep="\t")

    batch = queue[queue["validation_batch"] == args.validation_batch].copy()
    rescue_urls = dict(zip(rescue["best_assembly_accession"], rescue["protein_url"], strict=False))

    genome_cols = [
        "scientific_name",
        "best_assembly_accession",
        "ftp_path_refseq",
        "ftp_path_genbank",
        "assembly_name",
        "assembly_status",
        "annotation_report_url",
    ]
    batch = batch.merge(
        genome[genome_cols],
        on=["scientific_name", "best_assembly_accession"],
        how="left",
    )

    protein_urls = []
    protein_url_sources = []
    for _, row in batch.iterrows():
        accession = row["best_assembly_accession"]
        if accession in rescue_urls and isinstance(rescue_urls[accession], str) and rescue_urls[accession]:
            protein_urls.append(rescue_urls[accession])
            protein_url_sources.append("annotation_rescue_manifest")
            continue
        refseq_url = protein_url_from_ftp(row.get("ftp_path_refseq", ""))
        genbank_url = protein_url_from_ftp(row.get("ftp_path_genbank", ""))
        if refseq_url:
            protein_urls.append(refseq_url)
            protein_url_sources.append("genome_audit_refseq_ftp")
        elif genbank_url:
            protein_urls.append(genbank_url)
            protein_url_sources.append("genome_audit_genbank_ftp")
        else:
            protein_urls.append("")
            protein_url_sources.append("missing_ftp_path")
    batch["protein_url"] = protein_urls
    batch["protein_url_source"] = protein_url_sources

    fasta_manifest = (
        batch[
            [
                "scientific_name",
                "clade",
                "flight_status",
                "species_taxid",
                "best_assembly_accession",
                "genome_analysis_tier",
                "protein_url",
                "protein_url_source",
            ]
        ]
        .drop_duplicates()
        .sort_values(["clade", "scientific_name", "best_assembly_accession"])
    )
    fasta_manifest["download_priority"] = fasta_manifest["protein_url"].apply(
        lambda value: "download" if isinstance(value, str) and value else "missing_url"
    )

    args.batch_output.parent.mkdir(parents=True, exist_ok=True)
    batch.to_csv(args.batch_output, sep="\t", index=False)
    args.fasta_manifest_output.parent.mkdir(parents=True, exist_ok=True)
    fasta_manifest.to_csv(args.fasta_manifest_output, sep="\t", index=False)

    gene_counts = batch["human_gene_symbol"].value_counts().sort_index()
    status_counts = batch["final_candidate_status"].value_counts().sort_index()
    url_counts = fasta_manifest["download_priority"].value_counts().sort_index()
    lines = [
        "# Week 4 Transposon Validation Batch Report",
        "",
        f"Validation batch: {args.validation_batch}",
        f"Batch rows: {len(batch)}",
        f"Species: {batch['scientific_name'].nunique()}",
        f"Assemblies: {batch['best_assembly_accession'].nunique()}",
        f"Protein FASTA manifest rows: {len(fasta_manifest)}",
        "",
        "## Rows by Gene",
        "",
    ]
    for gene, count in gene_counts.items():
        lines.append(f"- {gene}: {count}")
    lines.extend(["", "## Rows by Candidate Status", ""])
    for status, count in status_counts.items():
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Protein FASTA URL Status", ""])
    for status, count in url_counts.items():
        lines.append(f"- {status}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This Week 4 validation batch should be checked before making mechanism claims. Rescue/unresolved batches test missingness and paralogy risk; NCBI cross-check batches test whether direct database candidates are robust to sequence-level validation.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.batch_output}, {args.fasta_manifest_output}, and {args.report}")


if __name__ == "__main__":
    main()
