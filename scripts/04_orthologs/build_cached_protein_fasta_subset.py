"""Build a cached-FASTA subset for validation when full downloads are incomplete."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DEFAULT_BATCH = pathlib.Path("data/processed/week4_transposon_ncbi_crosscheck_batch.tsv")
DEFAULT_MANIFEST = pathlib.Path("data/processed/week4_transposon_ncbi_crosscheck_protein_fasta_manifest.tsv")
DEFAULT_FASTA_DIR = pathlib.Path("data/raw/week4_protein_fastas")
DEFAULT_BATCH_OUTPUT = pathlib.Path("data/processed/week4_transposon_ncbi_crosscheck_cached_batch.tsv")
DEFAULT_LOG_OUTPUT = pathlib.Path("data/interim/week4_transposon_ncbi_crosscheck_cached_fasta_log.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week4_transposon_ncbi_crosscheck_cached_subset_report.md")


def local_fasta_path(fasta_dir: pathlib.Path, accession: str, protein_url: str) -> pathlib.Path:
    basename = protein_url.rstrip("/").split("/")[-1]
    return fasta_dir / accession / basename


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=pathlib.Path, default=DEFAULT_BATCH)
    parser.add_argument("--manifest", type=pathlib.Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--fasta-dir", type=pathlib.Path, default=DEFAULT_FASTA_DIR)
    parser.add_argument("--batch-output", type=pathlib.Path, default=DEFAULT_BATCH_OUTPUT)
    parser.add_argument("--log-output", type=pathlib.Path, default=DEFAULT_LOG_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    batch = pd.read_csv(args.batch, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    log_rows = []
    cached_accessions = set()
    for _, row in manifest.iterrows():
        path = local_fasta_path(args.fasta_dir, row["best_assembly_accession"], row["protein_url"])
        status = "cached" if path.exists() and path.stat().st_size > 0 else "missing"
        if status == "cached":
            cached_accessions.add(row["best_assembly_accession"])
        log_rows.append(
            {
                "scientific_name": row["scientific_name"],
                "clade": row["clade"],
                "species_taxid": row["species_taxid"],
                "best_assembly_accession": row["best_assembly_accession"],
                "protein_url": row["protein_url"],
                "protein_local_path": str(path),
                "download_status": status,
                "file_size_bytes": str(path.stat().st_size if path.exists() else 0),
                "error": "",
            }
        )
    cached_batch = batch[batch["best_assembly_accession"].isin(cached_accessions)].copy()
    args.batch_output.parent.mkdir(parents=True, exist_ok=True)
    cached_batch.to_csv(args.batch_output, sep="\t", index=False)
    args.log_output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(log_rows).to_csv(args.log_output, sep="\t", index=False)

    gene_counts = cached_batch["human_gene_symbol"].value_counts().sort_index()
    lines = [
        "# Cached Protein FASTA Subset Report",
        "",
        f"Manifest rows: {len(manifest)}",
        f"Cached assemblies: {len(cached_accessions)}",
        f"Cached batch rows: {len(cached_batch)}",
        "",
        "## Cached Batch Rows by Gene",
        "",
    ]
    for gene, count in gene_counts.items():
        lines.append(f"- {gene}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This subset is intended for pilot sequence validation while the full protein FASTA download continues or is retried.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.batch_output}, {args.log_output}, and {args.report}")


if __name__ == "__main__":
    main()
