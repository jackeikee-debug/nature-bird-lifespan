"""Build a manifest for assembly-level annotation rescue of low-coverage species."""

from __future__ import annotations

import argparse
import csv
import pathlib


DEFAULT_PANEL = pathlib.Path("data/processed/genome_maintenance_species_primary.tsv")
DEFAULT_FLAGS = pathlib.Path("data/processed/ortholog_mapping_qc_flags.tsv")
DEFAULT_GENES = pathlib.Path("data/processed/maintenance_gene_sets.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/annotation_rescue_manifest.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/annotation_rescue_manifest_report.md")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    fields = list(rows[0].keys()) if rows else [
        "scientific_name",
        "clade",
        "best_assembly_accession",
        "assembly_source",
        "gff_url",
        "protein_url",
        "genes_to_rescue",
        "rescue_reason",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def https_ftp_path(row: dict[str, str]) -> tuple[str, str]:
    if row.get("ftp_path_refseq"):
        return row["ftp_path_refseq"].replace("ftp://", "https://"), "refseq"
    if row.get("ftp_path_genbank"):
        return row["ftp_path_genbank"].replace("ftp://", "https://"), "genbank"
    return "", "missing_ftp"


def file_urls(ftp_path: str) -> tuple[str, str]:
    if not ftp_path:
        return "", ""
    basename = ftp_path.rstrip("/").split("/")[-1]
    return f"{ftp_path}/{basename}_genomic.gff.gz", f"{ftp_path}/{basename}_protein.faa.gz"


def rescue_species(flags: list[dict[str, str]]) -> set[str]:
    return {
        row["entity"]
        for row in flags
        if row["flag_type"] == "species_low_coverage"
        and row["priority"] == "high"
        and row["coverage_fraction"] == "0.000000"
    }


def genes_to_rescue(flags: list[dict[str, str]], genes: list[dict[str, str]]) -> str:
    return ",".join(row["human_gene_symbol"] for row in genes)


def build_manifest(panel: list[dict[str, str]], flags: list[dict[str, str]], genes: list[dict[str, str]]) -> list[dict[str, str]]:
    species = rescue_species(flags)
    priority_genes = genes_to_rescue(flags, genes)
    rows = []
    for row in panel:
        if row["scientific_name"] not in species:
            continue
        ftp_path, source = https_ftp_path(row)
        gff_url, protein_url = file_urls(ftp_path)
        rows.append(
            {
                "scientific_name": row["scientific_name"],
                "clade": row["clade"],
                "flight_status": row["flight_status"],
                "best_assembly_accession": row["best_assembly_accession"],
                "assembly_source": source,
                "ncbi_species_name": row["ncbi_species_name"],
                "species_taxid": row["species_taxid"],
                "gff_url": gff_url,
                "protein_url": protein_url,
                "genes_to_rescue": priority_genes,
                "rescue_reason": "zero_ncbi_gene_candidate_coverage",
                "recommended_first_pass": "download_gff_and_search_gene_attributes",
            }
        )
    rows.sort(key=lambda item: (item["clade"], item["scientific_name"]))
    return rows


def write_report(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    missing = [row for row in rows if not row["gff_url"]]
    genbank = [row for row in rows if row["assembly_source"] == "genbank"]
    refseq = [row for row in rows if row["assembly_source"] == "refseq"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Annotation Rescue Manifest Report",
                "",
                f"Species queued for annotation rescue: {len(rows)}",
                f"GenBank assemblies: {len(genbank)}",
                f"RefSeq assemblies: {len(refseq)}",
                f"Missing FTP paths: {len(missing)}",
                "",
                "## Interpretation",
                "These species had zero NCBI Gene symbol/taxid candidate coverage. The next rescue step searches assembly annotation files directly, because zero NCBI Gene coverage should not be interpreted as gene loss.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", type=pathlib.Path, default=DEFAULT_PANEL)
    parser.add_argument("--flags", type=pathlib.Path, default=DEFAULT_FLAGS)
    parser.add_argument("--genes", type=pathlib.Path, default=DEFAULT_GENES)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    rows = build_manifest(read_tsv(args.panel), read_tsv(args.flags), read_tsv(args.genes))
    write_tsv(args.output, rows)
    write_report(args.report, rows)
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
