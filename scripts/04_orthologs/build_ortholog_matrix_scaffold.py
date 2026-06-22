"""Create primary and sensitivity ortholog-matrix scaffolds.

The scaffold intentionally records query targets without inventing ortholog calls.
Later scripts can fill ortholog identifiers from NCBI, Ensembl Compara, OMA, or
OrthoDB while preserving the same row keys.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import Counter


DEFAULT_GENE_SET = pathlib.Path("data/processed/maintenance_gene_sets.tsv")
DEFAULT_PRIMARY_PANEL = pathlib.Path("data/processed/genome_maintenance_species_primary.tsv")
DEFAULT_SENSITIVITY_PANEL = pathlib.Path("data/processed/genome_maintenance_species_sensitivity.tsv")
DEFAULT_PRIMARY_OUTPUT = pathlib.Path("data/processed/ortholog_matrix_primary.tsv")
DEFAULT_SENSITIVITY_OUTPUT = pathlib.Path("data/processed/ortholog_matrix_sensitivity.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/ortholog_matrix_scaffold_report.md")


OUTPUT_FIELDS = [
    "genome_panel_version",
    "scientific_name",
    "clade",
    "flight_status",
    "maintenance_module",
    "human_gene_symbol",
    "source_set",
    "species_taxid",
    "best_assembly_accession",
    "genome_analysis_tier",
    "ortholog_query_status",
    "ortholog_gene_id",
    "ortholog_gene_symbol",
    "ortholog_status",
    "ortholog_source_database",
    "ortholog_source_url",
    "ortholog_confidence",
    "copy_number_estimate",
    "notes",
]


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def scaffold_rows(
    species_rows: list[dict[str, str]],
    gene_rows: list[dict[str, str]],
    panel_version: str,
) -> list[dict[str, str]]:
    rows = []
    for species in species_rows:
        for gene in gene_rows:
            rows.append(
                {
                    "genome_panel_version": panel_version,
                    "scientific_name": species["scientific_name"],
                    "clade": species["clade"],
                    "flight_status": species["flight_status"],
                    "maintenance_module": gene["maintenance_module"],
                    "human_gene_symbol": gene["human_gene_symbol"],
                    "source_set": gene["source_set"],
                    "species_taxid": species["species_taxid"],
                    "best_assembly_accession": species["best_assembly_accession"],
                    "genome_analysis_tier": species["genome_analysis_tier"],
                    "ortholog_query_status": "pending",
                    "ortholog_gene_id": "",
                    "ortholog_gene_symbol": "",
                    "ortholog_status": "",
                    "ortholog_source_database": "",
                    "ortholog_source_url": "",
                    "ortholog_confidence": "",
                    "copy_number_estimate": "",
                    "notes": "",
                }
            )
    rows.sort(
        key=lambda row: (
            row["genome_panel_version"],
            row["clade"],
            row["scientific_name"],
            row["maintenance_module"],
            row["human_gene_symbol"],
        )
    )
    return rows


def write_report(
    path: pathlib.Path,
    gene_rows: list[dict[str, str]],
    primary_species: list[dict[str, str]],
    sensitivity_species: list[dict[str, str]],
    primary_rows: list[dict[str, str]],
    sensitivity_rows: list[dict[str, str]],
) -> None:
    module_counts = Counter(row["maintenance_module"] for row in gene_rows)
    primary_clades = Counter(row["clade"] for row in primary_species)
    sensitivity_clades = Counter(row["clade"] for row in sensitivity_species)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Ortholog Matrix Scaffold Report",
                "",
                f"Seed genes: {len(gene_rows)}",
                f"Primary species: {len(primary_species)}",
                f"Sensitivity species: {len(sensitivity_species)}",
                f"Primary scaffold rows: {len(primary_rows)}",
                f"Sensitivity scaffold rows: {len(sensitivity_rows)}",
                "",
                "## Gene Modules",
                *[f"- {module}: {count}" for module, count in sorted(module_counts.items())],
                "",
                "## Primary Species by Clade",
                *[f"- {clade}: {count}" for clade, count in sorted(primary_clades.items())],
                "",
                "## Sensitivity Species by Clade",
                *[f"- {clade}: {count}" for clade, count in sorted(sensitivity_clades.items())],
                "",
                "## Interpretation",
                "These matrices are query scaffolds, not completed ortholog calls. All ortholog fields are pending until filled by a source-specific mapper.",
                "",
                "Primary matrix rows should support the first annotated-genome mechanism analysis. Sensitivity matrix rows should test whether broad patterns remain similar when Tier 3 assembly-only genomes are included.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-set", type=pathlib.Path, default=DEFAULT_GENE_SET)
    parser.add_argument("--primary-panel", type=pathlib.Path, default=DEFAULT_PRIMARY_PANEL)
    parser.add_argument("--sensitivity-panel", type=pathlib.Path, default=DEFAULT_SENSITIVITY_PANEL)
    parser.add_argument("--primary-output", type=pathlib.Path, default=DEFAULT_PRIMARY_OUTPUT)
    parser.add_argument("--sensitivity-output", type=pathlib.Path, default=DEFAULT_SENSITIVITY_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    genes = read_tsv(args.gene_set)
    primary_species = read_tsv(args.primary_panel)
    sensitivity_species = read_tsv(args.sensitivity_panel)
    primary_rows = scaffold_rows(primary_species, genes, "primary")
    sensitivity_rows = scaffold_rows(sensitivity_species, genes, "sensitivity")
    write_tsv(args.primary_output, primary_rows, OUTPUT_FIELDS)
    write_tsv(args.sensitivity_output, sensitivity_rows, OUTPUT_FIELDS)
    write_report(args.report, genes, primary_species, sensitivity_species, primary_rows, sensitivity_rows)
    print(f"Wrote {args.primary_output}, {args.sensitivity_output}, and {args.report}")


if __name__ == "__main__":
    main()
