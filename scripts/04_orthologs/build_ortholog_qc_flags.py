"""Build QC flags from ortholog coverage summaries."""

from __future__ import annotations

import argparse
import csv
import pathlib


DEFAULT_SPECIES = pathlib.Path("results/tables/ortholog_coverage_primary_by_species.tsv")
DEFAULT_GENE = pathlib.Path("results/tables/ortholog_coverage_primary_by_gene.tsv")
DEFAULT_MODULE_CLADE = pathlib.Path("results/tables/ortholog_coverage_primary_by_module_clade.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/ortholog_mapping_qc_flags.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/ortholog_mapping_qc_flags_report.md")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "flag_type",
        "priority",
        "entity",
        "clade",
        "maintenance_module",
        "coverage_fraction",
        "candidate_found",
        "rows_queried",
        "recommended_action",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def coverage(row: dict[str, str]) -> float:
    return float(row["coverage_fraction"])


def build_flags(
    species_rows: list[dict[str, str]],
    gene_rows: list[dict[str, str]],
    module_clade_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    flags = []
    for row in species_rows:
        cov = coverage(row)
        if cov == 0.0:
            priority = "high"
            action = "do_not_interpret_as_gene_loss; recover from annotation GFF/protein FASTA or cross-database ortholog source"
        elif cov < 0.75:
            priority = "medium"
            action = "manual taxonomy and annotation-source review before mechanism scoring"
        else:
            continue
        flags.append(
            {
                "flag_type": "species_low_coverage",
                "priority": priority,
                "entity": row["scientific_name"],
                "clade": row["clade"],
                "maintenance_module": "all",
                "coverage_fraction": row["coverage_fraction"],
                "candidate_found": row["candidate_found"],
                "rows_queried": row["rows_queried"],
                "recommended_action": action,
            }
        )

    for row in gene_rows:
        cov = coverage(row)
        if cov < 0.5:
            priority = "high"
        elif cov < 0.75:
            priority = "medium"
        else:
            continue
        flags.append(
            {
                "flag_type": "gene_low_coverage",
                "priority": priority,
                "entity": row["human_gene_symbol"],
                "clade": "all",
                "maintenance_module": row["maintenance_module"],
                "coverage_fraction": row["coverage_fraction"],
                "candidate_found": row["candidate_found"],
                "rows_queried": row["rows_queried"],
                "recommended_action": "check symbol aliases and cross-database ortholog support before treating missingness as biological",
            }
        )

    for row in module_clade_rows:
        cov = coverage(row)
        if cov < 0.5:
            priority = "high"
        elif cov < 0.75:
            priority = "medium"
        else:
            continue
        flags.append(
            {
                "flag_type": "module_clade_low_coverage",
                "priority": priority,
                "entity": f"{row['maintenance_module']}|{row['clade']}",
                "clade": row["clade"],
                "maintenance_module": row["maintenance_module"],
                "coverage_fraction": row["coverage_fraction"],
                "candidate_found": row["candidate_found"],
                "rows_queried": row["rows_queried"],
                "recommended_action": "avoid module-level biological interpretation until missing annotations are rescued",
            }
        )

    flags.sort(key=lambda item: (item["priority"] != "high", item["flag_type"], float(item["coverage_fraction"])))
    return flags


def write_report(path: pathlib.Path, flags: list[dict[str, str]]) -> None:
    high = [row for row in flags if row["priority"] == "high"]
    medium = [row for row in flags if row["priority"] == "medium"]
    zero_species = [
        row for row in flags if row["flag_type"] == "species_low_coverage" and row["coverage_fraction"] == "0.000000"
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Ortholog Mapping QC Flags Report",
                "",
                f"QC flags: {len(flags)}",
                f"High-priority flags: {len(high)}",
                f"Medium-priority flags: {len(medium)}",
                f"Zero-coverage species: {len(zero_species)}",
                "",
                "## Interpretation",
                "Low NCBI Gene candidate coverage should be treated primarily as a database/annotation retrieval problem at this stage. It is not evidence for gene loss unless confirmed by assembly-level searches and independent ortholog databases.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", type=pathlib.Path, default=DEFAULT_SPECIES)
    parser.add_argument("--gene", type=pathlib.Path, default=DEFAULT_GENE)
    parser.add_argument("--module-clade", type=pathlib.Path, default=DEFAULT_MODULE_CLADE)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    flags = build_flags(read_tsv(args.species), read_tsv(args.gene), read_tsv(args.module_clade))
    write_tsv(args.output, flags)
    write_report(args.report, flags)
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
