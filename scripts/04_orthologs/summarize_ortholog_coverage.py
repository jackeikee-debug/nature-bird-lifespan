"""Summarize candidate ortholog coverage by module, clade, species, and gene."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import defaultdict


DEFAULT_INPUT = pathlib.Path("data/processed/ortholog_matrix_primary_ncbi_gene_candidates.tsv")
DEFAULT_PREFIX = pathlib.Path("results/tables/ortholog_coverage_primary")
DEFAULT_REPORT = pathlib.Path("results/reports/ortholog_coverage_primary_report.md")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def is_found(row: dict[str, str], status_field: str, found_statuses: set[str]) -> bool:
    return row.get(status_field) in found_statuses


def group_summary(
    rows: list[dict[str, str]],
    group_fields: list[str],
    status_field: str,
    found_statuses: set[str],
) -> list[dict[str, str]]:
    grouped: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[tuple(row[field] for field in group_fields)].append(row)

    summary = []
    for key, group_rows in grouped.items():
        queried = [
            row
            for row in group_rows
            if row.get(status_field) != "not_queried_max_queries_reached"
        ]
        found = [row for row in queried if is_found(row, status_field, found_statuses)]
        missing = [
            row
            for row in queried
            if row.get(status_field) in {"no_ncbi_gene_candidate", "not_found"}
        ]
        errors = [row for row in queried if row.get(status_field) == "query_error"]
        denominator = len(queried)
        coverage = len(found) / denominator if denominator else 0.0
        out = {field: value for field, value in zip(group_fields, key)}
        out.update(
            {
                "rows_total": str(len(group_rows)),
                "rows_queried": str(denominator),
                "candidate_found": str(len(found)),
                "no_candidate": str(len(missing)),
                "query_error": str(len(errors)),
                "coverage_fraction": f"{coverage:.6f}",
            }
        )
        summary.append(out)
    summary.sort(key=lambda item: tuple(item[field] for field in group_fields))
    return summary


def write_report(
    path: pathlib.Path,
    rows: list[dict[str, str]],
    module_clade: list[dict[str, str]],
    species: list[dict[str, str]],
    genes: list[dict[str, str]],
    status_field: str,
    found_statuses: set[str],
) -> None:
    queried = [row for row in rows if row.get(status_field) != "not_queried_max_queries_reached"]
    found = [row for row in queried if is_found(row, status_field, found_statuses)]
    coverage = len(found) / len(queried) if queried else 0.0
    low_species = [
        row
        for row in species
        if int(row["rows_queried"]) > 0 and float(row["coverage_fraction"]) < 0.75
    ]
    low_genes = [
        row
        for row in genes
        if int(row["rows_queried"]) > 0 and float(row["coverage_fraction"]) < 0.75
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Ortholog Coverage Report",
                "",
                f"Input rows: {len(rows)}",
                f"Queried rows: {len(queried)}",
                f"Candidate hits: {len(found)}",
                f"Overall candidate coverage: {coverage:.3f}",
                f"Status field: {status_field}",
                f"Found statuses: {', '.join(sorted(found_statuses))}",
                "",
                "## Low-Coverage Species",
                *(
                    [
                        f"- {row['scientific_name']}: {row['coverage_fraction']} ({row['candidate_found']}/{row['rows_queried']})"
                        for row in low_species[:20]
                    ]
                    or ["- none below threshold"]
                ),
                "",
                "## Low-Coverage Genes",
                *(
                    [
                        f"- {row['human_gene_symbol']}: {row['coverage_fraction']} ({row['candidate_found']}/{row['rows_queried']})"
                        for row in low_genes[:20]
                    ]
                    or ["- none below threshold"]
                ),
                "",
                "## Interpretation",
                "Coverage is based on the configured candidate status field and found statuses. It is useful for feasibility triage but should not be interpreted as final ortholog conservation without cross-database validation.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-prefix", type=pathlib.Path, default=DEFAULT_PREFIX)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    parser.add_argument("--status-field", default="ortholog_query_status")
    parser.add_argument("--found-statuses", default="candidate_found")
    args = parser.parse_args()

    rows = read_tsv(args.input)
    found_statuses = {status.strip() for status in args.found_statuses.split(",") if status.strip()}
    module_clade = group_summary(
        rows, ["genome_panel_version", "maintenance_module", "clade"], args.status_field, found_statuses
    )
    species = group_summary(
        rows, ["genome_panel_version", "scientific_name", "clade"], args.status_field, found_statuses
    )
    genes = group_summary(
        rows, ["genome_panel_version", "human_gene_symbol", "maintenance_module"], args.status_field, found_statuses
    )

    write_tsv(
        args.output_prefix.with_name(args.output_prefix.name + "_by_module_clade.tsv"),
        module_clade,
        list(module_clade[0].keys()) if module_clade else [],
    )
    write_tsv(
        args.output_prefix.with_name(args.output_prefix.name + "_by_species.tsv"),
        species,
        list(species[0].keys()) if species else [],
    )
    write_tsv(
        args.output_prefix.with_name(args.output_prefix.name + "_by_gene.tsv"),
        genes,
        list(genes[0].keys()) if genes else [],
    )
    write_report(args.report, rows, module_clade, species, genes, args.status_field, found_statuses)
    print(f"Wrote coverage summaries with prefix {args.output_prefix} and {args.report}")


if __name__ == "__main__":
    main()
