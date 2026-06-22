"""Build a normalized seed table of human genome-maintenance genes."""

from __future__ import annotations

import argparse
import csv
import pathlib
import re


DEFAULT_INPUT = pathlib.Path("config/pathway_sets.yaml")
DEFAULT_OUTPUT = pathlib.Path("data/processed/maintenance_gene_sets.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/maintenance_gene_sets_report.md")


MODULE_RE = re.compile(r"^  ([A-Za-z0-9_]+):\s*$")
GENES_RE = re.compile(r"^    example_genes:\s*\[(.*)\]\s*$")


def parse_seed_yaml(path: pathlib.Path) -> list[dict[str, str]]:
    rows = []
    current_module = ""
    module_order = 0
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            module_match = MODULE_RE.match(line)
            if module_match:
                current_module = module_match.group(1)
                module_order += 1
                continue
            genes_match = GENES_RE.match(line)
            if not genes_match or not current_module:
                continue
            genes = [gene.strip() for gene in genes_match.group(1).split(",") if gene.strip()]
            for gene_order, gene in enumerate(genes, start=1):
                rows.append(
                    {
                        "human_gene_symbol": gene.upper(),
                        "maintenance_module": current_module,
                        "source_set": "curated_seed_v0",
                        "source_field": "config/pathway_sets.yaml:example_genes",
                        "module_order": str(module_order),
                        "gene_order_within_module": str(gene_order),
                        "seed_status": "seed_for_feasibility",
                        "expansion_status": "not_yet_expanded_with_reactome_go_genage",
                    }
                )
    if not rows:
        raise ValueError(f"No maintenance genes parsed from {path}")
    return rows


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_report(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    module_counts: dict[str, int] = {}
    for row in rows:
        module_counts[row["maintenance_module"]] = module_counts.get(row["maintenance_module"], 0) + 1

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Maintenance Gene Sets Report",
                "",
                f"Seed genes: {len(rows)}",
                f"Maintenance modules: {len(module_counts)}",
                "",
                "## Module Counts",
                *[f"- {module}: {count}" for module, count in sorted(module_counts.items())],
                "",
                "## Interpretation",
                "This is a conservative human-centered seed gene set for the Week 3 feasibility pass. It is intentionally small and should be expanded later with Reactome, GO, GenAge/HAGR, and UniProt evidence before manuscript-level pathway claims.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    rows = parse_seed_yaml(args.input)
    write_tsv(args.output, rows)
    write_report(args.report, rows)
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
