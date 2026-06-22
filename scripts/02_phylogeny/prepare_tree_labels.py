"""Prepare taxonomy and tree labels for PGLS inputs.

This script does not download trees. It creates auditable label files that can
be submitted to OpenTree, BirdTree, VertLife, or reptile tree sources and then
joined back to the lifespan residual table.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
from collections import Counter


DEFAULT_RESIDUALS = pathlib.Path("data/processed/lifespan_residuals.tsv")
DEFAULT_TAXONOMY_AUDIT = pathlib.Path("data/processed/taxonomy_audit.tsv")
DEFAULT_PGLS_SPECIES = pathlib.Path("data/processed/pgls_species.tsv")
DEFAULT_TREE_AUDIT = pathlib.Path("data/processed/tree_label_audit.tsv")
DEFAULT_INPUT_DIR = pathlib.Path("data/processed/phylogeny_inputs")
DEFAULT_REPORT = pathlib.Path("results/reports/phylogeny_prep_report.md")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_lines(path: pathlib.Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def tree_label(name: str) -> str:
    cleaned = re.sub(r"\s+", " ", name.strip())
    cleaned = re.sub(r"[^\w\s.-]", "", cleaned)
    return cleaned.replace(" ", "_")


def name_rank(name: str) -> str:
    parts = name.split()
    if len(parts) == 2:
        return "binomial"
    if len(parts) == 3:
        return "trinomial_or_subspecies"
    return "nonstandard"


def source_hint(row: dict[str, str]) -> str:
    clade = row["clade"]
    if clade == "Aves":
        return "BirdTree_or_OpenTree"
    if clade in {"Mammalia_Chiroptera", "Mammalia_nonChiroptera"}:
        return "VertLife_or_OpenTree"
    if clade == "Reptilia":
        return "OpenTree_or_published_reptile_supertree"
    return "OpenTree"


def risk_flags(row: dict[str, str], tax: dict[str, str]) -> str:
    flags = []
    rank = name_rank(row["scientific_name"])
    if rank != "binomial":
        flags.append(rank)
    if row.get("scientific_name") != row.get("anage_matched_name"):
        flags.append("anage_alias_or_subspecies")
    if tax.get("data_quality") in {"low", "questionable"}:
        flags.append(f"anage_quality_{tax.get('data_quality')}")
    if tax.get("sample_size") in {"tiny", "small"}:
        flags.append(f"sample_size_{tax.get('sample_size')}")
    if row.get("sexual_maturity_years", "") == "":
        flags.append("missing_sexual_maturity")
    return ";".join(flags) if flags else "none"


def build_rows(
    residual_rows: list[dict[str, str]],
    taxonomy_rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    tax_by_name = {row["scientific_name"]: row for row in taxonomy_rows}
    pgls_rows = []
    audit_rows = []

    for row in residual_rows:
        tax = tax_by_name.get(row["scientific_name"], {})
        label_name = row.get("anage_matched_name") or row["scientific_name"]
        label = tree_label(label_name)
        pgls_row = {
            "scientific_name": row["scientific_name"],
            "tree_label": label,
            "tree_search_name": label_name,
            "name_rank": name_rank(row["scientific_name"]),
            "clade": row["clade"],
            "flight_status": row["flight_status"],
            "source_hint": source_hint(row),
            "anage_order": tax.get("anage_order", ""),
            "anage_family": tax.get("anage_family", ""),
            "common_name": row.get("common_name", ""),
            "body_mass_g": row.get("body_mass_g", ""),
            "max_lifespan_years": row.get("max_lifespan_years", ""),
            "sexual_maturity_years": row.get("sexual_maturity_years", ""),
            "lifespan_residual_log10": row.get("lifespan_residual_log10", ""),
            "lifespan_residual_ratio": row.get("lifespan_residual_ratio", ""),
            "data_quality": tax.get("data_quality", ""),
            "sample_size": tax.get("sample_size", ""),
            "specimen_origin": tax.get("specimen_origin", ""),
            "risk_flags": risk_flags(row, tax),
        }
        pgls_rows.append(pgls_row)
        audit_rows.append(
            {
                "scientific_name": row["scientific_name"],
                "tree_search_name": label_name,
                "tree_label": label,
                "name_rank": pgls_row["name_rank"],
                "clade": row["clade"],
                "source_hint": pgls_row["source_hint"],
                "anage_matched_name": row.get("anage_matched_name", ""),
                "anage_order": tax.get("anage_order", ""),
                "anage_family": tax.get("anage_family", ""),
                "taxonomy_issue": tax.get("taxonomy_issue", ""),
                "data_quality": tax.get("data_quality", ""),
                "sample_size": tax.get("sample_size", ""),
                "risk_flags": pgls_row["risk_flags"],
                "tree_match_status": "not_checked",
                "ott_id": "",
                "matched_tree_tip": "",
                "manual_tree_label": "",
                "notes": "",
            }
        )

    pgls_rows.sort(key=lambda item: (item["clade"], item["tree_label"]))
    audit_rows.sort(key=lambda item: (item["clade"], item["tree_label"]))
    return pgls_rows, audit_rows


def write_phylogeny_inputs(input_dir: pathlib.Path, rows: list[dict[str, str]]) -> None:
    all_labels = [row["tree_search_name"] for row in rows]
    bird_labels = [row["tree_search_name"] for row in rows if row["clade"] == "Aves"]
    mammal_labels = [
        row["tree_search_name"]
        for row in rows
        if row["clade"] in {"Mammalia_Chiroptera", "Mammalia_nonChiroptera"}
    ]
    bat_labels = [
        row["tree_search_name"]
        for row in rows
        if row["clade"] == "Mammalia_Chiroptera"
    ]
    reptile_labels = [row["tree_search_name"] for row in rows if row["clade"] == "Reptilia"]

    write_lines(input_dir / "opentree_tnrs_names.txt", all_labels)
    write_lines(input_dir / "birdtree_species.txt", bird_labels)
    write_lines(input_dir / "vertlife_mammal_species.txt", mammal_labels)
    write_lines(input_dir / "bat_species.txt", bat_labels)
    write_lines(input_dir / "reptile_species.txt", reptile_labels)


def write_report(report: pathlib.Path, rows: list[dict[str, str]], audit_rows: list[dict[str, str]]) -> None:
    clade_counts = Counter(row["clade"] for row in rows)
    source_counts = Counter(row["source_hint"] for row in rows)
    risk_counts = Counter()
    for row in rows:
        for flag in row["risk_flags"].split(";"):
            risk_counts[flag] += 1

    high_risk = [row for row in audit_rows if row["risk_flags"] != "none"]
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "\n".join(
            [
                "# Phylogeny Prep Report",
                "",
                f"PGLS-ready species: {len(rows)}",
                "",
                "## Clade Counts",
                *[f"- {key}: {value}" for key, value in sorted(clade_counts.items())],
                "",
                "## Suggested Tree Sources",
                *[f"- {key}: {value}" for key, value in sorted(source_counts.items())],
                "",
                "## Risk Flags",
                *[f"- {key}: {value}" for key, value in sorted(risk_counts.items())],
                "",
                "## Files",
                "- `data/processed/pgls_species.tsv`",
                "- `data/processed/tree_label_audit.tsv`",
                "- `data/processed/phylogeny_inputs/opentree_tnrs_names.txt`",
                "- `data/processed/phylogeny_inputs/birdtree_species.txt`",
                "- `data/processed/phylogeny_inputs/vertlife_mammal_species.txt`",
                "- `data/processed/phylogeny_inputs/reptile_species.txt`",
                "",
                "## Manual Review",
                f"Rows with one or more risk flags: {len(high_risk)}",
                "Review these before pruning or joining any external tree.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--residuals", type=pathlib.Path, default=DEFAULT_RESIDUALS)
    parser.add_argument("--taxonomy-audit", type=pathlib.Path, default=DEFAULT_TAXONOMY_AUDIT)
    parser.add_argument("--pgls-output", type=pathlib.Path, default=DEFAULT_PGLS_SPECIES)
    parser.add_argument("--tree-audit-output", type=pathlib.Path, default=DEFAULT_TREE_AUDIT)
    parser.add_argument("--input-dir", type=pathlib.Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    pgls_rows, audit_rows = build_rows(read_tsv(args.residuals), read_tsv(args.taxonomy_audit))
    pgls_fields = [
        "scientific_name",
        "tree_label",
        "tree_search_name",
        "name_rank",
        "clade",
        "flight_status",
        "source_hint",
        "anage_order",
        "anage_family",
        "common_name",
        "body_mass_g",
        "max_lifespan_years",
        "sexual_maturity_years",
        "lifespan_residual_log10",
        "lifespan_residual_ratio",
        "data_quality",
        "sample_size",
        "specimen_origin",
        "risk_flags",
    ]
    audit_fields = [
        "scientific_name",
        "tree_search_name",
        "tree_label",
        "name_rank",
        "clade",
        "source_hint",
        "anage_matched_name",
        "anage_order",
        "anage_family",
        "taxonomy_issue",
        "data_quality",
        "sample_size",
        "risk_flags",
        "tree_match_status",
        "ott_id",
        "matched_tree_tip",
        "manual_tree_label",
        "notes",
    ]
    write_tsv(args.pgls_output, pgls_rows, pgls_fields)
    write_tsv(args.tree_audit_output, audit_rows, audit_fields)
    write_phylogeny_inputs(args.input_dir, pgls_rows)
    write_report(args.report, pgls_rows, audit_rows)
    print(f"Wrote {args.pgls_output}, {args.tree_audit_output}, and {args.report}")


if __name__ == "__main__":
    main()

