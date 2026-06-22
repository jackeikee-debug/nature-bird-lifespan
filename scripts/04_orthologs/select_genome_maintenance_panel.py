"""Select primary and sensitivity genome-maintenance analysis panels."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import Counter


DEFAULT_INPUT = pathlib.Path("data/processed/genome_availability_audit.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/genome_maintenance_species.tsv")
DEFAULT_PRIMARY_OUTPUT = pathlib.Path("data/processed/genome_maintenance_species_primary.tsv")
DEFAULT_SENSITIVITY_OUTPUT = pathlib.Path("data/processed/genome_maintenance_species_sensitivity.tsv")
DEFAULT_VERSION_TABLE = pathlib.Path("data/processed/genome_maintenance_panel_versions.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/genome_maintenance_panel_report.md")

PRIMARY_TIERS = {
    "tier1_refseq_annotated_chromosome",
    "tier2_annotated",
}

SENSITIVITY_TIERS = {
    "tier1_refseq_annotated_chromosome",
    "tier2_annotated",
    "tier3_assembly_only",
}


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def panel_role(row: dict[str, str]) -> str:
    if row["clade"] == "Aves":
        return "bird_mechanism_panel"
    if row["clade"] == "Mammalia_Chiroptera":
        return "bat_convergence_panel"
    if row["clade"] == "Mammalia_nonChiroptera":
        return "mammal_control_panel"
    if row["clade"] == "Reptilia":
        return "reptile_outgroup_panel"
    return "other"


def risk_label(row: dict[str, str]) -> str:
    tier = row["genome_analysis_tier"]
    if tier == "tier1_refseq_annotated_chromosome":
        return "low"
    if tier == "tier2_annotated":
        return "moderate"
    if tier == "tier3_assembly_only":
        return "high"
    return "not_recommended"


def next_step(row: dict[str, str], version: str) -> str:
    tier = row["genome_analysis_tier"]
    if version == "primary":
        return "ortholog_mapping_and_pathway_scoring"
    if tier == "tier3_assembly_only":
        return "assembly_presence_sensitivity_only; exclude_from_primary_claims"
    return "ortholog_mapping_and_pathway_scoring"


def inclusion_reason(row: dict[str, str], version: str) -> str:
    tier = row["genome_analysis_tier"]
    if version == "primary":
        if tier == "tier1_refseq_annotated_chromosome":
            return "RefSeq annotated chromosome/complete assembly"
        return "annotated assembly suitable for first-pass ortholog mapping"
    if tier == "tier3_assembly_only":
        return "assembly exists but lacks NCBI annotation report; sensitivity only"
    if tier == "tier1_refseq_annotated_chromosome":
        return "retained from primary panel"
    return "retained from primary panel"


def build_panel(rows: list[dict[str, str]], keep_tiers: set[str], version: str) -> list[dict[str, str]]:
    panel = []
    for row in rows:
        if row["genome_analysis_tier"] not in keep_tiers:
            continue
        selected = dict(row)
        selected["genome_panel_version"] = version
        selected["mechanism_panel_role"] = panel_role(row)
        selected["ortholog_priority"] = {
            "tier1_refseq_annotated_chromosome": "high",
            "tier2_annotated": "medium",
            "tier3_assembly_only": "low",
        }.get(row["genome_analysis_tier"], "defer")
        selected["genome_quality_risk"] = risk_label(row)
        selected["inclusion_reason"] = inclusion_reason(row, version)
        selected["recommended_next_step"] = next_step(row, version)
        panel.append(selected)
    panel.sort(key=lambda item: (item["clade"], item["genome_analysis_tier"], item["scientific_name"]))
    return panel


def summarize_panel(panel: list[dict[str, str]]) -> tuple[Counter[str], Counter[str], Counter[str], Counter[str]]:
    return (
        Counter(row["genome_analysis_tier"] for row in panel),
        Counter(row["clade"] for row in panel),
        Counter(row["mechanism_panel_role"] for row in panel),
        Counter(row["genome_quality_risk"] for row in panel),
    )


def version_rows(
    rows: list[dict[str, str]],
    primary: list[dict[str, str]],
    sensitivity: list[dict[str, str]],
) -> list[dict[str, str]]:
    primary_names = {row["scientific_name"] for row in primary}
    sensitivity_names = {row["scientific_name"] for row in sensitivity}
    table = []
    for row in rows:
        in_primary = row["scientific_name"] in primary_names
        in_sensitivity = row["scientific_name"] in sensitivity_names
        if in_primary:
            assignment = "primary_and_sensitivity"
        elif in_sensitivity:
            assignment = "sensitivity_only"
        else:
            assignment = "excluded_or_deferred"
        table.append(
            {
                "scientific_name": row["scientific_name"],
                "clade": row["clade"],
                "flight_status": row["flight_status"],
                "genome_analysis_tier": row["genome_analysis_tier"],
                "best_assembly_accession": row["best_assembly_accession"],
                "has_annotation_report": row["has_annotation_report"],
                "genome_panel_assignment": assignment,
                "include_primary": "yes" if in_primary else "no",
                "include_sensitivity": "yes" if in_sensitivity else "no",
                "reason": (
                    "primary ortholog/pathway panel"
                    if in_primary
                    else "tier3 sensitivity panel only"
                    if in_sensitivity
                    else row["manual_review_reason"] or "no suitable assembly for mechanism panel"
                ),
            }
        )
    table.sort(key=lambda item: (item["genome_panel_assignment"], item["clade"], item["scientific_name"]))
    return table


def counter_lines(counter: Counter[str]) -> list[str]:
    if not counter:
        return ["- none: 0"]
    return [f"- {key}: {value}" for key, value in sorted(counter.items())]


def write_report(
    rows: list[dict[str, str]],
    primary: list[dict[str, str]],
    sensitivity: list[dict[str, str]],
    report: pathlib.Path,
) -> None:
    primary_tiers, primary_clades, primary_roles, primary_risks = summarize_panel(primary)
    sens_tiers, sens_clades, sens_roles, sens_risks = summarize_panel(sensitivity)
    sensitivity_only = [
        row for row in sensitivity if row["genome_analysis_tier"] == "tier3_assembly_only"
    ]
    excluded_clade_counts = Counter(
        row["clade"] for row in rows if row["genome_analysis_tier"] not in SENSITIVITY_TIERS
    )
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "\n".join(
            [
                "# Genome Maintenance Panel Report",
                "",
                f"Input species from genome audit: {len(rows)}",
                f"Primary Tier 1/2 species: {len(primary)}",
                f"Sensitivity Tier 1/2/3 species: {len(sensitivity)}",
                f"Sensitivity-only Tier 3 species: {len(sensitivity_only)}",
                "",
                "## Primary Panel Tiers",
                *counter_lines(primary_tiers),
                "",
                "## Primary Panel Clades",
                *counter_lines(primary_clades),
                "",
                "## Primary Panel Roles",
                *counter_lines(primary_roles),
                "",
                "## Primary Panel Genome Quality Risk",
                *counter_lines(primary_risks),
                "",
                "## Sensitivity Panel Tiers",
                *counter_lines(sens_tiers),
                "",
                "## Sensitivity Panel Clades",
                *counter_lines(sens_clades),
                "",
                "## Sensitivity Panel Roles",
                *counter_lines(sens_roles),
                "",
                "## Sensitivity Panel Genome Quality Risk",
                *counter_lines(sens_risks),
                "",
                "## Excluded or Deferred After Sensitivity Panel",
                *counter_lines(excluded_clade_counts),
                "",
                "## Interpretation",
                "The primary panel should define the first ortholog and pathway-score claims. It keeps annotated Tier 1/2 genomes only.",
                "",
                "The sensitivity panel adds Tier 3 assembly-only species to test whether broad phylogenetic or clade patterns are stable when coverage is expanded. Tier 3 species should not be used as the basis for primary gene-loss, copy-number, or annotation-dependent claims.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--primary-output", type=pathlib.Path, default=DEFAULT_PRIMARY_OUTPUT)
    parser.add_argument("--sensitivity-output", type=pathlib.Path, default=DEFAULT_SENSITIVITY_OUTPUT)
    parser.add_argument("--version-table", type=pathlib.Path, default=DEFAULT_VERSION_TABLE)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    rows = read_tsv(args.input)
    primary = build_panel(rows, PRIMARY_TIERS, "primary")
    sensitivity = build_panel(rows, SENSITIVITY_TIERS, "sensitivity")
    versions = version_rows(rows, primary, sensitivity)
    fields = list(primary[0].keys()) if primary else list(rows[0].keys())
    sensitivity_fields = list(sensitivity[0].keys()) if sensitivity else fields
    version_fields = list(versions[0].keys()) if versions else []
    write_tsv(args.primary_output, primary, fields)
    write_tsv(args.output, primary, fields)
    write_tsv(args.sensitivity_output, sensitivity, sensitivity_fields)
    write_tsv(args.version_table, versions, version_fields)
    write_report(rows, primary, sensitivity, args.report)
    print(
        "Wrote "
        f"{args.primary_output}, {args.sensitivity_output}, "
        f"{args.version_table}, {args.output}, and {args.report}"
    )


if __name__ == "__main__":
    main()
