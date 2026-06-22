"""Build module-level genome-maintenance scores from validated candidate matrix."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import defaultdict


DEFAULT_INPUT = pathlib.Path("data/processed/ortholog_matrix_primary_diamond_validated_candidates.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/maintenance_scores_primary.tsv")
DEFAULT_WIDE_OUTPUT = pathlib.Path("data/processed/maintenance_scores_primary_wide.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/maintenance_scores_primary_report.md")

FOUND_STATUSES = {
    "ncbi_gene_candidate",
    "gff_rescue_candidate",
    "diamond_validated_protein_candidate",
}

CONFIDENCE_WEIGHTS = {
    "high": 1.0,
    "medium": 0.8,
    "low": 0.5,
    "": 0.0,
}


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def is_found(row: dict[str, str]) -> bool:
    return row["final_candidate_status"] in FOUND_STATUSES


def confidence(row: dict[str, str]) -> str:
    return row.get("final_candidate_confidence", "") or row.get("ortholog_confidence", "")


def confidence_weight(row: dict[str, str]) -> float:
    return CONFIDENCE_WEIGHTS.get(confidence(row), 0.0)


def source_group(row: dict[str, str]) -> str:
    status = row["final_candidate_status"]
    if status == "ncbi_gene_candidate":
        return "ncbi_gene"
    if status == "gff_rescue_candidate":
        return "gff_rescue"
    if status == "diamond_validated_protein_candidate":
        return "diamond_protein"
    return "not_found"


def build_scores(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row["scientific_name"], row["maintenance_module"])].append(row)

    scores = []
    species_meta = {
        row["scientific_name"]: {
            "clade": row["clade"],
            "flight_status": row["flight_status"],
            "genome_panel_version": row["genome_panel_version"],
            "genome_analysis_tier": row["genome_analysis_tier"],
        }
        for row in rows
    }
    for (species, module), group_rows in groups.items():
        found = [row for row in group_rows if is_found(row)]
        not_found = [row for row in group_rows if not is_found(row)]
        n_total = len(group_rows)
        weighted_sum = sum(confidence_weight(row) for row in found)
        counts = defaultdict(int)
        conf_counts = defaultdict(int)
        for row in group_rows:
            counts[source_group(row)] += 1
            if is_found(row):
                conf_counts[confidence(row) or "unlabeled"] += 1
        meta = species_meta[species]
        scores.append(
            {
                "scientific_name": species,
                "genome_panel_version": meta["genome_panel_version"],
                "clade": meta["clade"],
                "flight_status": meta["flight_status"],
                "genome_analysis_tier": meta["genome_analysis_tier"],
                "maintenance_module": module,
                "genes_total": str(n_total),
                "genes_found": str(len(found)),
                "genes_unresolved": str(len(not_found)),
                "coverage_fraction": f"{len(found) / n_total if n_total else 0.0:.6f}",
                "confidence_weighted_score": f"{weighted_sum / n_total if n_total else 0.0:.6f}",
                "ncbi_gene_count": str(counts["ncbi_gene"]),
                "gff_rescue_count": str(counts["gff_rescue"]),
                "diamond_protein_count": str(counts["diamond_protein"]),
                "not_found_count": str(counts["not_found"]),
                "high_confidence_count": str(conf_counts["high"]),
                "medium_confidence_count": str(conf_counts["medium"]),
                "low_confidence_count": str(conf_counts["low"]),
            }
        )
    scores.sort(key=lambda row: (row["scientific_name"], row["maintenance_module"]))
    return scores


def build_wide(scores: list[dict[str, str]]) -> list[dict[str, str]]:
    modules = sorted({row["maintenance_module"] for row in scores})
    by_species: dict[str, dict[str, str]] = {}
    for row in scores:
        species = row["scientific_name"]
        if species not in by_species:
            by_species[species] = {
                "scientific_name": species,
                "clade": row["clade"],
                "flight_status": row["flight_status"],
                "genome_analysis_tier": row["genome_analysis_tier"],
            }
        by_species[species][f"{row['maintenance_module']}_coverage"] = row["coverage_fraction"]
        by_species[species][f"{row['maintenance_module']}_score"] = row["confidence_weighted_score"]
    fields = ["scientific_name", "clade", "flight_status", "genome_analysis_tier"]
    for module in modules:
        fields.extend([f"{module}_coverage", f"{module}_score"])
    rows = list(by_species.values())
    rows.sort(key=lambda row: row["scientific_name"])
    for row in rows:
        for field in fields:
            row.setdefault(field, "")
    return rows, fields


def write_report(path: pathlib.Path, scores: list[dict[str, str]]) -> None:
    species = sorted({row["scientific_name"] for row in scores})
    modules = sorted({row["maintenance_module"] for row in scores})
    mean_by_module = []
    for module in modules:
        rows = [row for row in scores if row["maintenance_module"] == module]
        cov = sum(float(row["coverage_fraction"]) for row in rows) / len(rows)
        score = sum(float(row["confidence_weighted_score"]) for row in rows) / len(rows)
        mean_by_module.append((module, cov, score))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Maintenance Scores Primary Report",
                "",
                f"Species scored: {len(species)}",
                f"Modules scored: {len(modules)}",
                f"Long-format rows: {len(scores)}",
                "",
                "## Mean Module Scores",
                *[
                    f"- {module}: coverage={cov:.3f}, weighted_score={score:.3f}"
                    for module, cov, score in mean_by_module
                ],
                "",
                "## Interpretation",
                "Scores summarize candidate ortholog coverage, not expression or functional activity. Confidence-weighted scores give less weight to GFF rescue and lower-confidence protein rescue rows.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--wide-output", type=pathlib.Path, default=DEFAULT_WIDE_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    rows = read_tsv(args.input)
    scores = build_scores(rows)
    fields = list(scores[0].keys()) if scores else []
    write_tsv(args.output, scores, fields)
    wide_rows, wide_fields = build_wide(scores)
    write_tsv(args.wide_output, wide_rows, wide_fields)
    write_report(args.report, scores)
    print(f"Wrote {args.output}, {args.wide_output}, and {args.report}")


if __name__ == "__main__":
    main()
