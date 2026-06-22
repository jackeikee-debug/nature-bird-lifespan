"""Fit a first-pass lifespan residual model.

This intentionally starts with an ordinary least-squares model on log body mass.
PGLS and phylogenetic mixed models should replace or complement this once the
species tree is added.
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
import statistics
from collections import Counter


DEFAULT_INPUT = pathlib.Path("data/processed/species_master.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/lifespan_residuals.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week1_feasibility_report.md")


def parse_float(value: str) -> float | None:
    try:
        if value.strip() == "":
            return None
        parsed = float(value)
    except ValueError:
        return None
    if parsed <= 0:
        return None
    return parsed


def read_rows(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fit_simple_linear(x: list[float], y: list[float]) -> tuple[float, float]:
    x_mean = statistics.fmean(x)
    y_mean = statistics.fmean(y)
    numerator = sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(x, y))
    denominator = sum((xi - x_mean) ** 2 for xi in x)
    if denominator == 0:
        raise RuntimeError("Cannot fit model: all body masses are identical.")
    slope = numerator / denominator
    intercept = y_mean - slope * x_mean
    return intercept, slope


def write_outputs(
    rows: list[dict[str, str]],
    output: pathlib.Path,
    report: pathlib.Path,
) -> None:
    model_rows = []
    for row in rows:
        mass = parse_float(row.get("body_mass_g", ""))
        lifespan = parse_float(row.get("max_lifespan_years", ""))
        if mass is None or lifespan is None:
            continue
        model_rows.append((row, math.log10(mass), math.log10(lifespan)))

    if len(model_rows) < 8:
        raise RuntimeError(
            f"Need at least 8 species with mass and lifespan; found {len(model_rows)}."
        )

    intercept, slope = fit_simple_linear(
        [item[1] for item in model_rows],
        [item[2] for item in model_rows],
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys()) + [
        "log10_body_mass_g",
        "log10_max_lifespan_years",
        "predicted_log10_lifespan",
        "lifespan_residual_log10",
        "lifespan_residual_ratio",
    ]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row, log_mass, log_life in model_rows:
            predicted = intercept + slope * log_mass
            residual = log_life - predicted
            enriched = dict(row)
            enriched.update(
                {
                    "log10_body_mass_g": f"{log_mass:.6f}",
                    "log10_max_lifespan_years": f"{log_life:.6f}",
                    "predicted_log10_lifespan": f"{predicted:.6f}",
                    "lifespan_residual_log10": f"{residual:.6f}",
                    "lifespan_residual_ratio": f"{10 ** residual:.6f}",
                }
            )
            writer.writerow(enriched)

    report.parent.mkdir(parents=True, exist_ok=True)
    by_flight: dict[str, list[float]] = {}
    by_clade: dict[str, list[float]] = {}
    for row, log_mass, log_life in model_rows:
        predicted = intercept + slope * log_mass
        by_flight.setdefault(row["flight_status"], []).append(log_life - predicted)
        by_clade.setdefault(row["clade"], []).append(log_life - predicted)

    flight_lines = [
        f"- {group}: n={len(values)}, mean residual={statistics.fmean(values):.4f}"
        for group, values in sorted(by_flight.items())
    ]
    clade_lines = [
        f"- {group}: n={len(values)}, mean residual={statistics.fmean(values):.4f}"
        for group, values in sorted(by_clade.items())
    ]
    missing_anage = [
        row["scientific_name"] for row in rows if row.get("anage_record_present") == "no"
    ]
    missing_model = [
        row["scientific_name"]
        for row in rows
        if not parse_float(row.get("body_mass_g", ""))
        or not parse_float(row.get("max_lifespan_years", ""))
    ]
    inclusion_counts = Counter(row.get("inclusion_role", "") for row in rows)
    inclusion_lines = [
        f"- {role}: n={count}" for role, count in sorted(inclusion_counts.items())
    ]
    report.write_text(
        "\n".join(
            [
                "# Week 1 Feasibility Report",
                "",
                f"Total species in master table: {len(rows)}",
                f"Species with model-ready AnAge data: {len(model_rows)}",
                f"Model: log10(max_lifespan_years) = {intercept:.4f} + {slope:.4f} * log10(body_mass_g)",
                "",
                "## Flight Status Residual Summary",
                *flight_lines,
                "",
                "## Clade Residual Summary",
                *clade_lines,
                "",
                "## Inclusion Role Counts",
                *inclusion_lines,
                "",
                "## Missing or Incomplete Records",
                f"- Missing AnAge match: {', '.join(missing_anage) if missing_anage else 'none'}",
                f"- Missing model fields: {', '.join(missing_model) if missing_model else 'none'}",
                "",
                "## Next Checks",
                "- Manually audit the expanded 240-species seed list.",
                "- Curate habitat, diet, genome availability, and annotation quality for auto-selected species.",
                "- Add taxonomy harmonization and phylogenetic tree labels.",
                "- Replace OLS with PGLS once the tree is available.",
                "- Audit body-mass and maximum-longevity outliers manually.",
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

    rows = read_rows(args.input)
    write_outputs(rows, args.output, args.report)
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
