"""Fit OLS models on the exact species subset recovered in the tree."""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
from collections import defaultdict

import numpy as np


DEFAULT_INPUT = pathlib.Path("data/processed/pgls_trait_table.tsv")
DEFAULT_OUTPUT = pathlib.Path("results/tables/tree_subset_ols_models.tsv")
DEFAULT_RESIDUALS = pathlib.Path("data/processed/tree_subset_ols_residuals.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/tree_subset_ols_report.md")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_float(value: str) -> float | None:
    try:
        if value.strip() == "":
            return None
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def design_matrix(rows: list[dict[str, str]], predictors: list[str]) -> tuple[np.ndarray, list[str]]:
    columns = [np.ones(len(rows))]
    names = ["intercept"]
    for predictor in predictors:
        if predictor == "log10_body_mass_g":
            columns.append(np.array([log10_body_mass(row) for row in rows], dtype=float))
            names.append(predictor)
        else:
            levels = sorted({row[predictor] for row in rows})
            baseline = levels[0]
            for level in levels[1:]:
                columns.append(np.array([1.0 if row[predictor] == level else 0.0 for row in rows]))
                names.append(f"{predictor}={level}")
    return np.column_stack(columns), names


def fit_ols(rows: list[dict[str, str]], predictors: list[str]) -> dict:
    y = np.array([log10_lifespan(row) for row in rows], dtype=float)
    x, names = design_matrix(rows, predictors)
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    fitted = x @ beta
    residuals = y - fitted
    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")
    n = len(rows)
    p = x.shape[1]
    adj_r2 = 1 - (1 - r2) * (n - 1) / (n - p) if n > p + 1 else float("nan")
    return {
        "predictors": predictors,
        "coef_names": names,
        "coef": beta,
        "fitted": fitted,
        "residuals": residuals,
        "n": n,
        "p": p,
        "r2": r2,
        "adj_r2": adj_r2,
        "rmse": math.sqrt(ss_res / max(1, n - p)),
    }


def log10_body_mass(row: dict[str, str]) -> float | None:
    existing = parse_float(row.get("log10_body_mass_g", ""))
    if existing is not None:
        return existing
    value = parse_float(row.get("body_mass_g", ""))
    return math.log10(value) if value and value > 0 else None


def log10_lifespan(row: dict[str, str]) -> float | None:
    existing = parse_float(row.get("log10_max_lifespan_years", ""))
    if existing is not None:
        return existing
    value = parse_float(row.get("max_lifespan_years", ""))
    return math.log10(value) if value and value > 0 else None


def write_models(models: list[tuple[str, dict]], output: pathlib.Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["model", "n", "p", "r2", "adj_r2", "rmse", "terms", "coefficients"]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for name, model in models:
            writer.writerow(
                {
                    "model": name,
                    "n": model["n"],
                    "p": model["p"],
                    "r2": f"{model['r2']:.6f}",
                    "adj_r2": f"{model['adj_r2']:.6f}",
                    "rmse": f"{model['rmse']:.6f}",
                    "terms": ";".join(model["coef_names"]),
                    "coefficients": ";".join(f"{value:.8g}" for value in model["coef"]),
                }
            )


def write_residuals(rows: list[dict[str, str]], model: dict, output: pathlib.Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) + [
        "tree_subset_predicted_log10_lifespan",
        "tree_subset_residual_log10",
        "tree_subset_residual_ratio",
    ]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row, fitted, residual in zip(rows, model["fitted"], model["residuals"]):
            enriched = dict(row)
            enriched["tree_subset_predicted_log10_lifespan"] = f"{fitted:.6f}"
            enriched["tree_subset_residual_log10"] = f"{residual:.6f}"
            enriched["tree_subset_residual_ratio"] = f"{10 ** residual:.6f}"
            writer.writerow(enriched)


def write_report(rows: list[dict[str, str]], models: list[tuple[str, dict]], report: pathlib.Path) -> None:
    clade_counts = defaultdict(int)
    flight_counts = defaultdict(int)
    for row in rows:
        clade_counts[row["clade"]] += 1
        flight_counts[row["flight_status"]] += 1

    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "\n".join(
            [
                "# Tree-Subset OLS Report",
                "",
                f"Species in OpenTree-induced subtree trait table: {len(rows)}",
                "",
                "## Clade Counts",
                *[f"- {key}: {value}" for key, value in sorted(clade_counts.items())],
                "",
                "## Flight Status Counts",
                *[f"- {key}: {value}" for key, value in sorted(flight_counts.items())],
                "",
                "## OLS Models",
                *[
                    f"- {name}: n={model['n']}, adj_r2={model['adj_r2']:.4f}, rmse={model['rmse']:.4f}"
                    for name, model in models
                ],
                "",
                "## Caveat",
                "These are OLS checks on the same species subset that will enter PGLS. They are not phylogenetically corrected.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--residuals-output", type=pathlib.Path, default=DEFAULT_RESIDUALS)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    rows = [
        row
        for row in read_tsv(args.input)
        if parse_float(row.get("log10_body_mass_g", "")) is not None
        or log10_body_mass(row) is not None
    ]
    rows = [row for row in rows if log10_lifespan(row) is not None]
    models = [
        ("model_a_mass", fit_ols(rows, ["log10_body_mass_g"])),
        ("model_b_mass_flight_clade", fit_ols(rows, ["log10_body_mass_g", "flight_status", "clade"])),
    ]
    write_models(models, args.output)
    write_residuals(rows, models[-1][1], args.residuals_output)
    write_report(rows, models, args.report)
    print(f"Wrote {args.output}, {args.residuals_output}, and {args.report}")


if __name__ == "__main__":
    main()
