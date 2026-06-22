"""Summarize maintenance PGLS sensitivity results into module priorities."""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
from collections import defaultdict


DEFAULT_INPUT = pathlib.Path("results/tables/maintenance_pgls_sensitivity.tsv")
DEFAULT_OUTPUT = pathlib.Path("results/tables/maintenance_pgls_sensitivity_summary.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/maintenance_pgls_sensitivity_summary_report.md")


CORE_MODEL = "mass_clade_module"
FOUNDATION_SUBSETS = {
    "all_primary",
    "coverage_ge_0.60",
    "exclude_human",
    "leave_out_Mammalia_Chiroptera",
    "leave_out_Mammalia_nonChiroptera",
    "leave_out_Reptilia",
}
STRESS_SUBSETS = {
    "coverage_ge_0.75",
    "exclude_top_abs_residual_5",
    "leave_out_Aves",
    "tier1_only",
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


def as_float(value: str) -> float:
    try:
        if value in {"", "NA", "NaN"}:
            return math.nan
        return float(value)
    except ValueError:
        return math.nan


def priority_label(foundational_hits: int, stress_hits: int, positive_fraction: float) -> str:
    if foundational_hits >= 5 and positive_fraction >= 0.8:
        if stress_hits >= 2:
            return "high_robust"
        return "high_but_stress_sensitive"
    if foundational_hits >= 3 and positive_fraction >= 0.7:
        return "moderate"
    return "exploratory"


def summarize(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["model"] != CORE_MODEL or row["error"]:
            continue
        groups[row["maintenance_module"]].append(row)

    summaries = []
    for module, module_rows in groups.items():
        tested = len(module_rows)
        positive = [row for row in module_rows if as_float(row["module_estimate"]) > 0]
        significant = [
            row
            for row in module_rows
            if as_float(row["module_p_bh_by_subset_model"]) < 0.05
        ]
        foundational = [
            row
            for row in significant
            if row["subset"] in FOUNDATION_SUBSETS
        ]
        stress = [
            row
            for row in significant
            if row["subset"] in STRESS_SUBSETS
        ]
        all_primary = next((row for row in module_rows if row["subset"] == "all_primary"), {})
        leave_out_aves = next((row for row in module_rows if row["subset"] == "leave_out_Aves"), {})
        coverage75 = next((row for row in module_rows if row["subset"] == "coverage_ge_0.75"), {})
        tier1 = next((row for row in module_rows if row["subset"] == "tier1_only"), {})
        pos_frac = len(positive) / tested if tested else 0.0
        summaries.append(
            {
                "maintenance_module": module,
                "model": CORE_MODEL,
                "subsets_tested": str(tested),
                "positive_estimate_subsets": str(len(positive)),
                "positive_fraction": f"{pos_frac:.3f}",
                "bh_significant_subsets": str(len(significant)),
                "foundation_significant_subsets": str(len(foundational)),
                "stress_significant_subsets": str(len(stress)),
                "all_primary_estimate": all_primary.get("module_estimate", ""),
                "all_primary_p": all_primary.get("module_p", ""),
                "all_primary_bh": all_primary.get("module_p_bh_by_subset_model", ""),
                "leave_out_aves_p": leave_out_aves.get("module_p", ""),
                "coverage_ge_0.75_p": coverage75.get("module_p", ""),
                "tier1_only_p": tier1.get("module_p", ""),
                "priority": priority_label(len(foundational), len(stress), pos_frac),
            }
        )
    summaries.sort(
        key=lambda row: (
            {
                "high_robust": 0,
                "high_but_stress_sensitive": 1,
                "moderate": 2,
                "exploratory": 3,
            }[row["priority"]],
            -int(row["foundation_significant_subsets"]),
            as_float(row["all_primary_p"]),
        )
    )
    return summaries


def write_report(path: pathlib.Path, summaries: list[dict[str, str]]) -> None:
    top = summaries[:7]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Maintenance PGLS Sensitivity Summary Report",
                "",
                f"Modules summarized: {len(summaries)}",
                "",
                "## Module Priorities",
                *[
                    f"- {row['maintenance_module']}: {row['priority']}, foundation_hits={row['foundation_significant_subsets']}, stress_hits={row['stress_significant_subsets']}, all_primary_p={row['all_primary_p']}"
                    for row in top
                ],
                "",
                "## Interpretation",
                "Foundation subsets test whether the signal survives common exclusions. Stress subsets intentionally remove informative coverage or clades and are expected to be underpowered. A high-but-stress-sensitive label means the module is promising but likely dependent on bird/Tier2 annotation rescue structure.",
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

    summaries = summarize(read_tsv(args.input))
    fields = list(summaries[0].keys()) if summaries else []
    write_tsv(args.output, summaries, fields)
    write_report(args.report, summaries)
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
