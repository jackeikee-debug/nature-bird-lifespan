"""Build transposon gene leave-one-out scores from the full sequence matrix."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import defaultdict

import numpy as np
import pandas as pd


DEFAULT_MATRIX = pathlib.Path("data/processed/ortholog_matrix_primary_week4_full_sequence_validated.tsv")
DEFAULT_TRAITS = pathlib.Path("data/processed/pgls_trait_table.tsv")
DEFAULT_PGLS = pathlib.Path("data/processed/pgls_first_pass_residuals.tsv")
DEFAULT_LONG = pathlib.Path("data/processed/week4_transposon_gene_leave_one_out_scores.tsv")
DEFAULT_WIDE = pathlib.Path("data/processed/week4_transposon_gene_leave_one_out_scores_wide.tsv")
DEFAULT_MERGED = pathlib.Path("data/processed/week4_transposon_gene_leave_one_out_lifespan.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week4_transposon_gene_leave_one_out_scores_report.md")

TRANSPOSON_GENES = ["MOV10", "PIWIL1", "PIWIL2", "SETDB1", "TRIM28"]
STRICT_FOUND = {"week4_sequence_supported_candidate"}
WEAK_FOUND = {"week4_sequence_supported_candidate", "week4_sequence_weak_candidate"}
CONFIDENCE_WEIGHTS = {"high": 1.0, "medium": 0.8, "low": 0.5, "": 0.0}


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def found(row: dict[str, str], weak_inclusive: bool) -> bool:
    status = row.get("week4_candidate_status", "")
    return status in (WEAK_FOUND if weak_inclusive else STRICT_FOUND)


def confidence(row: dict[str, str]) -> str:
    return row.get("week4_candidate_confidence", "") or row.get("final_candidate_confidence", "")


def build_scores(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    trans = [row for row in rows if row["maintenance_module"] == "transposon_suppression"]
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in trans:
        groups[row["scientific_name"]].append(row)

    variants = [("all_5_genes_strict", "", False), ("all_5_genes_weak_inclusive", "", True)]
    for gene in TRANSPOSON_GENES:
        variants.append((f"drop_{gene}_strict", gene, False))
        variants.append((f"drop_{gene}_weak_inclusive", gene, True))

    score_rows = []
    for species, group_rows in groups.items():
        meta = group_rows[0]
        for variant, dropped_gene, weak_inclusive in variants:
            selected = [
                row for row in group_rows
                if not dropped_gene or row["human_gene_symbol"] != dropped_gene
            ]
            found_rows = [row for row in selected if found(row, weak_inclusive)]
            weighted_sum = sum(CONFIDENCE_WEIGHTS.get(confidence(row), 0.0) for row in found_rows)
            score_rows.append(
                {
                    "score_variant": variant,
                    "dropped_gene": dropped_gene,
                    "weak_inclusive": "yes" if weak_inclusive else "no",
                    "scientific_name": species,
                    "clade": meta["clade"],
                    "flight_status": meta["flight_status"],
                    "genome_analysis_tier": meta["genome_analysis_tier"],
                    "genes_total": str(len(selected)),
                    "genes_found": str(len(found_rows)),
                    "coverage_fraction": f"{len(found_rows) / len(selected) if selected else 0.0:.6f}",
                    "transposon_suppression_score": f"{weighted_sum / len(selected) if selected else 0.0:.6f}",
                }
            )
    score_rows.sort(key=lambda row: (row["score_variant"], row["scientific_name"]))
    return score_rows


def build_wide(scores: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[str]]:
    rows = []
    fields = [
        "score_variant",
        "dropped_gene",
        "weak_inclusive",
        "scientific_name",
        "clade",
        "flight_status",
        "genome_analysis_tier",
        "transposon_suppression_coverage",
        "transposon_suppression_score",
    ]
    for row in scores:
        rows.append(
            {
                "score_variant": row["score_variant"],
                "dropped_gene": row["dropped_gene"],
                "weak_inclusive": row["weak_inclusive"],
                "scientific_name": row["scientific_name"],
                "clade": row["clade"],
                "flight_status": row["flight_status"],
                "genome_analysis_tier": row["genome_analysis_tier"],
                "transposon_suppression_coverage": row["coverage_fraction"],
                "transposon_suppression_score": row["transposon_suppression_score"],
            }
        )
    return rows, fields


def build_merged(wide_path: pathlib.Path, traits_path: pathlib.Path, pgls_path: pathlib.Path) -> pd.DataFrame:
    wide = pd.read_csv(wide_path, sep="\t")
    traits = pd.read_csv(traits_path, sep="\t")
    pgls = pd.read_csv(pgls_path, sep="\t")
    keep_traits = [
        "scientific_name",
        "opentree_tip_label",
        "body_mass_g",
        "max_lifespan_years",
        "lifespan_residual_log10",
        "lifespan_residual_ratio",
        "data_quality",
        "sample_size",
        "specimen_origin",
    ]
    merged = wide.merge(traits[keep_traits], on="scientific_name", how="left")
    merged["log10_body_mass_g"] = np.log10(pd.to_numeric(merged["body_mass_g"], errors="coerce"))
    merged = merged.merge(pgls, on="opentree_tip_label", how="left")
    return merged


def write_report(path: pathlib.Path, scores: list[dict[str, str]], merged: pd.DataFrame) -> None:
    df = pd.DataFrame(scores)
    df["score"] = pd.to_numeric(df["transposon_suppression_score"], errors="coerce")
    summary = df.groupby("score_variant", as_index=False).agg(
        mean_score=("score", "mean"),
        min_score=("score", "min"),
        max_score=("score", "max"),
    )
    lines = [
        "# Week 4 Transposon Gene Leave-One-Out Scores Report",
        "",
        f"Species per variant: {merged['scientific_name'].nunique()}",
        f"Score variants: {df['score_variant'].nunique()}",
        f"Merged rows: {len(merged)}",
        "",
        "## Mean Scores",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- {row['score_variant']}: mean={row['mean_score']:.3f}, "
            f"min={row['min_score']:.3f}, max={row['max_score']:.3f}"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "These scores test whether the transposon signal is driven by one seed gene. Strict variants count only reciprocal-supported candidates; weak-inclusive variants also count reciprocal-weak candidates.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, default=DEFAULT_MATRIX)
    parser.add_argument("--traits", type=pathlib.Path, default=DEFAULT_TRAITS)
    parser.add_argument("--pgls", type=pathlib.Path, default=DEFAULT_PGLS)
    parser.add_argument("--long-output", type=pathlib.Path, default=DEFAULT_LONG)
    parser.add_argument("--wide-output", type=pathlib.Path, default=DEFAULT_WIDE)
    parser.add_argument("--merged-output", type=pathlib.Path, default=DEFAULT_MERGED)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    matrix = read_tsv(args.matrix)
    scores = build_scores(matrix)
    write_tsv(args.long_output, scores, list(scores[0].keys()))
    wide, fields = build_wide(scores)
    write_tsv(args.wide_output, wide, fields)
    merged = build_merged(args.wide_output, args.traits, args.pgls)
    args.merged_output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.merged_output, sep="\t", index=False)
    write_report(args.report, scores, merged)
    print(f"Wrote {args.long_output}, {args.wide_output}, {args.merged_output}, and {args.report}")


if __name__ == "__main__":
    main()
