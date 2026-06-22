"""Build Week 4 sequence-validated transposon score variants."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import defaultdict

import numpy as np
import pandas as pd


DEFAULT_MATRIX = pathlib.Path("data/processed/ortholog_matrix_primary_week4_sequence_validated.tsv")
DEFAULT_TRAITS = pathlib.Path("data/processed/pgls_trait_table.tsv")
DEFAULT_PGLS = pathlib.Path("data/processed/pgls_first_pass_residuals.tsv")
DEFAULT_LONG = pathlib.Path("data/processed/maintenance_scores_week4_sequence_validated.tsv")
DEFAULT_WIDE = pathlib.Path("data/processed/maintenance_scores_week4_sequence_validated_wide.tsv")
DEFAULT_MERGED = pathlib.Path("data/processed/maintenance_lifespan_week4_sequence_validated.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/maintenance_scores_week4_sequence_validated_report.md")

BASE_FOUND = {
    "ncbi_gene_candidate",
    "gff_rescue_candidate",
    "diamond_validated_protein_candidate",
}
STRICT_TRANSPOSON = {
    "ncbi_gene_candidate",
    "week4_sequence_supported_candidate",
}
WEAK_TRANSPOSON = {
    "ncbi_gene_candidate",
    "week4_sequence_supported_candidate",
    "week4_sequence_weak_candidate",
}
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


def base_confidence(row: dict[str, str]) -> str:
    return row.get("final_candidate_confidence", "") or row.get("ortholog_confidence", "")


def week4_confidence(row: dict[str, str]) -> str:
    return row.get("week4_candidate_confidence", "") or base_confidence(row)


def status_for_variant(row: dict[str, str], variant: str) -> str:
    if row["maintenance_module"] != "transposon_suppression":
        return row["final_candidate_status"]
    if variant in {"transposon_sequence_strict", "transposon_sequence_weak_inclusive"}:
        return row.get("week4_candidate_status", row["final_candidate_status"])
    return row["final_candidate_status"]


def found_for_variant(row: dict[str, str], variant: str) -> bool:
    status = status_for_variant(row, variant)
    if row["maintenance_module"] != "transposon_suppression":
        return status in BASE_FOUND
    if variant == "transposon_sequence_strict":
        return status in STRICT_TRANSPOSON
    if variant == "transposon_sequence_weak_inclusive":
        return status in WEAK_TRANSPOSON
    return status in BASE_FOUND


def confidence_for_variant(row: dict[str, str], variant: str) -> str:
    if row["maintenance_module"] == "transposon_suppression" and variant.startswith("transposon_sequence"):
        return week4_confidence(row)
    return base_confidence(row)


def build_scores(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row["scientific_name"], row["maintenance_module"])].append(row)

    variants = {
        "week3_all_validated": "Week 3 validated source logic.",
        "transposon_sequence_strict": "Week 3 logic for non-transposon modules; transposon uses NCBI direct plus reciprocal-supported candidates only.",
        "transposon_sequence_weak_inclusive": "Strict sequence-validated transposon score plus reciprocal-weak candidates.",
    }
    out = []
    for variant, description in variants.items():
        for (species, module), group_rows in groups.items():
            n_total = len(group_rows)
            found_rows = [row for row in group_rows if found_for_variant(row, variant)]
            weighted_sum = sum(
                CONFIDENCE_WEIGHTS.get(confidence_for_variant(row, variant), 0.0)
                for row in found_rows
            )
            meta = group_rows[0]
            status_counts = defaultdict(int)
            for row in found_rows:
                status_counts[status_for_variant(row, variant)] += 1
            out.append(
                {
                    "score_variant": variant,
                    "variant_description": description,
                    "scientific_name": species,
                    "genome_panel_version": meta["genome_panel_version"],
                    "clade": meta["clade"],
                    "flight_status": meta["flight_status"],
                    "genome_analysis_tier": meta["genome_analysis_tier"],
                    "maintenance_module": module,
                    "genes_total": str(n_total),
                    "genes_found": str(len(found_rows)),
                    "genes_unresolved": str(n_total - len(found_rows)),
                    "coverage_fraction": f"{len(found_rows) / n_total if n_total else 0.0:.6f}",
                    "confidence_weighted_score": f"{weighted_sum / n_total if n_total else 0.0:.6f}",
                    "ncbi_gene_count": str(status_counts["ncbi_gene_candidate"]),
                    "gff_rescue_count": str(status_counts["gff_rescue_candidate"]),
                    "diamond_protein_count": str(status_counts["diamond_validated_protein_candidate"]),
                    "week4_sequence_supported_count": str(status_counts["week4_sequence_supported_candidate"]),
                    "week4_sequence_weak_count": str(status_counts["week4_sequence_weak_candidate"]),
                }
            )
    out.sort(key=lambda row: (row["score_variant"], row["scientific_name"], row["maintenance_module"]))
    return out


def build_wide(scores: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[str]]:
    modules = sorted({row["maintenance_module"] for row in scores})
    by_key: dict[tuple[str, str], dict[str, str]] = {}
    for row in scores:
        key = (row["score_variant"], row["scientific_name"])
        if key not in by_key:
            by_key[key] = {
                "score_variant": row["score_variant"],
                "variant_description": row["variant_description"],
                "scientific_name": row["scientific_name"],
                "clade": row["clade"],
                "flight_status": row["flight_status"],
                "genome_analysis_tier": row["genome_analysis_tier"],
            }
        by_key[key][f"{row['maintenance_module']}_coverage"] = row["coverage_fraction"]
        by_key[key][f"{row['maintenance_module']}_score"] = row["confidence_weighted_score"]
    fields = [
        "score_variant",
        "variant_description",
        "scientific_name",
        "clade",
        "flight_status",
        "genome_analysis_tier",
    ]
    for module in modules:
        fields.extend([f"{module}_coverage", f"{module}_score"])
    wide = list(by_key.values())
    for row in wide:
        for field in fields:
            row.setdefault(field, "")
    wide.sort(key=lambda row: (row["score_variant"], row["scientific_name"]))
    return wide, fields


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
    trans = df[df["maintenance_module"] == "transposon_suppression"].copy()
    trans["score"] = pd.to_numeric(trans["confidence_weighted_score"], errors="coerce")
    trans["coverage"] = pd.to_numeric(trans["coverage_fraction"], errors="coerce")
    summary = trans.groupby("score_variant", as_index=False).agg(
        mean_coverage=("coverage", "mean"),
        mean_score=("score", "mean"),
        sequence_supported=("week4_sequence_supported_count", lambda x: pd.to_numeric(x, errors="coerce").sum()),
        sequence_weak=("week4_sequence_weak_count", lambda x: pd.to_numeric(x, errors="coerce").sum()),
    )
    lines = [
        "# Week 4 Sequence-Validated Maintenance Scores Report",
        "",
        f"Species per variant: {merged['scientific_name'].nunique()}",
        f"Merged rows: {len(merged)}",
        "",
        "## Transposon Summary",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- {row['score_variant']}: mean_coverage={row['mean_coverage']:.3f}, "
            f"mean_score={row['mean_score']:.3f}, sequence_supported={int(row['sequence_supported'])}, "
            f"sequence_weak={int(row['sequence_weak'])}"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "These scores keep the non-transposon modules comparable to Week 3 while testing stricter transposon-suppression scoring after reciprocal DIAMOND validation.",
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

    rows = read_tsv(args.matrix)
    scores = build_scores(rows)
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
