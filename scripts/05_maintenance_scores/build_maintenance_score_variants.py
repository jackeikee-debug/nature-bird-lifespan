"""Build Week 4 maintenance-score variants for rescue/quality sensitivity tests."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import defaultdict

import numpy as np
import pandas as pd


DEFAULT_MATRIX = pathlib.Path("data/processed/ortholog_matrix_primary_diamond_validated_candidates.tsv")
DEFAULT_TRAITS = pathlib.Path("data/processed/pgls_trait_table.tsv")
DEFAULT_PGLS = pathlib.Path("data/processed/pgls_first_pass_residuals.tsv")
DEFAULT_LONG = pathlib.Path("data/processed/maintenance_scores_week4_variants.tsv")
DEFAULT_WIDE = pathlib.Path("data/processed/maintenance_scores_week4_variants_wide.tsv")
DEFAULT_MERGED = pathlib.Path("data/processed/maintenance_lifespan_week4_variants.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/maintenance_score_variants_week4_report.md")

CONFIDENCE_WEIGHTS = {
    "high": 1.0,
    "medium": 0.8,
    "low": 0.5,
    "": 0.0,
}

VARIANTS = {
    "all_validated": {
        "description": "NCBI Gene, GFF rescue, and DIAMOND-validated protein candidates.",
        "statuses": {
            "ncbi_gene_candidate",
            "gff_rescue_candidate",
            "diamond_validated_protein_candidate",
        },
        "min_confidence": "low",
    },
    "high_confidence_only": {
        "description": "Only high-confidence candidates from any accepted source.",
        "statuses": {
            "ncbi_gene_candidate",
            "gff_rescue_candidate",
            "diamond_validated_protein_candidate",
        },
        "min_confidence": "high",
    },
    "ncbi_only": {
        "description": "Only direct NCBI Gene symbol/taxid candidates; excludes GFF and protein rescue.",
        "statuses": {"ncbi_gene_candidate"},
        "min_confidence": "low",
    },
    "no_protein_rescue": {
        "description": "NCBI Gene plus GFF rescue candidates; excludes DIAMOND protein-only rescue.",
        "statuses": {"ncbi_gene_candidate", "gff_rescue_candidate"},
        "min_confidence": "low",
    },
}

CONFIDENCE_RANK = {"": 0, "low": 1, "medium": 2, "high": 3}


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def confidence(row: dict[str, str]) -> str:
    return row.get("final_candidate_confidence", "") or row.get("ortholog_confidence", "")


def passes_variant(row: dict[str, str], variant: dict[str, object]) -> bool:
    status = row["final_candidate_status"]
    if status not in variant["statuses"]:
        return False
    min_conf = str(variant["min_confidence"])
    return CONFIDENCE_RANK.get(confidence(row), 0) >= CONFIDENCE_RANK[min_conf]


def confidence_weight(row: dict[str, str], found: bool) -> float:
    if not found:
        return 0.0
    return CONFIDENCE_WEIGHTS.get(confidence(row), 0.0)


def build_scores(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row["scientific_name"], row["maintenance_module"])].append(row)

    out = []
    for variant_name, variant in VARIANTS.items():
        for (species, module), group_rows in groups.items():
            found_flags = [passes_variant(row, variant) for row in group_rows]
            n_total = len(group_rows)
            genes_found = sum(found_flags)
            weighted_sum = sum(
                confidence_weight(row, found)
                for row, found in zip(group_rows, found_flags, strict=True)
            )
            meta = group_rows[0]
            status_counts = defaultdict(int)
            confidence_counts = defaultdict(int)
            for row, found in zip(group_rows, found_flags, strict=True):
                if found:
                    status_counts[row["final_candidate_status"]] += 1
                    confidence_counts[confidence(row) or "unlabeled"] += 1
            out.append(
                {
                    "score_variant": variant_name,
                    "variant_description": str(variant["description"]),
                    "scientific_name": species,
                    "genome_panel_version": meta["genome_panel_version"],
                    "clade": meta["clade"],
                    "flight_status": meta["flight_status"],
                    "genome_analysis_tier": meta["genome_analysis_tier"],
                    "maintenance_module": module,
                    "genes_total": str(n_total),
                    "genes_found": str(genes_found),
                    "genes_unresolved": str(n_total - genes_found),
                    "coverage_fraction": f"{genes_found / n_total if n_total else 0.0:.6f}",
                    "confidence_weighted_score": f"{weighted_sum / n_total if n_total else 0.0:.6f}",
                    "ncbi_gene_count": str(status_counts["ncbi_gene_candidate"]),
                    "gff_rescue_count": str(status_counts["gff_rescue_candidate"]),
                    "diamond_protein_count": str(status_counts["diamond_validated_protein_candidate"]),
                    "high_confidence_count": str(confidence_counts["high"]),
                    "medium_confidence_count": str(confidence_counts["medium"]),
                    "low_confidence_count": str(confidence_counts["low"]),
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
    lines = [
        "# Week 4 Maintenance Score Variants Report",
        "",
        f"Score variants: {len(VARIANTS)}",
        f"Species per variant: {merged['scientific_name'].nunique()}",
        f"Merged rows: {len(merged)}",
        "",
        "## Variant Definitions",
        "",
    ]
    for name, spec in VARIANTS.items():
        lines.append(f"- `{name}`: {spec['description']}")
    lines.extend(["", "## Mean Module Scores by Variant", ""])
    summary = (
        df.assign(score=pd.to_numeric(df["confidence_weighted_score"], errors="coerce"))
        .groupby(["score_variant", "maintenance_module"], as_index=False)["score"]
        .mean()
        .sort_values(["score_variant", "maintenance_module"])
    )
    for _, row in summary.iterrows():
        lines.append(f"- {row['score_variant']} / {row['maintenance_module']}: mean_score={row['score']:.3f}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "These variants test whether Week 3 module-lifespan signals are driven by annotation rescue or lower-confidence rows. `all_validated` should reproduce the Week 3 primary score logic. `ncbi_only` is the strictest source-control variant, while `high_confidence_only` is the strictest confidence-control variant.",
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
    wide_rows, wide_fields = build_wide(scores)
    write_tsv(args.wide_output, wide_rows, wide_fields)
    merged = build_merged(args.wide_output, args.traits, args.pgls)
    args.merged_output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.merged_output, sep="\t", index=False)
    write_report(args.report, scores, merged)
    print(f"Wrote {args.long_output}, {args.wide_output}, {args.merged_output}, and {args.report}")


if __name__ == "__main__":
    main()
