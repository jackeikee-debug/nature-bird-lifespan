"""Build W2 repeat/chromatin submodule scores."""

from __future__ import annotations

import argparse
import pathlib
from collections import defaultdict

import numpy as np
import pandas as pd


FOUND_STATUSES = {
    "ncbi_gene_candidate",
    "gff_rescue_candidate",
    "diamond_validated_protein_candidate",
    "week4_sequence_supported_candidate",
}
CONFIDENCE_WEIGHTS = {"high": 1.0, "medium": 0.8, "low": 0.5, "": 0.0, np.nan: 0.0}
W2_GROUPS = {
    "strict_ready",
    "strict_sequence_supported",
    "domain_supported_paralog_guard",
    "domain_supported_manual_upgrade_candidate",
    "crossdb_confirm",
}
TARGET_MODULES = {
    "transposon_repeat_suppression",
    "chromatin_repression_heterochromatin",
}


def found(row: pd.Series) -> bool:
    for col in ["week4_candidate_status", "final_candidate_status", "combined_candidate_status"]:
        value = str(row.get(col, ""))
        if value and value != "nan":
            return value in FOUND_STATUSES
    return False


def confidence(row: pd.Series) -> str:
    for col in ["week4_candidate_confidence", "final_candidate_confidence", "ortholog_confidence"]:
        value = row.get(col, "")
        if isinstance(value, str) and value and value != "nan":
            return value
    return ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, required=True)
    parser.add_argument("--eligibility", type=pathlib.Path, required=True)
    parser.add_argument("--species", type=pathlib.Path, required=True)
    parser.add_argument("--lifespan", type=pathlib.Path, required=True)
    parser.add_argument("--long-output", type=pathlib.Path, required=True)
    parser.add_argument("--wide-output", type=pathlib.Path, required=True)
    parser.add_argument("--merged-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t")
    eligibility = pd.read_csv(args.eligibility, sep="\t")
    species = pd.read_csv(args.species, sep="\t")
    lifespan = pd.read_csv(args.lifespan, sep="\t")
    lifespan = lifespan[lifespan["score_variant"] == "phase2_W2_crossdb_sensitivity"].copy()

    panel = eligibility[
        eligibility["v2_scoring_group"].isin(W2_GROUPS)
        & eligibility["maintenance_module_v2"].isin(TARGET_MODULES)
    ].copy()
    gene_meta = panel.set_index("human_gene_symbol")
    rows = []
    for _, sp in species.iterrows():
        sp_matrix = matrix[matrix["scientific_name"] == sp["scientific_name"]]
        by_gene = {row["human_gene_symbol"]: row for _, row in sp_matrix.iterrows()}
        for (module, submodule), sub_panel in panel.groupby(["maintenance_module_v2", "submodule_v2"]):
            genes = sorted(sub_panel["human_gene_symbol"].unique())
            genes_found = 0
            weighted = 0.0
            status_counts: defaultdict[str, int] = defaultdict(int)
            for gene in genes:
                row = by_gene.get(gene)
                if row is None:
                    status_counts["missing_matrix_row"] += 1
                    continue
                status = str(row.get("week4_candidate_status", "")) or str(row.get("final_candidate_status", ""))
                status_counts[status] += 1
                if found(row):
                    genes_found += 1
                    weighted += CONFIDENCE_WEIGHTS.get(confidence(row), 0.0)
            total = len(genes)
            rows.append(
                {
                    "score_variant": "phase2_W2_crossdb_sensitivity",
                    "scientific_name": sp["scientific_name"],
                    "clade": sp["clade"],
                    "flight_status": sp["flight_status"],
                    "genome_analysis_tier": sp["genome_analysis_tier"],
                    "maintenance_module_v2": module,
                    "submodule_v2": submodule,
                    "genes_total": total,
                    "genes_found": genes_found,
                    "genes_unresolved": total - genes_found,
                    "coverage_fraction": genes_found / total if total else 0.0,
                    "confidence_weighted_score": weighted / total if total else 0.0,
                    "ncbi_gene_count": status_counts["ncbi_gene_candidate"],
                    "week4_sequence_supported_count": status_counts["week4_sequence_supported_candidate"],
                    "missing_matrix_rows": status_counts["missing_matrix_row"],
                }
            )
    long = pd.DataFrame(rows).sort_values(["scientific_name", "maintenance_module_v2", "submodule_v2"])

    wide_rows = []
    for species_name, sub in long.groupby("scientific_name"):
        first = sub.iloc[0]
        rec = {
            "score_variant": "phase2_W2_crossdb_sensitivity",
            "scientific_name": species_name,
            "clade": first["clade"],
            "flight_status": first["flight_status"],
            "genome_analysis_tier": first["genome_analysis_tier"],
        }
        for _, row in sub.iterrows():
            key = row["submodule_v2"]
            rec[f"{key}_score"] = row["confidence_weighted_score"]
            rec[f"{key}_coverage"] = row["coverage_fraction"]
            rec[f"{key}_genes_total"] = row["genes_total"]
        wide_rows.append(rec)
    wide = pd.DataFrame(wide_rows).sort_values("scientific_name")
    merged = wide.merge(
        lifespan.drop(columns=["score_variant", "clade", "flight_status", "genome_analysis_tier"], errors="ignore"),
        on="scientific_name",
        how="left",
    )

    args.long_output.parent.mkdir(parents=True, exist_ok=True)
    args.wide_output.parent.mkdir(parents=True, exist_ok=True)
    args.merged_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(args.long_output, sep="\t", index=False)
    wide.to_csv(args.wide_output, sep="\t", index=False)
    merged.to_csv(args.merged_output, sep="\t", index=False)

    summary = (
        long.groupby(["maintenance_module_v2", "submodule_v2"], as_index=False)
        .agg(
            genes_total=("genes_total", "first"),
            mean_coverage=("coverage_fraction", "mean"),
            mean_score=("confidence_weighted_score", "mean"),
        )
        .sort_values(["maintenance_module_v2", "submodule_v2"])
    )
    lines = [
        "# Phase 2 W2 Submodule Scores",
        "",
        "## Submodule Coverage",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- {row['maintenance_module_v2']} / {row['submodule_v2']}: "
            f"genes={int(row['genes_total'])}, mean_coverage={row['mean_coverage']:.3f}, "
            f"mean_score={row['mean_score']:.3f}"
        )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.long_output}, {args.wide_output}, and {args.merged_output}")


if __name__ == "__main__":
    main()
