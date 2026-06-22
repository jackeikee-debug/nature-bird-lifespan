"""Build interim P2.4 module scores from priority-1 expanded matrix."""

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


def found(row: pd.Series) -> bool:
    for col in ["week4_candidate_status", "final_candidate_status", "combined_candidate_status"]:
        value = str(row.get(col, ""))
        if value:
            return value in FOUND_STATUSES
    return False


def confidence(row: pd.Series) -> str:
    for col in ["week4_candidate_confidence", "final_candidate_confidence", "ortholog_confidence"]:
        value = row.get(col, "")
        if isinstance(value, str) and value:
            return value
    return ""


def variant_panels(eligibility: pd.DataFrame) -> dict[str, tuple[str, set[str]]]:
    baseline = set(eligibility.loc[eligibility["base_v2_scoring_group"] == "strict_ready", "human_gene_symbol"])
    strict = set(eligibility.loc[eligibility["strict_score_allowed"].astype(bool), "human_gene_symbol"])
    priority1_sensitivity = set(
        eligibility.loc[
            eligibility["v2_scoring_group"].isin(
                {
                    "strict_ready",
                    "strict_sequence_supported",
                    "domain_supported_paralog_guard",
                    "domain_supported_manual_upgrade_candidate",
                }
            ),
            "human_gene_symbol",
        ]
    )
    w2_crossdb_sensitivity = set(
        eligibility.loc[
            eligibility["v2_scoring_group"].isin(
                {
                    "strict_ready",
                    "strict_sequence_supported",
                    "domain_supported_paralog_guard",
                    "domain_supported_manual_upgrade_candidate",
                    "crossdb_confirm",
                }
            ),
            "human_gene_symbol",
        ]
    )
    return {
        "phase2_seed41_baseline": (
            "Original species-complete 41-gene strict-ready panel.",
            baseline,
        ),
        "phase2_strict48_sequence_updated": (
            "41 strict-ready genes plus 7 DIAMOND/BLAST sequence-supported strict upgrades.",
            strict,
        ),
        "phase2_priority1_domain_sensitivity": (
            "Strict48 plus 9 InterProScan domain-supported priority-1 sensitivity genes.",
            priority1_sensitivity,
        ),
        "phase2_W2_crossdb_sensitivity": (
            "Priority1-domain sensitivity panel plus W2 crossdb_confirm full-species NCBI Gene expansion.",
            w2_crossdb_sensitivity,
        ),
        "phase2_W3_DNA_repair_background_sensitivity": (
            "W2 crossdb sensitivity panel plus W3 standard-mapping DNA repair genes for background-control expansion.",
            w2_crossdb_sensitivity
            | set(
                eligibility.loc[
                    eligibility["maintenance_module_v2"].eq("DNA_repair_replication_stress")
                    & eligibility["v2_scoring_group"].eq("standard_mapping_pending"),
                    "human_gene_symbol",
                ]
            ),
        ),
        "phase2_W3_DNA_repair_proteostasis_background_sensitivity": (
            "W2 crossdb sensitivity panel plus W3 standard-mapping DNA repair and proteostasis genes for stronger background-control expansion.",
            w2_crossdb_sensitivity
            | set(
                eligibility.loc[
                    eligibility["maintenance_module_v2"].isin(
                        {"DNA_repair_replication_stress", "proteostasis_autophagy_mitophagy"}
                    )
                    & eligibility["v2_scoring_group"].eq("standard_mapping_pending"),
                    "human_gene_symbol",
                ]
            ),
        ),
        "phase2_W3_DNA_repair_proteostasis_cancer_background_sensitivity": (
            "W2 crossdb sensitivity panel plus W3 standard-mapping DNA repair, proteostasis, and cancer/senescence genes.",
            w2_crossdb_sensitivity
            | set(
                eligibility.loc[
                    eligibility["maintenance_module_v2"].isin(
                        {
                            "DNA_repair_replication_stress",
                            "proteostasis_autophagy_mitophagy",
                            "cancer_surveillance_senescence",
                        }
                    )
                    & eligibility["v2_scoring_group"].eq("standard_mapping_pending"),
                    "human_gene_symbol",
                ]
            ),
        ),
        "phase2_W3_full_background_sensitivity": (
            "W2 crossdb sensitivity panel plus all W3 standard-mapping background modules.",
            w2_crossdb_sensitivity
            | set(
                eligibility.loc[
                    eligibility["maintenance_module_v2"].isin(
                        {
                            "DNA_repair_replication_stress",
                            "proteostasis_autophagy_mitophagy",
                            "cancer_surveillance_senescence",
                            "inflammation_innate_immune_restraint",
                        }
                    )
                    & eligibility["v2_scoring_group"].eq("standard_mapping_pending"),
                    "human_gene_symbol",
                ]
            ),
        ),
    }


def build_scores(matrix: pd.DataFrame, eligibility: pd.DataFrame, species: pd.DataFrame) -> pd.DataFrame:
    gene_meta = eligibility.set_index("human_gene_symbol")
    variants = variant_panels(eligibility)
    rows = []
    for variant, (description, genes) in variants.items():
        panel = eligibility[eligibility["human_gene_symbol"].isin(genes)].copy()
        for _, sp in species.iterrows():
            sp_matrix = matrix[matrix["scientific_name"] == sp["scientific_name"]]
            for module, module_panel in panel.groupby("maintenance_module_v2"):
                module_genes = set(module_panel["human_gene_symbol"])
                module_rows = sp_matrix[sp_matrix["human_gene_symbol"].isin(module_genes)]
                by_gene = {row["human_gene_symbol"]: row for _, row in module_rows.iterrows()}
                genes_found = 0
                weighted = 0.0
                status_counts: defaultdict[str, int] = defaultdict(int)
                missing_matrix_rows = 0
                for gene in sorted(module_genes):
                    row = by_gene.get(gene)
                    if row is None:
                        missing_matrix_rows += 1
                        status_counts["missing_matrix_row"] += 1
                        continue
                    status = str(row.get("week4_candidate_status", "")) or str(row.get("final_candidate_status", ""))
                    status_counts[status] += 1
                    if found(row):
                        genes_found += 1
                        weighted += CONFIDENCE_WEIGHTS.get(confidence(row), 0.0)
                total = len(module_genes)
                rows.append(
                    {
                        "score_variant": variant,
                        "variant_description": description,
                        "scientific_name": sp["scientific_name"],
                        "clade": sp["clade"],
                        "flight_status": sp["flight_status"],
                        "genome_analysis_tier": sp["genome_analysis_tier"],
                        "maintenance_module_v2": module,
                        "genes_total": total,
                        "genes_found": genes_found,
                        "genes_unresolved": total - genes_found,
                        "missing_matrix_rows": missing_matrix_rows,
                        "coverage_fraction": genes_found / total if total else 0.0,
                        "confidence_weighted_score": weighted / total if total else 0.0,
                        "ncbi_gene_count": status_counts["ncbi_gene_candidate"],
                        "gff_rescue_count": status_counts["gff_rescue_candidate"],
                        "diamond_protein_count": status_counts["diamond_validated_protein_candidate"],
                        "week4_sequence_supported_count": status_counts["week4_sequence_supported_candidate"],
                        "priority1_no_ncbi_gene_count": status_counts[
                            "priority1_expansion_no_ncbi_gene_candidate"
                        ],
                    }
                )
    return pd.DataFrame(rows).sort_values(["score_variant", "scientific_name", "maintenance_module_v2"])


def build_wide(scores: pd.DataFrame) -> pd.DataFrame:
    modules = sorted(scores["maintenance_module_v2"].unique())
    rows = []
    for (variant, species_name), sub in scores.groupby(["score_variant", "scientific_name"]):
        first = sub.iloc[0]
        rec = {
            "score_variant": variant,
            "variant_description": first["variant_description"],
            "scientific_name": species_name,
            "clade": first["clade"],
            "flight_status": first["flight_status"],
            "genome_analysis_tier": first["genome_analysis_tier"],
        }
        for _, row in sub.iterrows():
            module = row["maintenance_module_v2"]
            rec[f"{module}_coverage"] = row["coverage_fraction"]
            rec[f"{module}_score"] = row["confidence_weighted_score"]
            rec[f"{module}_missing_matrix_rows"] = row["missing_matrix_rows"]
        for module in modules:
            rec.setdefault(f"{module}_coverage", "")
            rec.setdefault(f"{module}_score", "")
            rec.setdefault(f"{module}_missing_matrix_rows", "")
        rows.append(rec)
    return pd.DataFrame(rows).sort_values(["score_variant", "scientific_name"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, required=True)
    parser.add_argument("--eligibility", type=pathlib.Path, required=True)
    parser.add_argument("--species", type=pathlib.Path, required=True)
    parser.add_argument("--long-output", type=pathlib.Path, required=True)
    parser.add_argument("--wide-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t")
    eligibility = pd.read_csv(args.eligibility, sep="\t")
    species = pd.read_csv(args.species, sep="\t")
    scores = build_scores(matrix, eligibility, species)
    wide = build_wide(scores)

    args.long_output.parent.mkdir(parents=True, exist_ok=True)
    args.wide_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    scores.to_csv(args.long_output, sep="\t", index=False)
    wide.to_csv(args.wide_output, sep="\t", index=False)

    trans = scores[scores["maintenance_module_v2"] == "transposon_repeat_suppression"].copy()
    summary = trans.groupby("score_variant", as_index=False).agg(
        genes_total=("genes_total", "first"),
        mean_coverage=("coverage_fraction", "mean"),
        mean_score=("confidence_weighted_score", "mean"),
        total_missing_matrix_rows=("missing_matrix_rows", "sum"),
    )
    lines = [
        "# Phase 2 Priority-1 Expanded Scores Report",
        "",
        f"Species scored: {species['scientific_name'].nunique()}",
        f"Score variants: {scores['score_variant'].nunique()}",
        f"Long rows: {len(scores)}",
        "",
        "## Transposon Module Summary",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- {row['score_variant']}: genes_total={int(row['genes_total'])}, "
            f"mean_coverage={row['mean_coverage']:.3f}, mean_score={row['mean_score']:.3f}, "
            f"missing_matrix_rows={int(row['total_missing_matrix_rows'])}"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "These are interim P2.4 scores. The strict48 variant now has species-level rows for the 7 sequence-supported upgrades. The priority1-domain sensitivity variant includes domain-supported genes but remains a sensitivity analysis, not a main absence-claim panel.",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.long_output}")


if __name__ == "__main__":
    main()
