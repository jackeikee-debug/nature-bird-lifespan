"""Matched random gene-set tests for Phase 2 W2 repeat/chromatin scores."""

from __future__ import annotations

import argparse
import pathlib
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats


FOUND_STATUSES = {
    "ncbi_gene_candidate",
    "gff_rescue_candidate",
    "diamond_validated_protein_candidate",
    "week4_sequence_supported_candidate",
}
CONFIDENCE_WEIGHTS = {"high": 1.0, "medium": 0.8, "low": 0.5, "": 0.0}
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


def gene_weight(row: pd.Series) -> float:
    if not found(row):
        return 0.0
    return CONFIDENCE_WEIGHTS.get(confidence(row), 0.0)


def fit_xy(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    slope, intercept, r_value, p_value, stderr = stats.linregress(x, y)
    return {
        "n": float(len(y)),
        "estimate": float(slope),
        "p": float(p_value),
        "r": float(r_value),
        "r2": float(r_value**2),
    }


def empirical_p(null_values: np.ndarray, observed: float, side: str = "greater") -> float:
    if side == "greater":
        return float((np.sum(null_values >= observed) + 1) / (len(null_values) + 1))
    return float((np.sum(np.abs(null_values) >= abs(observed)) + 1) / (len(null_values) + 1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, required=True)
    parser.add_argument("--eligibility", type=pathlib.Path, required=True)
    parser.add_argument("--lifespan", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--null-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--n-permutations", type=int, default=5000)
    parser.add_argument("--coverage-window", type=float, default=0.08)
    parser.add_argument("--nearest-matching", action="store_true")
    parser.add_argument("--max-attempt-multiplier", type=int, default=30)
    parser.add_argument("--score-variant", default="phase2_W2_crossdb_sensitivity")
    parser.add_argument("--include-w3-dna-repair", action="store_true")
    parser.add_argument("--include-w3-proteostasis", action="store_true")
    parser.add_argument("--include-w3-cancer", action="store_true")
    parser.add_argument("--include-w3-inflammation", action="store_true")
    parser.add_argument("--seed", type=int, default=20260612)
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t")
    eligibility = pd.read_csv(args.eligibility, sep="\t")
    traits = pd.read_csv(args.lifespan, sep="\t")
    traits = traits[traits["score_variant"] == args.score_variant].copy()

    panel_mask = eligibility["v2_scoring_group"].isin(W2_GROUPS)
    if args.include_w3_dna_repair:
        panel_mask = panel_mask | (
            eligibility["maintenance_module_v2"].eq("DNA_repair_replication_stress")
            & eligibility["v2_scoring_group"].eq("standard_mapping_pending")
        )
    if args.include_w3_proteostasis:
        panel_mask = panel_mask | (
            eligibility["maintenance_module_v2"].eq("proteostasis_autophagy_mitophagy")
            & eligibility["v2_scoring_group"].eq("standard_mapping_pending")
        )
    if args.include_w3_cancer:
        panel_mask = panel_mask | (
            eligibility["maintenance_module_v2"].eq("cancer_surveillance_senescence")
            & eligibility["v2_scoring_group"].eq("standard_mapping_pending")
        )
    if args.include_w3_inflammation:
        panel_mask = panel_mask | (
            eligibility["maintenance_module_v2"].eq("inflammation_innate_immune_restraint")
            & eligibility["v2_scoring_group"].eq("standard_mapping_pending")
        )
    panel = eligibility[panel_mask].copy()
    genes = sorted(panel["human_gene_symbol"].unique())
    species = sorted(traits["scientific_name"].unique())

    matrix = matrix[matrix["human_gene_symbol"].isin(genes)].copy()
    by_pair = {
        (row["scientific_name"], row["human_gene_symbol"]): row
        for _, row in matrix.iterrows()
    }
    gene_scores: dict[str, np.ndarray] = {}
    gene_coverage = {}
    for gene in genes:
        vals = []
        for sp in species:
            row = by_pair.get((sp, gene))
            vals.append(gene_weight(row) if row is not None else 0.0)
        arr = np.array(vals, dtype=float)
        gene_scores[gene] = arr
        gene_coverage[gene] = float((arr > 0).mean())

    y_by_species = traits.set_index("scientific_name").loc[species, "pgls_model_c_mass_clade_residual"]
    y = pd.to_numeric(y_by_species, errors="coerce").to_numpy(dtype=float)
    valid_y = np.isfinite(y)

    gene_to_module = panel.set_index("human_gene_symbol")["maintenance_module_v2"].to_dict()
    rng = np.random.default_rng(args.seed)
    rows = []
    null_rows = []

    for module in sorted(TARGET_MODULES):
        target_genes = sorted(panel.loc[panel["maintenance_module_v2"] == module, "human_gene_symbol"].unique())
        target_size = len(target_genes)
        target_score = np.vstack([gene_scores[g] for g in target_genes]).mean(axis=0)
        target_coverage_mean = float(np.mean([gene_coverage[g] for g in target_genes]))
        obs = fit_xy(target_score[valid_y], y[valid_y])

        universe = [g for g in genes if gene_to_module.get(g) not in TARGET_MODULES]
        gene_cov = np.array([gene_coverage[g] for g in universe], dtype=float)
        attempts = 0
        candidates = []
        max_attempts = args.n_permutations * args.max_attempt_multiplier
        while attempts < max_attempts:
            attempts += 1
            sample = list(rng.choice(universe, size=target_size, replace=False))
            mean_cov = float(np.mean([gene_coverage[g] for g in sample]))
            if not args.nearest_matching and abs(mean_cov - target_coverage_mean) > args.coverage_window:
                continue
            score = np.vstack([gene_scores[g] for g in sample]).mean(axis=0)
            fit = fit_xy(score[valid_y], y[valid_y])
            candidates.append(
                {
                    "target_module": module,
                    "sample_size": target_size,
                    "mean_gene_coverage": mean_cov,
                    "coverage_distance": abs(mean_cov - target_coverage_mean),
                    "estimate": fit["estimate"],
                    "p": fit["p"],
                    "r": fit["r"],
                    "r2": fit["r2"],
                    "sampled_genes": ",".join(sample),
                }
            )

        if args.nearest_matching:
            candidates = sorted(candidates, key=lambda row: row["coverage_distance"])[: args.n_permutations]
        else:
            candidates = candidates[: args.n_permutations]
        accepted = len(candidates)
        if accepted == 0:
            raise RuntimeError(
                f"No matched random sets accepted for {module}. "
                "Rerun with --nearest-matching or a wider --coverage-window."
            )
        null_estimates = [row["estimate"] for row in candidates]
        null_p = [row["p"] for row in candidates]
        null_r = [row["r"] for row in candidates]
        null_r2 = [row["r2"] for row in candidates]
        null_cov = [row["mean_gene_coverage"] for row in candidates]
        for idx, row in enumerate(candidates, start=1):
            row = dict(row)
            row["permutation_id"] = idx
            null_rows.append(row)

        null_estimates_arr = np.array(null_estimates)
        null_r_arr = np.array(null_r)
        null_r2_arr = np.array(null_r2)
        rows.append(
            {
                "target_module": module,
                "target_gene_count": target_size,
                "target_mean_gene_coverage": target_coverage_mean,
                "accepted_permutations": accepted,
                "attempts": attempts,
                "coverage_window": args.coverage_window,
                "nearest_matching": args.nearest_matching,
                "matched_mean_coverage": float(np.mean(null_cov)),
                "matched_coverage_distance_mean": float(np.mean([abs(v - target_coverage_mean) for v in null_cov])),
                "observed_estimate": obs["estimate"],
                "observed_p": obs["p"],
                "observed_r": obs["r"],
                "observed_r2": obs["r2"],
                "null_estimate_mean": float(np.mean(null_estimates_arr)),
                "null_estimate_sd": float(np.std(null_estimates_arr, ddof=1)),
                "null_estimate_p95": float(np.quantile(null_estimates_arr, 0.95)),
                "null_r_mean": float(np.mean(null_r_arr)),
                "null_r_sd": float(np.std(null_r_arr, ddof=1)),
                "null_r_p95": float(np.quantile(null_r_arr, 0.95)),
                "null_r2_mean": float(np.mean(null_r2_arr)),
                "null_r2_p95": float(np.quantile(null_r2_arr, 0.95)),
                "empirical_p_estimate_greater": empirical_p(null_estimates_arr, obs["estimate"]),
                "empirical_p_r_greater": empirical_p(null_r_arr, obs["r"]),
                "empirical_p_abs_estimate": empirical_p(null_estimates_arr, obs["estimate"], side="abs"),
                "empirical_p_r2_greater": empirical_p(null_r2_arr, obs["r2"]),
            }
        )

    summary = pd.DataFrame(rows)
    null_df = pd.DataFrame(null_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.null_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output, sep="\t", index=False)
    null_df.to_csv(args.null_output, sep="\t", index=False)

    lines = [
        "# Phase 2 W2 Matched Random Gene-Set Tests",
        "",
        f"Score variant: {args.score_variant}",
        f"Background universe: W2-scored genes excluding repeat/chromatin target modules.",
        f"Include W3 DNA repair genes: {args.include_w3_dna_repair}",
        f"Include W3 proteostasis genes: {args.include_w3_proteostasis}",
        f"Include W3 cancer/senescence genes: {args.include_w3_cancer}",
        f"Include W3 inflammation genes: {args.include_w3_inflammation}",
        f"Permutations requested per target: {args.n_permutations}",
        f"Coverage matching window: +/- {args.coverage_window}",
        f"Nearest matching mode: {args.nearest_matching}",
        "",
        "## Results",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- {row['target_module']}: observed_estimate={row['observed_estimate']:.4f}, "
            f"observed_r={row['observed_r']:.4f}, observed_r2={row['observed_r2']:.4f}, "
            f"empirical_p_r={row['empirical_p_r_greater']:.4g}, "
            f"empirical_p_r2={row['empirical_p_r2_greater']:.4g}, "
            f"target_cov={row['target_mean_gene_coverage']:.3f}, matched_cov={row['matched_mean_coverage']:.3f}, "
            f"accepted={int(row['accepted_permutations'])}"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This is a coverage-matched random-set stress test using the currently W2-scored gene universe. It tests specificity against similarly observable maintenance genes, not against all possible genes in the genome.",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}, {args.null_output}, and {args.report}")


if __name__ == "__main__":
    main()
