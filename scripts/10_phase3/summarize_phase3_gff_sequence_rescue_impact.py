"""Summarize model impact of Phase 3 GFF sequence rescue."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


VARIANT = "phase2_W3_full_background_sensitivity"
MODULE = "transposon_repeat_suppression"


def read(path: pathlib.Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def module_rank(pgls: pd.DataFrame, model: str, module: str) -> tuple[int | None, float | None, float | None]:
    sub = pgls[(pgls["score_variant"] == VARIANT) & (pgls["model"] == model)].copy()
    sub = numeric(sub, ["module_p", "module_estimate"])
    sub = sub[sub["error"] == ""].sort_values("module_p")
    sub["rank"] = range(1, len(sub) + 1)
    hit = sub[sub["maintenance_module"] == module]
    if hit.empty:
        return None, None, None
    row = hit.iloc[0]
    return int(row["rank"]), float(row["module_estimate"]), float(row["module_p"])


def interpretation_text(
    base_rank: int | None,
    base_p: float | None,
    new_rank: int | None,
    new_p: float | None,
    coverage_delta: float,
) -> str:
    if base_rank is None or new_rank is None or base_p is None or new_p is None:
        return (
            "The rescue changed coverage, but one of the compared residual models did not return a transposon/repeat rank. "
            "Treat this as a workflow audit result rather than biological evidence."
        )
    if new_rank < base_rank:
        rank_phrase = "improves in rank"
    elif new_rank > base_rank:
        rank_phrase = "drops in rank"
    else:
        rank_phrase = "keeps the same rank"
    if new_p < base_p:
        p_phrase = "with a smaller nominal P value"
    elif new_p > base_p:
        p_phrase = "with a larger nominal P value"
    else:
        p_phrase = "with an unchanged nominal P value"
    if abs(coverage_delta) < 0.005:
        coverage_phrase = "Coverage changes are minimal"
    else:
        coverage_phrase = "Coverage improves in the intended low-coverage species"
    return (
        f"{coverage_phrase}. In the all-module residual model the transposon/repeat module remains positive and {rank_phrase}, "
        f"{p_phrase}. This supports a cautious Phase 3 conclusion: orthology rescue is technically feasible, but the primary "
        "manuscript claim should remain validation-focused until broader rescue and external database confirmation are complete."
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-scores", type=pathlib.Path, required=True)
    parser.add_argument("--rescued-scores", type=pathlib.Path, required=True)
    parser.add_argument("--baseline-pgls", type=pathlib.Path, required=True)
    parser.add_argument("--rescued-pgls", type=pathlib.Path, required=True)
    parser.add_argument("--baseline-highcov", type=pathlib.Path, required=True)
    parser.add_argument("--rescued-highcov", type=pathlib.Path, required=True)
    parser.add_argument("--overlay-summary", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    base_scores = numeric(read(args.baseline_scores), ["coverage_fraction", "confidence_weighted_score"])
    rescued_scores = numeric(read(args.rescued_scores), ["coverage_fraction", "confidence_weighted_score"])
    base_trans = base_scores[(base_scores["score_variant"] == VARIANT) & (base_scores["maintenance_module_v2"] == MODULE)]
    rescued_trans = rescued_scores[(rescued_scores["score_variant"] == VARIANT) & (rescued_scores["maintenance_module_v2"] == MODULE)]
    merged_scores = base_trans[
        ["scientific_name", "coverage_fraction", "confidence_weighted_score"]
    ].merge(
        rescued_trans[["scientific_name", "coverage_fraction", "confidence_weighted_score"]],
        on="scientific_name",
        suffixes=("_baseline", "_rescued"),
    )
    merged_scores["coverage_delta"] = merged_scores["coverage_fraction_rescued"] - merged_scores["coverage_fraction_baseline"]
    merged_scores["score_delta"] = merged_scores["confidence_weighted_score_rescued"] - merged_scores["confidence_weighted_score_baseline"]

    baseline_pgls = read(args.baseline_pgls)
    rescued_pgls = read(args.rescued_pgls)
    rows = []
    for label, pgls in [("baseline", baseline_pgls), ("rescued", rescued_pgls)]:
        for model in ["residual_module", "mass_clade_module"]:
            rank, estimate, pval = module_rank(pgls, model, MODULE)
            rows.append(
                {
                    "comparison_layer": "all_module_pgls",
                    "dataset": label,
                    "model": model,
                    "metric": "transposon_rank_estimate_p",
                    "rank": rank,
                    "estimate": estimate,
                    "p": pval,
                    "mean_transposon_coverage": float(base_trans["coverage_fraction"].mean() if label == "baseline" else rescued_trans["coverage_fraction"].mean()),
                    "notes": "",
                }
            )

    base_high = numeric(read(args.baseline_highcov), ["estimate", "p", "n", "birds"])
    rescued_high = numeric(read(args.rescued_highcov), ["estimate", "p", "n", "birds"])
    for label, high in [("baseline", base_high), ("rescued", rescued_high)]:
        sub = high[(high["model"] == "pgls_clade_residual_transposon") & (high["subset"].isin(["all_primary", "transposon_coverage_ge_0_50", "birds_only_transposon_coverage_ge_0_50"]))]
        for _, row in sub.iterrows():
            rows.append(
                {
                    "comparison_layer": "high_coverage_subset",
                    "dataset": label,
                    "model": row["subset"],
                    "metric": "estimate_p_n_birds",
                    "rank": "",
                    "estimate": row["estimate"],
                    "p": row["p"],
                    "mean_transposon_coverage": "",
                    "notes": f"n={int(row['n'])};birds={int(row['birds'])}",
                }
            )

    output = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)

    overlay = read(args.overlay_summary)
    updated_rows = int((overlay["update_status"] == "updated").sum()) if not overlay.empty else 0
    low_species = merged_scores.sort_values("coverage_delta", ascending=False).head(8)
    base_res_rank, base_res_est, base_res_p = module_rank(baseline_pgls, "residual_module", MODULE)
    new_res_rank, new_res_est, new_res_p = module_rank(rescued_pgls, "residual_module", MODULE)
    coverage_delta = float(rescued_trans["coverage_fraction"].mean() - base_trans["coverage_fraction"].mean())
    lines = [
        "# Phase 3 Sequence Rescue Impact Report",
        "",
        f"Strict sequence rows applied to matrix: {updated_rows}",
        f"Mean transposon coverage baseline: {base_trans['coverage_fraction'].mean():.3f}",
        f"Mean transposon coverage rescued: {rescued_trans['coverage_fraction'].mean():.3f}",
        f"Mean transposon coverage delta: {coverage_delta:.3f}",
        "",
        "## Main Residual Model Impact",
        "",
        f"- Baseline transposon rank: {base_res_rank}, estimate={base_res_est:.4g}, p={base_res_p:.4g}.",
        f"- Rescued transposon rank: {new_res_rank}, estimate={new_res_est:.4g}, p={new_res_p:.4g}.",
        "",
        "## Largest Species Coverage Gains",
    ]
    for _, row in low_species.iterrows():
        lines.append(
            f"- {row['scientific_name']}: coverage {row['coverage_fraction_baseline']:.3f} -> "
            f"{row['coverage_fraction_rescued']:.3f} (delta {row['coverage_delta']:.3f})"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            interpretation_text(base_res_rank, base_res_p, new_res_rank, new_res_p, coverage_delta),
            "",
            "## Outputs",
            f"- impact table: `{args.output}`",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
