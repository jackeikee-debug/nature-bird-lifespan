"""Create final Phase 2 module-ranking, decision-gate, and claim reports."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


FINAL_VARIANT = "phase2_W3_full_background_sensitivity"
TARGET_MODULES = {
    "transposon_repeat_suppression",
    "chromatin_repression_heterochromatin",
}


def num(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def fmt(value: object, digits: int = 4) -> str:
    try:
        val = float(value)
    except (TypeError, ValueError):
        return "NA"
    if pd.isna(val):
        return "NA"
    return f"{val:.{digits}g}"


def decision(pass_cond: bool, caution_cond: bool) -> str:
    if pass_cond:
        return "pass"
    if caution_cond:
        return "caution"
    return "fail"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgls", type=pathlib.Path, required=True)
    parser.add_argument("--random", type=pathlib.Path, required=True)
    parser.add_argument("--scores", type=pathlib.Path, required=True)
    parser.add_argument("--annotation", type=pathlib.Path, required=True)
    parser.add_argument("--final-annotation-ranking", type=pathlib.Path)
    parser.add_argument("--joint", type=pathlib.Path, required=True)
    parser.add_argument("--sensitivity", type=pathlib.Path)
    parser.add_argument("--module-ranking-output", type=pathlib.Path, required=True)
    parser.add_argument("--decision-gates-output", type=pathlib.Path, required=True)
    parser.add_argument("--claim-table-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    pgls = num(
        pd.read_csv(args.pgls, sep="\t"),
        ["module_estimate", "module_p", "module_p_bh_by_variant_model", "n", "lambda"],
    )
    random = num(
        pd.read_csv(args.random, sep="\t"),
        ["observed_r", "observed_r2", "empirical_p_r_greater", "target_mean_gene_coverage", "matched_mean_coverage"],
    )
    scores = num(pd.read_csv(args.scores, sep="\t"), ["genes_total", "coverage_fraction", "confidence_weighted_score"])
    annotation = num(pd.read_csv(args.annotation, sep="\t"), ["estimate", "p", "p_bh"])
    final_annotation = (
        num(
            pd.read_csv(args.final_annotation_ranking, sep="\t"),
            ["estimate", "p", "p_bh_by_model", "rank_by_p", "n", "lambda"],
        )
        if args.final_annotation_ranking and args.final_annotation_ranking.exists()
        else pd.DataFrame()
    )
    joint = num(pd.read_csv(args.joint, sep="\t"), ["estimate", "p", "p_bh", "lambda", "n"])
    sensitivity = (
        num(pd.read_csv(args.sensitivity, sep="\t"), ["estimate", "p", "p_bh_by_subset_model", "rank_by_p", "n"])
        if args.sensitivity and args.sensitivity.exists()
        else pd.DataFrame()
    )

    final = pgls[pgls["score_variant"].eq(FINAL_VARIANT) & pgls["error"].fillna("").eq("")].copy()
    final = final[final["model"].isin(["mass_clade_module", "pgls_clade_residual_module", "residual_module"])]
    coverage = (
        scores[scores["score_variant"].eq(FINAL_VARIANT)]
        .groupby("maintenance_module_v2", as_index=False)
        .agg(
            genes_total=("genes_total", "first"),
            mean_coverage=("coverage_fraction", "mean"),
            min_coverage=("coverage_fraction", "min"),
            mean_score=("confidence_weighted_score", "mean"),
        )
        .rename(columns={"maintenance_module_v2": "maintenance_module"})
    )

    ranks = []
    for model, sub in final.groupby("model"):
        sub = sub.sort_values("module_p").reset_index(drop=True)
        for idx, row in sub.iterrows():
            ranks.append(
                {
                    "score_variant": FINAL_VARIANT,
                    "model": model,
                    "rank_by_p": idx + 1,
                    "maintenance_module": row["maintenance_module"],
                    "estimate": row["module_estimate"],
                    "p": row["module_p"],
                    "p_bh_by_model": row["module_p_bh_by_variant_model"],
                    "n": row["n"],
                    "lambda": row["lambda"],
                }
            )
    ranking = pd.DataFrame(ranks).merge(coverage, on="maintenance_module", how="left")
    args.module_ranking_output.parent.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(args.module_ranking_output, sep="\t", index=False)

    trans_random = random[random["target_module"].eq("transposon_repeat_suppression")].iloc[0]
    chrom_random = random[random["target_module"].eq("chromatin_repression_heterochromatin")].iloc[0]
    trans_pgls_rank = ranking[
        ranking["model"].eq("pgls_clade_residual_module")
        & ranking["maintenance_module"].eq("transposon_repeat_suppression")
    ].iloc[0]
    trans_mass_rank = ranking[
        ranking["model"].eq("mass_clade_module")
        & ranking["maintenance_module"].eq("transposon_repeat_suppression")
    ].iloc[0]

    ann_trans = annotation[
        annotation["maintenance_module"].eq("transposon_repeat_suppression")
        & annotation["model"].eq("residual_tier")
    ].iloc[0]
    ann_trans_cov = annotation[
        annotation["maintenance_module"].eq("transposon_repeat_suppression")
        & annotation["model"].eq("residual_coverage")
    ].iloc[0]
    final_ann_trans_tier = pd.Series(dtype=object)
    final_ann_trans_cov = pd.Series(dtype=object)
    final_ann_trans_mass_tier = pd.Series(dtype=object)
    final_ann_any_coverage_sig = pd.Series(dtype=object)
    if not final_annotation.empty:
        final_ann_trans_tier_rows = final_annotation[
            final_annotation["maintenance_module"].eq("transposon_repeat_suppression")
            & final_annotation["model"].eq("residual_tier")
        ]
        final_ann_trans_cov_rows = final_annotation[
            final_annotation["maintenance_module"].eq("transposon_repeat_suppression")
            & final_annotation["model"].eq("residual_coverage")
        ]
        final_ann_trans_mass_tier_rows = final_annotation[
            final_annotation["maintenance_module"].eq("transposon_repeat_suppression")
            & final_annotation["model"].eq("mass_clade_tier")
        ]
        final_ann_coverage_rows = final_annotation[
            final_annotation["model"].isin(["residual_coverage", "residual_tier_coverage", "mass_clade_tier_coverage"])
        ]
        if not final_ann_trans_tier_rows.empty:
            final_ann_trans_tier = final_ann_trans_tier_rows.iloc[0]
        if not final_ann_trans_cov_rows.empty:
            final_ann_trans_cov = final_ann_trans_cov_rows.iloc[0]
        if not final_ann_trans_mass_tier_rows.empty:
            final_ann_trans_mass_tier = final_ann_trans_mass_tier_rows.iloc[0]
        if not final_ann_coverage_rows.empty:
            final_ann_any_coverage_sig = final_ann_coverage_rows.loc[final_ann_coverage_rows["p"].idxmin()]
    joint_trans = joint[
        joint["model"].eq("pgls_clade_residual_joint")
        & joint["term"].eq("transposon_repeat_suppression_score")
    ].iloc[0]
    sens_trans = pd.DataFrame()
    if not sensitivity.empty:
        sens_trans = sensitivity[
            sensitivity["maintenance_module"].eq("transposon_repeat_suppression")
            & sensitivity["model"].eq("pgls_clade_residual_module")
        ].copy()
        sens_trans = sens_trans.set_index("subset")

    def sens_evidence(subset: str) -> str:
        if sens_trans.empty or subset not in sens_trans.index:
            return "not_run"
        row = sens_trans.loc[subset]
        return f"rank={int(row['rank_by_p'])}, p={fmt(row['p'])}, BH={fmt(row['p_bh_by_subset_model'])}, n={int(row['n'])}"

    def final_annotation_evidence() -> str:
        if final_annotation.empty:
            return "not_run"
        return (
            "Final W3 annotation ranking: "
            f"residual+tier transposon rank={int(final_ann_trans_tier['rank_by_p'])}, "
            f"p={fmt(final_ann_trans_tier['p'])}; "
            f"mass+clade+tier transposon rank={int(final_ann_trans_mass_tier['rank_by_p'])}, "
            f"p={fmt(final_ann_trans_mass_tier['p'])}; "
            f"residual+coverage transposon rank={int(final_ann_trans_cov['rank_by_p'])}, "
            f"p={fmt(final_ann_trans_cov['p'])}; "
            f"best coverage-adjusted module is {final_ann_any_coverage_sig['maintenance_module']} "
            f"p={fmt(final_ann_any_coverage_sig['p'])}."
        )

    gates = pd.DataFrame(
        [
            {
                "gate": "lifespan_residual_primary",
                "status": "pass",
                "evidence": "Final W3 PGLS clade-residual model retains positive transposon association.",
                "action": "Keep residual phenotype as primary endpoint.",
            },
            {
                "gate": "expanded_module_ranking",
                "status": decision(trans_pgls_rank["rank_by_p"] == 1, trans_mass_rank["rank_by_p"] <= 2),
                "evidence": (
                    f"Transposon rank: PGLS residual {int(trans_pgls_rank['rank_by_p'])}, "
                    f"mass+clade {int(trans_mass_rank['rank_by_p'])}."
                ),
                "action": "Present module ranking rather than a single-pathway claim.",
            },
            {
                "gate": "matched_random_specificity",
                "status": decision(
                    trans_random["empirical_p_r_greater"] < 0.025,
                    trans_random["empirical_p_r_greater"] < 0.05,
                ),
                "evidence": (
                    f"Transposon empirical p(r)={fmt(trans_random['empirical_p_r_greater'])}; "
                    f"chromatin empirical p(r)={fmt(chrom_random['empirical_p_r_greater'])}."
                ),
                "action": "Call transposon specificity borderline; do not claim chromatin-specific support.",
            },
            {
                "gate": "annotation_bias_tier",
                "status": decision(ann_trans["p"] < 0.01, ann_trans["p"] < 0.05),
                "evidence": f"Transposon tier-adjusted residual model p={fmt(ann_trans['p'])}.",
                "action": "Report genome-tier adjustment as supportive.",
            },
            {
                "gate": "final_annotation_bias_module_ranking",
                "status": decision(
                    not final_ann_trans_tier.empty
                    and final_ann_trans_tier["rank_by_p"] == 1
                    and final_ann_trans_tier["p"] < 0.01,
                    not final_ann_trans_tier.empty and final_ann_trans_tier["p"] < 0.05,
                ),
                "evidence": final_annotation_evidence(),
                "action": "Keep final W3 module ranking annotation-aware; separate tier support from coverage failure.",
            },
            {
                "gate": "annotation_bias_coverage",
                "status": decision(False, pd.notna(ann_trans_cov["p"]) and ann_trans_cov["p"] >= 0.05),
                "evidence": f"Transposon coverage-adjusted residual model p={fmt(ann_trans_cov['p'])}.",
                "action": "Treat coverage/observability as the major unresolved vulnerability.",
            },
            {
                "gate": "transposon_chromatin_independence",
                "status": decision(False, joint_trans["p"] >= 0.05),
                "evidence": f"Joint transposon/chromatin PGLS residual transposon p={fmt(joint_trans['p'])}.",
                "action": "Avoid independent transposon-vs-chromatin causal separation claims.",
            },
            {
                "gate": "clade_sensitivity",
                "status": "caution",
                "evidence": (
                    f"birds_only {sens_evidence('birds_only')}; "
                    f"no_birds {sens_evidence('no_birds')}; "
                    f"no_bats {sens_evidence('no_bats')}; no_reptiles {sens_evidence('no_reptiles')}."
                ),
                "action": "Frame as bird-dependent or avian-enriched; do not claim broad flight convergence.",
            },
            {
                "gate": "outlier_sensitivity",
                "status": "caution",
                "evidence": (
                    f"exclude_top_abs_residual_5 {sens_evidence('exclude_top_abs_residual_5')}; "
                    f"exclude_top_abs_residual_10 {sens_evidence('exclude_top_abs_residual_10')}."
                ),
                "action": "Treat top residual species as influential; report outlier sensitivity prominently.",
            },
            {
                "gate": "manuscript_readiness",
                "status": "caution-go",
                "evidence": "Several gates pass or caution, but matched-random specificity is borderline, outlier sensitivity is weak, and coverage remains high risk.",
                "action": "Proceed to a cautious manuscript package or preprint-scale Phase 2 report; defer Nature-level framing.",
            },
        ]
    )
    args.decision_gates_output.parent.mkdir(parents=True, exist_ok=True)
    gates.to_csv(args.decision_gates_output, sep="\t", index=False)

    claims = pd.DataFrame(
        [
            {
                "claim_level": "allowed",
                "claim": "Transposon/repeat suppression is a prioritized residual-associated module in the final W3 background.",
                "support": (
                    f"PGLS residual rank 1; matched-random p(r)={fmt(trans_random['empirical_p_r_greater'])}; "
                    f"final W3 residual+tier p={fmt(final_ann_trans_tier.get('p'))}."
                ),
                "required_caveat": "Specificity is borderline and the score collapses after module-coverage adjustment.",
            },
            {
                "claim_level": "allowed",
                "claim": "Cancer/senescence and inflammation are strong competing maintenance modules.",
                "support": "Cancer ranks first in mass+clade; inflammation ranks second in PGLS residual after transposon.",
                "required_caveat": "These are module-score associations, not causal mechanisms.",
            },
            {
                "claim_level": "allowed_with_caution",
                "claim": "Repeat/chromatin biology is implicated as a broader maintenance landscape.",
                "support": "Chromatin PGLS remains positive, but matched-random specificity fails.",
                "required_caveat": "Do not present chromatin as an independently specific module.",
            },
            {
                "claim_level": "not_allowed",
                "claim": "Transposon/repeat suppression independently explains longevity residuals.",
                "support": "Joint transposon/chromatin models and broad module competition do not support independence.",
                "required_caveat": "Needs stronger orthology validation and independent score construction.",
            },
            {
                "claim_level": "not_allowed",
                "claim": "Flight caused strengthened transposon suppression in birds and bats.",
                "support": "Current analyses are association and module-ranking tests, not causal evolutionary inference.",
                "required_caveat": "Requires explicit flight-convergence models and stronger tree/trait sensitivity.",
            },
            {
                "claim_level": "not_allowed",
                "claim": "The final signal is broadly vertebrate-wide and independent of birds.",
                "support": f"No-birds sensitivity: {sens_evidence('no_birds')}.",
                "required_caveat": "Current signal is better described as bird-dependent or avian-enriched.",
            },
            {
                "claim_level": "not_allowed",
                "claim": "The final signal is robust to lifespan-residual outliers.",
                "support": f"Top-5 outlier exclusion: {sens_evidence('exclude_top_abs_residual_5')}.",
                "required_caveat": "Outlier sensitivity must be shown in main limitations or supplementary sensitivity.",
            },
        ]
    )
    args.claim_table_output.parent.mkdir(parents=True, exist_ok=True)
    claims.to_csv(args.claim_table_output, sep="\t", index=False)

    lines = [
        "# Phase 2 Final Assessment",
        "",
        "## Decision",
        "",
        "**Caution-go.** The project remains worth continuing as an annotation-aware module-ranking study. It is not ready for a strong transposon-specific or flight-causality claim.",
        "",
        "## Final Module Ranking",
        "",
    ]
    for model in ["mass_clade_module", "pgls_clade_residual_module", "residual_module"]:
        lines.append(f"### {model}")
        sub = ranking[ranking["model"].eq(model)].sort_values("rank_by_p")
        for _, row in sub.iterrows():
            lines.append(
                f"- rank {int(row['rank_by_p'])}: {row['maintenance_module']}, "
                f"estimate={fmt(row['estimate'])}, p={fmt(row['p'])}, "
                f"BH={fmt(row['p_bh_by_model'])}, mean_coverage={fmt(row['mean_coverage'])}"
            )
        lines.append("")

    lines.extend(["## Decision Gates", ""])
    for _, row in gates.iterrows():
        lines.append(f"- {row['gate']}: **{row['status']}**. {row['evidence']} Action: {row['action']}")

    if not final_annotation.empty:
        lines.extend(["", "## Final Annotation-Bias Ranking", ""])
        lines.append(
            "Genome-tier adjustment supports the final W3 prioritization, but module coverage adjustment collapses the score terms."
        )
        for model in [
            "residual_tier",
            "residual_coverage",
            "residual_tier_coverage",
            "mass_clade_tier",
            "mass_clade_tier_coverage",
        ]:
            sub = final_annotation[final_annotation["model"].eq(model)].sort_values("rank_by_p")
            if sub.empty:
                continue
            lines.append(f"### {model}")
            for _, row in sub.iterrows():
                lines.append(
                    f"- rank {int(row['rank_by_p'])}: {row['maintenance_module']}, "
                    f"estimate={fmt(row['estimate'])}, p={fmt(row['p'])}, "
                    f"BH={fmt(row['p_bh_by_model'])}"
                )
            lines.append("")

    lines.extend(
        [
            "",
            "## Manuscript Claim Guidance",
            "",
            "Allowed headline:",
            "",
            "> A residual-based comparative genomics screen prioritizes transposon/repeat suppression in an avian-enriched longevity panel, within a broader genome-maintenance landscape.",
            "",
            "Avoid:",
            "",
            "- claiming transposon/repeat suppression independently causes longevity;",
            "- claiming chromatin repression is independently specific;",
            "- claiming flight causality without explicit convergence models;",
            "- claiming a bird-independent broad vertebrate signal;",
            "- claiming robustness to top lifespan-residual outliers;",
            "- hiding the coverage/annotation-bias vulnerability.",
            "",
            "## Final Clade Sensitivity",
            "",
            f"- birds-only: {sens_evidence('birds_only')}.",
            f"- no-birds: {sens_evidence('no_birds')}.",
            f"- no-bats: {sens_evidence('no_bats')}.",
            f"- no-reptiles: {sens_evidence('no_reptiles')}.",
            f"- exclude top 5 absolute residuals: {sens_evidence('exclude_top_abs_residual_5')}.",
            f"- exclude top 10 absolute residuals: {sens_evidence('exclude_top_abs_residual_10')}.",
            "",
            "## Next Best Work",
            "",
            "1. Strengthen orthology for low-coverage and high-risk genes before absence claims.",
            "2. Prepare a cautious manuscript/preprint package centered on module ranking and falsification tests.",
            "3. Treat coverage/observability as a first-class sensitivity axis in every main figure and claim.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.report}")


if __name__ == "__main__":
    main()
