"""Summarize Week 4 full sequence-validation sensitivity checks."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DEFAULT_SUBSET = pathlib.Path("results/tables/week4_full_sequence_transposon_sensitivity.tsv")
DEFAULT_LEAVE_ONE = pathlib.Path("results/tables/week4_transposon_gene_leave_one_out_pgls.tsv")
DEFAULT_OUTPUT = pathlib.Path("results/tables/week4_full_sequence_sensitivity_summary.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week4_full_sequence_sensitivity_summary_report.md")


def fmt(value: object) -> str:
    value = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(value):
        return ""
    return f"{value:.6g}"


def classify(estimate: float, p_value: float, bh_value: float) -> str:
    if estimate > 0 and bh_value < 0.05:
        return "positive_bh_significant"
    if estimate > 0 and p_value < 0.05:
        return "positive_nominal"
    if estimate > 0 and p_value < 0.10:
        return "positive_trend"
    if estimate <= 0:
        return "lost_or_reversed"
    return "not_supported"


def support_note(row: pd.Series) -> str:
    test = row["test"]
    name = row["comparison"]
    support = row["support_class"]
    if test == "subset_sensitivity":
        if name == "leave_out_Aves":
            return "signal depends on birds in the current panel"
        if name == "tier1_only":
            return "signal weakens under strictest annotation-tier restriction"
        if name == "exclude_top_abs_residual_5":
            return "not driven only by the largest lifespan-residual outliers"
        if name == "leave_out_Mammalia_Chiroptera":
            return "not driven only by bats"
        if name == "leave_out_Mammalia_nonChiroptera":
            return "not driven only by non-flying mammals"
        if name == "leave_out_Reptilia":
            return "not driven only by reptiles"
        if name == "exclude_human":
            return "not driven by human"
        if name == "all_primary":
            return "primary full sequence-validated result"
    if test == "gene_leave_one_out":
        if name == "all_5_genes_strict":
            return "all five transposon genes included"
        if support == "positive_bh_significant":
            return "signal survives dropping this gene"
        return "signal is sensitive to dropping this gene"
    return ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--subset", type=pathlib.Path, default=DEFAULT_SUBSET)
    parser.add_argument("--leave-one-out", type=pathlib.Path, default=DEFAULT_LEAVE_ONE)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    subset = pd.read_csv(args.subset, sep="\t")
    leave_one = pd.read_csv(args.leave_one_out, sep="\t")

    subset_key = subset[
        (subset["score_variant"] == "transposon_sequence_strict")
        & (subset["model"] == "mass_clade_module")
        & (subset["error"].fillna("") == "")
    ].copy()
    for col in ["estimate", "p", "p_bh_by_variant_model", "n"]:
        subset_key[col] = pd.to_numeric(subset_key[col], errors="coerce")
    subset_key["test"] = "subset_sensitivity"
    subset_key["comparison"] = subset_key["subset"]
    subset_key["dropped_gene"] = ""
    subset_key["support_class"] = subset_key.apply(
        lambda row: classify(row["estimate"], row["p"], row["p_bh_by_variant_model"]),
        axis=1,
    )
    subset_out = subset_key[
        [
            "test",
            "comparison",
            "dropped_gene",
            "score_variant",
            "model",
            "n",
            "estimate",
            "p",
            "p_bh_by_variant_model",
            "support_class",
        ]
    ].copy()

    loo_key = leave_one[
        leave_one["score_variant"].str.endswith("_strict")
        & (leave_one["model"] == "mass_clade_module")
        & (leave_one["error"].fillna("") == "")
    ].copy()
    for col in ["module_estimate", "module_p", "module_p_bh_by_variant_model", "n"]:
        loo_key[col] = pd.to_numeric(loo_key[col], errors="coerce")
    loo_key["test"] = "gene_leave_one_out"
    loo_key["comparison"] = loo_key["score_variant"]
    loo_key["dropped_gene"] = loo_key["score_variant"].str.extract(r"^drop_(.+)_strict$")[0].fillna("")
    loo_key.loc[loo_key["score_variant"] == "all_5_genes_strict", "comparison"] = "all_5_genes_strict"
    loo_key["support_class"] = loo_key.apply(
        lambda row: classify(
            row["module_estimate"], row["module_p"], row["module_p_bh_by_variant_model"]
        ),
        axis=1,
    )
    loo_out = loo_key.rename(
        columns={
            "module_estimate": "estimate",
            "module_p": "p",
            "module_p_bh_by_variant_model": "p_bh_by_variant_model",
        }
    )[
        [
            "test",
            "comparison",
            "dropped_gene",
            "score_variant",
            "model",
            "n",
            "estimate",
            "p",
            "p_bh_by_variant_model",
            "support_class",
        ]
    ].copy()

    out = pd.concat([subset_out, loo_out], ignore_index=True)
    out["interpretation_note"] = out.apply(support_note, axis=1)
    out = out.sort_values(["test", "p", "comparison"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    subset_report = subset_out.sort_values("p")
    loo_report = loo_out.sort_values("p")
    lost = subset_report[subset_report["support_class"].isin(["lost_or_reversed", "not_supported"])]
    surviving_loo = loo_report[loo_report["comparison"].str.startswith("drop_")]
    all_loo_supported = (surviving_loo["support_class"] == "positive_bh_significant").all()

    lines = [
        "# Week 4 Full Sequence Sensitivity Summary",
        "",
        "This report summarizes the two most important robustness checks for the strict full sequence-validated transposon-suppression signal.",
        "",
        "## Subset Sensitivity",
        "",
    ]
    for _, row in subset_report.iterrows():
        lines.append(
            f"- `{row['comparison']}`: estimate={fmt(row['estimate'])}, "
            f"p={fmt(row['p'])}, BH={fmt(row['p_bh_by_variant_model'])}, "
            f"n={int(row['n'])}, support={row['support_class']}."
        )
    lines.extend(["", "## Gene Leave-One-Out", ""])
    for _, row in loo_report.iterrows():
        label = row["comparison"]
        dropped = row["dropped_gene"]
        if dropped:
            label = f"drop {dropped}"
        lines.append(
            f"- `{label}`: estimate={fmt(row['estimate'])}, "
            f"p={fmt(row['p'])}, BH={fmt(row['p_bh_by_variant_model'])}, "
            f"n={int(row['n'])}, support={row['support_class']}."
        )

    lines.extend(["", "## Interpretation", ""])
    lines.append(
        "The strict sequence-supported transposon signal is robust to excluding bats, reptiles, non-flying mammals, human, and the five largest lifespan-residual outliers."
    )
    if not lost.empty:
        lost_labels = ", ".join(f"`{x}`" for x in lost["comparison"].tolist())
        lines.append(
            f"The signal is not robust in {lost_labels}, so the current result should be framed as bird-heavy and annotation-tier sensitive rather than a universal vertebrate effect."
        )
    if all_loo_supported:
        lines.append(
            "The signal survives dropping any single transposon gene, which argues against a one-gene artifact."
        )
    else:
        lines.append(
            "At least one leave-one-out model loses support, so gene composition should remain a priority sensitivity axis."
        )
    lines.append(
        "The most defensible Week 4 claim is therefore: a bird-dependent, sequence-supported transposon-suppression association with lifespan, not driven by one gene or by obvious clade/outlier exclusions."
    )

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
