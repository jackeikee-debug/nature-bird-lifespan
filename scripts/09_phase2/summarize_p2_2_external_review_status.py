"""Summarize Phase 2 P2.2 external review status and strict scoring implications."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def decision(row: pd.Series) -> str:
    if row["symbol_like_hits"] >= 3:
        return "can_move_to_crossdb_confirmation"
    if row["broad_hits"] >= 2:
        return "manual_review_then_crossdb_confirmation"
    if row["no_hits"] == row["rows_protein"]:
        return "exclude_from_strict_scoring_pending_domain_validation"
    return "external_review_required"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-plan", type=pathlib.Path, required=True)
    parser.add_argument("--protein-summary", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    plan = pd.read_csv(args.gene_plan, sep="\t")
    protein = pd.read_csv(args.protein_summary, sep="\t")
    merged = plan.merge(protein, on="human_gene_symbol", how="left", suffixes=("_pilot", "_protein"))
    merged["p2_2_external_decision"] = merged.apply(decision, axis=1)
    merged["strict_scoring_allowed_now"] = merged["p2_2_external_decision"].isin(
        {"can_move_to_crossdb_confirmation"}
    )
    merged = merged.sort_values(
        [
            "strict_scoring_allowed_now",
            "external_review_priority",
            "symbol_like_fraction",
            "human_gene_symbol",
        ]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)

    decision_counts = merged["p2_2_external_decision"].value_counts().sort_index()
    strict_now = int(merged["strict_scoring_allowed_now"].sum())
    blocked = merged[~merged["strict_scoring_allowed_now"]]
    lines = [
        "# Phase 2 P2.2 External Review Status Report",
        "",
        "## Summary",
        "",
        f"External-review genes: {len(merged)}",
        f"Allowed to move toward strict confirmation now: {strict_now}",
        f"Blocked from strict scoring pending external/domain validation: {len(blocked)}",
        "",
        "## Decision Counts",
        "",
    ]
    for label, count in decision_counts.items():
        lines.append(f"- {label}: {count}")
    lines.extend(
        [
            "",
            "## Blocked Genes",
            "",
        ]
    )
    for _, row in blocked.iterrows():
        lines.append(
            f"- {row['human_gene_symbol']}: {row['p2_2_external_decision']} "
            f"(protein symbol-like hits {int(row['symbol_like_hits'])}/{int(row['rows_protein'])})"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The current external review protects the project from overclaiming. Zero NCBI Gene and zero NCBI Protein symbol-like evidence in the pilot birds should not be scored as biological absence. APOBEC3-family genes and PIWIL3/PIWIL4 require lineage-aware domain or family-level validation before they can support any repeat-suppression claim.",
            "",
            "## Next Action",
            "",
            "Proceed with cross-database/domain confirmation for the few partially supported chromatin genes, while keeping APOBEC3, PIWIL3/4, DNMT3L, TREX1, and ZCCHC3 out of strict v2 scoring until stronger evidence is available.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.report}")


if __name__ == "__main__":
    main()
