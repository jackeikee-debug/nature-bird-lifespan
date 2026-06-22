"""Build row-level Phase 3 rescue decisions for priority-1 batch01."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def join_unique(values: pd.Series) -> str:
    vals = [str(v) for v in values if str(v) and str(v) != "nan"]
    return ";".join(sorted(set(vals)))


def decision(row: pd.Series) -> tuple[str, str]:
    status = row.get("ncbi_protein_review_status", "")
    reciprocal_tools = int(row.get("reciprocal_supported_tools", 0) or 0)
    forward_tools = int(row.get("forward_supported_tools", 0) or 0)
    weak_tools = int(row.get("weak_same_gene_tools", 0) or 0)
    no_same_tools = int(row.get("no_same_gene_tools", 0) or 0)
    top_genes = set(str(row.get("top_reference_genes", "")).split(";")) - {""}
    gene = str(row.get("human_gene_symbol", ""))

    if status == "no_protein_hit" or not status:
        return "not_rescued_no_protein_hit", "NCBI Protein search found no candidate sequence."
    if reciprocal_tools >= 2:
        return "sequence_supported_rescue", "Both DIAMOND and BLASTP support reciprocal same-gene assignment."
    if reciprocal_tools == 1 and forward_tools >= 1:
        return "sequence_supported_manual_review", "One tool is reciprocal and another supports same-gene forward matching."
    if reciprocal_tools == 1:
        return "sequence_supported_one_tool_review", "One tool supports reciprocal same-gene assignment; retain for manual review."
    if weak_tools >= 2:
        return "partial_same_gene_support_not_strict", "Both tools find same-gene hits but below coverage thresholds, consistent with partial proteins."
    if no_same_tools >= 1 and top_genes and gene not in top_genes:
        return "reject_probable_paralog_or_wrong_gene", f"Best reference hit is {join_unique(pd.Series(list(top_genes)))} rather than {gene}."
    if status == "protein_broad_hit":
        return "protein_broad_hit_unresolved", "Protein hit exists but reciprocal sequence support is insufficient."
    return "manual_review_required", "Insufficient evidence for automatic rescue decision."


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=pathlib.Path, required=True)
    parser.add_argument("--protein-review", type=pathlib.Path, required=True)
    parser.add_argument("--validation", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    batch = pd.read_csv(args.batch, sep="\t", dtype=str).fillna("")
    protein = pd.read_csv(args.protein_review, sep="\t", dtype=str).fillna("")
    if args.validation.exists() and args.validation.stat().st_size > 0:
        validation = pd.read_csv(args.validation, sep="\t", dtype=str).fillna("")
    else:
        validation = pd.DataFrame()

    key = ["human_gene_symbol", "scientific_name"]
    protein_cols = key + [
        "ncbi_protein_review_status",
        "ncbi_protein_hit_count",
        "ncbi_protein_accessions_top",
        "ncbi_protein_titles_top",
    ]
    merged = batch.merge(protein[protein_cols], on=key, how="left")

    if validation.empty:
        agg = pd.DataFrame(columns=key)
    else:
        validation["is_reciprocal"] = validation["validation_call"].eq("reciprocal_same_gene_supported")
        validation["is_forward"] = validation["validation_call"].eq("same_gene_forward_supported_reciprocal_weak")
        validation["is_weak"] = validation["validation_call"].eq("weak_same_gene_support")
        validation["is_no_same"] = validation["validation_call"].eq("no_same_gene_reference_hit")
        agg = validation.groupby(key, as_index=False).agg(
            validation_tools=("tool", join_unique),
            reciprocal_supported_tools=("is_reciprocal", "sum"),
            forward_supported_tools=("is_forward", "sum"),
            weak_same_gene_tools=("is_weak", "sum"),
            no_same_gene_tools=("is_no_same", "sum"),
            candidate_accessions_validated=("candidate_protein_accession", join_unique),
            top_reference_genes=("top_reference_gene", join_unique),
            reciprocal_best_genes=("reciprocal_best_gene", join_unique),
            validation_calls=("validation_call", join_unique),
            median_pident=("pident", lambda x: pd.to_numeric(x, errors="coerce").median()),
            median_qcovhsp=("qcovhsp", lambda x: pd.to_numeric(x, errors="coerce").median()),
            median_scovhsp=("scovhsp", lambda x: pd.to_numeric(x, errors="coerce").median()),
        )

    merged = merged.merge(agg, on=key, how="left")
    for col in [
        "reciprocal_supported_tools",
        "forward_supported_tools",
        "weak_same_gene_tools",
        "no_same_gene_tools",
    ]:
        if col not in merged.columns:
            merged[col] = 0
        merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0).astype(int)
    for col in [
        "validation_tools",
        "candidate_accessions_validated",
        "top_reference_genes",
        "reciprocal_best_genes",
        "validation_calls",
    ]:
        if col not in merged.columns:
            merged[col] = ""
        merged[col] = merged[col].fillna("")

    decisions = merged.apply(decision, axis=1, result_type="expand")
    merged["phase3_rescue_decision"] = decisions[0]
    merged["phase3_rescue_reason"] = decisions[1]
    merged["can_count_as_rescued_for_coverage"] = merged["phase3_rescue_decision"].isin(
        [
            "sequence_supported_rescue",
            "sequence_supported_manual_review",
            "sequence_supported_one_tool_review",
            "partial_same_gene_support_not_strict",
        ]
    )
    merged["can_count_as_strict_sequence"] = merged["phase3_rescue_decision"].eq("sequence_supported_rescue")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)

    summary = (
        merged.groupby(["phase3_rescue_decision"], as_index=False)
        .agg(
            rows=("human_gene_symbol", "count"),
            genes=("human_gene_symbol", "nunique"),
            species=("scientific_name", "nunique"),
        )
        .sort_values("rows", ascending=False)
    )
    summary.to_csv(args.summary_output, sep="\t", index=False)

    species_summary = merged.groupby("scientific_name", as_index=False).agg(
        rows=("human_gene_symbol", "count"),
        rescued_rows=("can_count_as_rescued_for_coverage", "sum"),
        strict_rows=("can_count_as_strict_sequence", "sum"),
        no_hit_rows=("phase3_rescue_decision", lambda s: int((s == "not_rescued_no_protein_hit").sum())),
        rejected_rows=("phase3_rescue_decision", lambda s: int((s == "reject_probable_paralog_or_wrong_gene").sum())),
    )
    batch_ids = sorted(set(map(str, merged.get("phase3_batch_id", pd.Series(["batch"])))) - {""})
    batch_label = batch_ids[0] if len(batch_ids) == 1 else "combined_batch"
    batch_label = batch_label.replace("phase3_priority1_transposon_bird_", "").upper()
    lines = [
        f"# Phase 3 Priority-1 {batch_label} Rescue Decision Report",
        "",
        f"Rows assessed: {len(merged)}",
        f"Rows that can count as coverage rescue: {int(merged['can_count_as_rescued_for_coverage'].sum())}",
        f"Rows that can count as strict sequence rescue: {int(merged['can_count_as_strict_sequence'].sum())}",
        "",
        "## Decision Counts",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(f"- {row['phase3_rescue_decision']}: {int(row['rows'])}")
    lines.extend(["", "## Species Summary", ""])
    for _, row in species_summary.iterrows():
        lines.append(
            f"- {row['scientific_name']}: rescued={int(row['rescued_rows'])}/{int(row['rows'])}, "
            f"strict={int(row['strict_rows'])}, rejected={int(row['rejected_rows'])}, no_hit={int(row['no_hit_rows'])}"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This batch shows whether low-coverage birds are recoverable through protein-level evidence. Strict rescue requires reciprocal support from both DIAMOND and BLASTP; partial proteins and one-tool cases can improve coverage sensitivity but should not be used for strict absence or mechanism claims.",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
