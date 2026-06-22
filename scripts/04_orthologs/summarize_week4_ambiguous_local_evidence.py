"""Summarize local evidence for ambiguous Week 4 transposon rows."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DEFAULT_QUEUE = pathlib.Path("data/processed/week4_external_orthology_queue.tsv")
DEFAULT_RECIP = pathlib.Path("results/tables/week4_transposon_reciprocal_diamond_hits.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/week4_ambiguous_local_evidence_summary.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week4_ambiguous_local_evidence_summary_report.md")


def top_reciprocal_by_gene(recip: pd.DataFrame) -> pd.DataFrame:
    recip = recip.copy()
    for col in ["bitscore", "pident", "qcovhsp", "scovhsp"]:
        recip[col] = pd.to_numeric(recip[col], errors="coerce")
    rows = []
    for (accession, query_gene, target), group in recip.groupby(
        ["best_assembly_accession", "query_gene", "target_protein_id"], dropna=False
    ):
        ordered = group.sort_values("bitscore", ascending=False)
        top = ordered.iloc[0]
        expected = group[group["reciprocal_gene"] == query_gene].sort_values("bitscore", ascending=False)
        expected_top = expected.iloc[0] if not expected.empty else None
        rows.append(
            {
                "best_assembly_accession": accession,
                "human_gene_symbol": query_gene,
                "forward_target_protein_id": target,
                "reciprocal_top_gene": top["reciprocal_gene"],
                "reciprocal_top_reference": top["sseqid"],
                "reciprocal_top_bitscore": top["bitscore"],
                "reciprocal_top_pident": top["pident"],
                "reciprocal_top_qcovhsp": top["qcovhsp"],
                "reciprocal_top_scovhsp": top["scovhsp"],
                "expected_gene_best_reference": expected_top["sseqid"] if expected_top is not None else "",
                "expected_gene_best_bitscore": expected_top["bitscore"] if expected_top is not None else pd.NA,
                "top_minus_expected_bitscore": (
                    top["bitscore"] - expected_top["bitscore"] if expected_top is not None else pd.NA
                ),
            }
        )
    return pd.DataFrame(rows)


def classify(row: pd.Series) -> str:
    gene = row["human_gene_symbol"]
    status = row["week4_diamond_status"]
    reciprocal_top = row.get("reciprocal_top_gene", "")
    margin = pd.to_numeric(row.get("top_minus_expected_bitscore"), errors="coerce")
    qcov = pd.to_numeric(row.get("forward_qcovhsp"), errors="coerce")
    scov = pd.to_numeric(row.get("forward_scovhsp"), errors="coerce")
    pident = pd.to_numeric(row.get("forward_pident"), errors="coerce")

    if gene == "PIWIL2" and reciprocal_top == "PIWIL1":
        if pd.notna(margin) and margin > 100:
            return "likely_PIWIL1_paralog_not_PIWIL2"
        return "PIWI_family_ambiguous"
    if gene == "TRIM28":
        if status == "not_validated":
            return "no_local_sequence_support"
        if pd.notna(qcov) and pd.notna(scov) and qcov < 50 and scov < 50 and pd.notna(pident) and pident >= 30:
            return "partial_TRIM28_like_low_coverage"
        return "TRIM28_like_needs_domain_check"
    if status == "reciprocal_weak":
        return "weak_reciprocal_support_needs_external_check"
    if status == "weak_forward_support":
        return "weak_forward_support_needs_external_check"
    return "ambiguous_needs_external_check"


def recommendation(row: pd.Series) -> str:
    cls = row["local_evidence_class"]
    if cls == "likely_PIWIL1_paralog_not_PIWIL2":
        return "do_not_count_as_PIWIL2_without_family_tree_or_external_database_support"
    if cls == "PIWI_family_ambiguous":
        return "check_PIWI_family_tree_OMA_OrthoDB_before_counting"
    if cls == "partial_TRIM28_like_low_coverage":
        return "hold_out_from_primary_positive_set_check_domain_architecture"
    if cls == "no_local_sequence_support":
        return "keep_unresolved_unless_external_database_supports"
    return "cross_check_OMA_OrthoDB_or_manual_domain_evidence"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=pathlib.Path, default=DEFAULT_QUEUE)
    parser.add_argument("--reciprocal", type=pathlib.Path, default=DEFAULT_RECIP)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    queue = pd.read_csv(args.queue, sep="\t")
    recip = pd.read_csv(args.reciprocal, sep="\t")
    recip_summary = top_reciprocal_by_gene(recip)
    merged = queue.merge(
        recip_summary,
        on=["best_assembly_accession", "human_gene_symbol", "forward_target_protein_id"],
        how="left",
        suffixes=("", "_local"),
    )
    if "reciprocal_top_reference_local" in merged.columns:
        merged["local_reciprocal_top_reference"] = merged["reciprocal_top_reference_local"]
    else:
        merged["local_reciprocal_top_reference"] = ""
    merged["local_evidence_class"] = merged.apply(classify, axis=1)
    merged["recommended_handling"] = merged.apply(recommendation, axis=1)
    preferred = [
        "external_validation_priority",
        "human_gene_symbol",
        "scientific_name",
        "clade",
        "flight_status",
        "genome_analysis_tier",
        "species_taxid",
        "best_assembly_accession",
        "final_candidate_status",
        "week4_diamond_status",
        "week4_diamond_reason",
        "forward_target_protein_id",
        "forward_pident",
        "forward_qcovhsp",
        "forward_scovhsp",
        "forward_bitscore",
        "reciprocal_gene",
        "reciprocal_top_reference",
        "reciprocal_top_gene",
        "local_reciprocal_top_reference",
        "reciprocal_top_bitscore",
        "expected_gene_best_reference",
        "expected_gene_best_bitscore",
        "top_minus_expected_bitscore",
        "local_evidence_class",
        "recommended_handling",
        "suggested_external_sources",
        "external_validation_question",
    ]
    extras = [col for col in merged.columns if col not in preferred and not col.endswith("_local")]
    merged = merged[[col for col in preferred if col in merged.columns] + extras]
    merged = merged.sort_values(
        ["local_evidence_class", "human_gene_symbol", "scientific_name"]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)

    class_counts = merged["local_evidence_class"].value_counts().sort_index()
    gene_class = merged.groupby(["human_gene_symbol", "local_evidence_class"]).size().reset_index(name="rows")
    lines = [
        "# Week 4 Ambiguous Local Evidence Summary Report",
        "",
        f"Rows summarized: {len(merged)}",
        f"Species: {merged['scientific_name'].nunique()}",
        "",
        "## Evidence Classes",
        "",
    ]
    for cls, count in class_counts.items():
        lines.append(f"- {cls}: {count}")
    lines.extend(["", "## Gene x Evidence Class", ""])
    for _, rec in gene_class.iterrows():
        lines.append(f"- {rec['human_gene_symbol']} / {rec['local_evidence_class']}: {rec['rows']}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "PIWIL2 ambiguity is mostly paralog ambiguity: many candidates reciprocally resolve to PIWIL1 with a large bitscore margin. TRIM28 ambiguity is mostly low-coverage partial similarity, so these rows should remain outside the primary positive set until domain or external orthology evidence supports them.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
