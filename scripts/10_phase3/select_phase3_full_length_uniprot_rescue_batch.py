"""Select Phase 3 rows for external full-length UniProt sequence rescue."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


TARGET_CLASSES = {
    "dnmt1_partial_fragment_weak_support",
    "mbd2_mbd3_short_fragment_ambiguity",
    "mbd3_partial_fragment_weak_support",
    "gff_absent_but_uniprot_supported",
    "gff_annotation_without_protein_id",
}

SUPPORT_CALLS = {
    "uniprot_gene_exact_support",
    "uniprot_symbol_like_support",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family-review", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    review = pd.read_csv(args.family_review, sep="\t", dtype=str).fillna("")
    batch = review[
        review["family_domain_review_class"].isin(TARGET_CLASSES)
        & review["uniprot_crosscheck_call"].isin(SUPPORT_CALLS)
        & (review["uniprot_accession_top"] != "")
    ].copy()

    batch["phase3_uniprot_rescue_scope"] = "external_full_length_sequence_sensitivity"
    batch["phase3_uniprot_rescue_reason"] = batch["family_domain_review_class"].map(
        {
            "dnmt1_partial_fragment_weak_support": "Local GFF protein is partial/weak; test whether UniProt provides a longer same-gene sequence.",
            "mbd2_mbd3_short_fragment_ambiguity": "Local MBD2/MBD3 fragment is short; test whether UniProt sequence resolves same-gene support.",
            "mbd3_partial_fragment_weak_support": "Local MBD3 fragment is partial/weak; test whether UniProt provides stronger same-gene support.",
            "gff_absent_but_uniprot_supported": "Local GFF route was absent, but UniProt has target-species support.",
            "gff_annotation_without_protein_id": "Local GFF annotation lacked protein_id, but UniProt has target-species support.",
        }
    )
    batch = batch.sort_values(["family_domain_review_class", "scientific_name", "human_gene_symbol"])

    summary = (
        batch.groupby(["family_domain_review_class", "human_gene_symbol", "uniprot_crosscheck_call"], as_index=False)
        .agg(rows=("scientific_name", "count"), species=("scientific_name", "nunique"))
        .sort_values(["family_domain_review_class", "human_gene_symbol", "uniprot_crosscheck_call"])
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    batch.to_csv(args.output, sep="\t", index=False)
    summary.to_csv(args.summary_output, sep="\t", index=False)

    excluded = review[
        review["family_domain_review_class"].isin(TARGET_CLASSES)
        & ~(
            review["uniprot_crosscheck_call"].isin(SUPPORT_CALLS)
            & (review["uniprot_accession_top"] != "")
        )
    ].copy()
    lines = [
        "# Phase 3 Full-Length UniProt Rescue Batch Report",
        "",
        f"Candidate review classes considered: {len(review[review['family_domain_review_class'].isin(TARGET_CLASSES)])}",
        f"Rows selected for UniProt sequence fetch: {len(batch)}",
        f"Rows excluded because UniProt target support/accession was missing: {len(excluded)}",
        "",
        "## Selected Rows by Review Class",
    ]
    for cls, count in batch["family_domain_review_class"].value_counts().sort_index().items():
        lines.append(f"- {cls}: {count}")
    if batch.empty:
        lines.append("- none: 0")
    lines.extend(
        [
            "",
            "## Interpretation",
            "This is an external-sequence sensitivity branch. Rows selected here should not replace assembly/GFF strict evidence unless reciprocal same-gene validation passes and the manuscript labels the evidence as UniProt-supported.",
            "",
            "## Outputs",
            f"- batch: `{args.output}`",
            f"- summary: `{args.summary_output}`",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
