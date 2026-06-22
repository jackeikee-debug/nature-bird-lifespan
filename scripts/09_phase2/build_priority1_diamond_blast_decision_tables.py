"""Build strict/domain/manual tables from DIAMOND/BLAST reciprocal validation."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-summary", type=pathlib.Path, required=True)
    parser.add_argument("--row-validation", type=pathlib.Path, required=True)
    parser.add_argument("--strict-output", type=pathlib.Path, required=True)
    parser.add_argument("--domain-output", type=pathlib.Path, required=True)
    parser.add_argument("--manual-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    genes = pd.read_csv(args.gene_summary, sep="\t")
    rows = pd.read_csv(args.row_validation, sep="\t")

    strict = genes[
        genes["diamond_blast_upgrade_decision"] == "strict_upgrade_candidate_diamond_blast_supported"
    ].copy()
    strict["strict_upgrade_basis"] = "NCBI_Gene_plus_DIAMOND_and_BLASTP_reciprocal_same_gene_support"
    strict["strict_upgrade_allowed"] = True
    strict["absence_claim_allowed"] = False

    domain = genes[
        genes["diamond_blast_upgrade_decision"] == "domain_paralog_check_required_after_reciprocal_support"
    ].copy()
    domain["domain_check_basis"] = "DIAMOND_BLAST_reciprocal_support_but_TDRD_paralog_family"
    domain["strict_upgrade_allowed"] = False
    domain["absence_claim_allowed"] = False

    manual = genes[
        ~genes["diamond_blast_upgrade_decision"].isin(
            {
                "strict_upgrade_candidate_diamond_blast_supported",
                "domain_paralog_check_required_after_reciprocal_support",
            }
        )
    ].copy()
    manual["manual_review_basis"] = "insufficient_species_count_or_mixed_reciprocal_support"
    manual["strict_upgrade_allowed"] = False
    manual["absence_claim_allowed"] = False

    args.strict_output.parent.mkdir(parents=True, exist_ok=True)
    args.domain_output.parent.mkdir(parents=True, exist_ok=True)
    args.manual_output.parent.mkdir(parents=True, exist_ok=True)
    strict.to_csv(args.strict_output, sep="\t", index=False)
    domain.to_csv(args.domain_output, sep="\t", index=False)
    manual.to_csv(args.manual_output, sep="\t", index=False)

    row_counts = rows.groupby(["tool", "validation_call"]).size().reset_index(name="rows")
    lines = [
        "# Phase 2 Priority 1 DIAMOND/BLAST Decision Tables Report",
        "",
        "## Summary",
        "",
        f"Genes assessed: {genes['human_gene_symbol'].nunique()}",
        f"Strict-upgrade candidates: {strict['human_gene_symbol'].nunique()}",
        f"Domain/paralog check genes: {domain['human_gene_symbol'].nunique()}",
        f"Manual-review/no-absence genes: {manual['human_gene_symbol'].nunique()}",
        "",
        "## Strict Upgrade Candidates",
        "",
    ]
    for gene in strict["human_gene_symbol"].tolist():
        lines.append(f"- {gene}")
    lines.extend(["", "## Domain/Paralog Check Queue", ""])
    for gene in domain["human_gene_symbol"].tolist():
        lines.append(f"- {gene}")
    lines.extend(["", "## Manual Review / No Absence Claim", ""])
    for gene in manual["human_gene_symbol"].tolist():
        lines.append(f"- {gene}")
    lines.extend(["", "## Row-Level Calls", ""])
    for _, row in row_counts.iterrows():
        lines.append(f"- {row['tool']} {row['validation_call']}: {row['rows']}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "These tables supersede the lightweight Biopython sequence-screening tables for priority 1 strict-upgrade decisions. Only the strict-upgrade candidates should be moved into the strict v2 gene set; TDRD-family rows remain domain/paralog-check pending, and mixed-support genes remain protected from absence claims.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.strict_output}")


if __name__ == "__main__":
    main()
