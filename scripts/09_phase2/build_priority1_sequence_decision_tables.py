"""Build strict-upgrade and review tables from priority 1 sequence confirmation."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-summary", type=pathlib.Path, required=True)
    parser.add_argument("--row-confirmation", type=pathlib.Path, required=True)
    parser.add_argument("--strict-output", type=pathlib.Path, required=True)
    parser.add_argument("--domain-output", type=pathlib.Path, required=True)
    parser.add_argument("--manual-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    genes = pd.read_csv(args.gene_summary, sep="\t")
    rows = pd.read_csv(args.row_confirmation, sep="\t")

    strict = genes[
        genes["sequence_upgrade_decision"] == "strict_upgrade_candidate_sequence_supported"
    ].copy()
    strict["strict_upgrade_basis"] = "NCBI_Gene_plus_NCBI_protein_sequence_similarity_to_human_reference"
    strict["strict_upgrade_allowed"] = True

    domain = genes[
        genes["sequence_upgrade_decision"] == "domain_architecture_check_required_before_strict_upgrade"
    ].copy()
    domain["domain_check_basis"] = "sequence_supported_but_high_paralog_family"
    domain["strict_upgrade_allowed"] = False

    manual = genes[
        ~genes["sequence_upgrade_decision"].isin(
            {
                "strict_upgrade_candidate_sequence_supported",
                "domain_architecture_check_required_before_strict_upgrade",
            }
        )
    ].copy()
    manual["absence_claim_allowed"] = False
    manual["manual_review_basis"] = "insufficient_sequence_coverage_or_mixed_support"

    args.strict_output.parent.mkdir(parents=True, exist_ok=True)
    args.domain_output.parent.mkdir(parents=True, exist_ok=True)
    args.manual_output.parent.mkdir(parents=True, exist_ok=True)
    strict.to_csv(args.strict_output, sep="\t", index=False)
    domain.to_csv(args.domain_output, sep="\t", index=False)
    manual.to_csv(args.manual_output, sep="\t", index=False)

    row_counts = rows["sequence_confirmation_call"].value_counts().sort_index()
    lines = [
        "# Phase 2 Priority 1 Sequence Decision Tables Report",
        "",
        "## Summary",
        "",
        f"Genes assessed: {genes['human_gene_symbol'].nunique()}",
        f"Strict-upgrade candidates: {strict['human_gene_symbol'].nunique()}",
        f"Domain/paralog check queue genes: {domain['human_gene_symbol'].nunique()}",
        f"Manual-review/no-absence genes: {manual['human_gene_symbol'].nunique()}",
        "",
        "## Row-Level Calls",
        "",
    ]
    for call, count in row_counts.items():
        lines.append(f"- {call}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Only the strict-upgrade table can be considered for moving priority 1 genes into strict v2. Domain/paralog and manual-review genes remain excluded from strict upgrade and must not be interpreted as absences.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.strict_output}")


if __name__ == "__main__":
    main()
