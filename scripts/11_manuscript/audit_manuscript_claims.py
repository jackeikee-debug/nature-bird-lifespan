#!/usr/bin/env python
"""Audit key manuscript numbers against final analysis tables."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Claim:
    claim_id: str
    value: str
    source: str
    manuscript_pattern: str
    note: str


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def row_by(rows: list[dict[str, str]], **criteria: str) -> dict[str, str]:
    matches = [row for row in rows if all(row.get(key) == value for key, value in criteria.items())]
    if len(matches) != 1:
        raise ValueError(f"Expected one row for {criteria}, found {len(matches)}")
    return matches[0]


def fmt4(value: str) -> str:
    return f"{float(value):.4f}"


def build_claims(root: Path) -> list[Claim]:
    waterfall_path = root / "results/tables/phase3_evidence_waterfall_counts.tsv"
    impact_path = root / "results/tables/phase3_gff_cds_plus_uniprot_sequence_rescue_model_impact.tsv"
    overlay_path = root / "results/tables/phase3_gff_cds_plus_uniprot_combined_overlay_summary.tsv"
    bird_path = root / "results/tables/phase3_gff_cds_plus_uniprot_sequence_rescued_bird_enriched_models.tsv"
    tree_path = root / "results/tables/gene_family_tree_validation_summary.tsv"
    trim_path = root / "results/tables/gene_family_trimmed_alignment_qc.tsv"

    waterfall = read_tsv(waterfall_path)
    counts = {row["waterfall_step"]: row["rows"] for row in waterfall}
    impact = read_tsv(impact_path)
    overlay = read_tsv(overlay_path)
    bird = read_tsv(bird_path)
    trees = read_tsv(tree_path)
    trims = read_tsv(trim_path)

    baseline = row_by(impact, comparison_layer="all_module_pgls", dataset="baseline", model="residual_module")
    rescued = row_by(impact, comparison_layer="all_module_pgls", dataset="rescued", model="residual_module")
    coverage50 = row_by(impact, comparison_layer="high_coverage_subset", dataset="rescued", model="transposon_coverage_ge_0_50")
    birds50 = row_by(impact, comparison_layer="high_coverage_subset", dataset="rescued", model="birds_only_transposon_coverage_ge_0_50")
    interaction = row_by(
        bird,
        model="residual_bird_interaction",
        term="transposon_repeat_suppression_score:bird_statusbird",
    )
    dnmt_full = row_by(trees, family="DNMT_family", tree_variant="full_alignment")
    dnmt_trim = row_by(trees, family="DNMT_family", tree_variant="trimmed_gappy70")
    mbd_full = row_by(trees, family="MBD_family", tree_variant="full_alignment")
    mbd_trim = row_by(trees, family="MBD_family", tree_variant="trimmed_gappy70")
    dnmt_qc = row_by(trims, family="DNMT_family")
    mbd_qc = row_by(trims, family="MBD_family")

    return [
        Claim("validation_input_rows", counts["input_priority1_rows"], str(waterfall_path), r"\b140\b", "High-priority validation rows"),
        Claim("gff_annotation_rows", counts["gff_annotation_rescue"], str(waterfall_path), r"\b91\b", "Rows with GFF annotation support"),
        Claim("strict_gff_rows", counts["local_gff_protein_strict"], str(waterfall_path), r"\b58\b", "Strict local GFF protein rows"),
        Claim("strict_local_total", counts["strict_local_sequence_total"], str(waterfall_path), r"\b62\b", "Total strict local rows"),
        Claim("strict_any_total", counts["strict_any_sequence_total"], str(waterfall_path), r"\b65\b", "Strict rows including external sensitivity"),
        Claim("partial_family_rows", counts["partial_or_family_not_scoreable"], str(waterfall_path), r"\b16\b", "Partial/family not-scoreable rows"),
        Claim("unresolved_rows", counts["not_found_or_unresolved"], str(waterfall_path), r"\b45\b", "Unresolved rows"),
        Claim("matrix_updates", str(len(overlay)), str(overlay_path), r"\bseven strict sequence-supported updates\b", "Validation-overlay updates"),
        Claim("baseline_residual_estimate", fmt4(baseline["estimate"]), str(impact_path), re.escape(fmt4(baseline["estimate"])), "Baseline residual estimate"),
        Claim("rescued_residual_estimate", fmt4(rescued["estimate"]), str(impact_path), re.escape(fmt4(rescued["estimate"])), "Validation-overlay residual estimate"),
        Claim("coverage50_p", fmt4(coverage50["p"]), str(impact_path), re.escape(fmt4(coverage50["p"])), "Coverage >=0.50 residual-model P value"),
        Claim("birds50_p", fmt4(birds50["p"]), str(impact_path), re.escape(fmt4(birds50["p"])), "Birds-only coverage >=0.50 P value"),
        Claim("bird_interaction_p", fmt4(interaction["p"]), str(bird_path), re.escape(fmt4(interaction["p"])), "Residual bird-interaction P value"),
        Claim("dnmt_sequences", dnmt_full["sequences"], str(tree_path), r"DNMT family analysis contained 33 sequences", "DNMT tree sequence count"),
        Claim("dnmt_trimmed_columns", dnmt_qc["trimmed_columns"], str(trim_path), r"860 of 2,067 alignment columns", "DNMT retained columns"),
        Claim("mbd_sequences", mbd_full["sequences"], str(tree_path), r"MBD analysis contained 22 sequences", "MBD tree sequence count"),
        Claim("mbd_trimmed_columns", mbd_qc["trimmed_columns"], str(trim_path), r"233 of 437 alignment columns", "MBD retained columns"),
        Claim("dnmt_full_convergence", "0.993", str(tree_path), r"0\.993", "DNMT full-alignment convergence diagnostic"),
        Claim("dnmt_trimmed_convergence", "0.998", str(tree_path), r"0\.998", "DNMT trimmed-alignment convergence diagnostic"),
        Claim("mbd_full_nonconvergence", "not_converged", str(tree_path), r"full-alignment bootstrap analysis did not converge", "MBD full-alignment warning"),
        Claim("mbd_trimmed_convergence", "0.996", str(tree_path), r"0\.996", "MBD trimmed-alignment convergence diagnostic"),
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manuscript", default="docs/manuscript_draft_v1.md")
    parser.add_argument("--table-output", default="results/tables/manuscript_claim_audit.tsv")
    parser.add_argument("--report", default="results/reports/manuscript_claim_audit_report.md")
    args = parser.parse_args()

    root = Path.cwd()
    manuscript_path = root / args.manuscript
    manuscript = manuscript_path.read_text(encoding="utf-8")
    claims = build_claims(root)
    rows = []
    for claim in claims:
        present = bool(re.search(claim.manuscript_pattern, manuscript, flags=re.IGNORECASE))
        rows.append(
            {
                "claim_id": claim.claim_id,
                "expected_value": claim.value,
                "source_file": str(Path(claim.source).relative_to(root)).replace("\\", "/"),
                "manuscript_pattern": claim.manuscript_pattern,
                "status": "pass" if present else "fail",
                "note": claim.note,
            }
        )

    table_path = root / args.table_output
    table_path.parent.mkdir(parents=True, exist_ok=True)
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    failures = [row for row in rows if row["status"] == "fail"]
    report_path = root / args.report
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_lines = [
        "# Manuscript Claim Audit",
        "",
        f"- Manuscript: `{args.manuscript}`",
        f"- Claims checked: {len(rows)}",
        f"- Passed: {len(rows) - len(failures)}",
        f"- Failed: {len(failures)}",
        "",
        "The audit checks whether key final-table values and gene-tree caveats are represented in the manuscript draft. It does not replace scientific review of wording or model interpretation.",
    ]
    if failures:
        report_lines.extend(["", "## Missing or stale claims", ""])
        report_lines.extend(f"- `{row['claim_id']}`: expected `{row['expected_value']}`" for row in failures)
    else:
        report_lines.extend(["", "All audited claims are consistent with the current final tables."])
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(f"Audited {len(rows)} claims: {len(failures)} failures")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
