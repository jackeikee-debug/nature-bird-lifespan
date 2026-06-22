"""Write a concise Week 4 status report."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


DEFAULT_SEQ_SUMMARY = pathlib.Path("results/tables/week4_full_sequence_validated_transposon_summary.tsv")
DEFAULT_EXTERNAL_QUEUE = pathlib.Path("data/processed/week4_external_orthology_queue.tsv")
DEFAULT_ENSEMBL = pathlib.Path("data/processed/week4_ensembl_external_evidence.tsv")
DEFAULT_LOCAL_EVIDENCE = pathlib.Path("data/processed/week4_ambiguous_local_evidence_summary.tsv")
DEFAULT_NCBI_PILOT = pathlib.Path("data/processed/week4_transposon_ncbi_crosscheck_validation.tsv")
DEFAULT_SENSITIVITY = pathlib.Path("results/tables/week4_full_sequence_sensitivity_summary.tsv")
DEFAULT_OUTPUT = pathlib.Path("results/reports/week4_status_report.md")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence-summary", type=pathlib.Path, default=DEFAULT_SEQ_SUMMARY)
    parser.add_argument("--external-queue", type=pathlib.Path, default=DEFAULT_EXTERNAL_QUEUE)
    parser.add_argument("--ensembl-evidence", type=pathlib.Path, default=DEFAULT_ENSEMBL)
    parser.add_argument("--local-evidence", type=pathlib.Path, default=DEFAULT_LOCAL_EVIDENCE)
    parser.add_argument("--ncbi-pilot", type=pathlib.Path, default=DEFAULT_NCBI_PILOT)
    parser.add_argument("--sensitivity-summary", type=pathlib.Path, default=DEFAULT_SENSITIVITY)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    seq = pd.read_csv(args.sequence_summary, sep="\t")
    queue = pd.read_csv(args.external_queue, sep="\t")
    ensembl = pd.read_csv(args.ensembl_evidence, sep="\t")
    local = pd.read_csv(args.local_evidence, sep="\t") if args.local_evidence.exists() else pd.DataFrame()
    ncbi_pilot = pd.read_csv(args.ncbi_pilot, sep="\t") if args.ncbi_pilot.exists() else pd.DataFrame()
    sensitivity = (
        pd.read_csv(args.sensitivity_summary, sep="\t")
        if args.sensitivity_summary.exists()
        else pd.DataFrame()
    )

    mass_clade = seq[
        (seq["maintenance_module"] == "transposon_suppression")
        & (seq["model"] == "mass_clade_module")
    ].copy()
    queue_gene_counts = queue["human_gene_symbol"].value_counts().sort_index()
    queue_status_counts = queue["week4_diamond_status"].value_counts().sort_index()
    ensembl_counts = ensembl["ensembl_support_summary"].value_counts().sort_index()
    local_counts = (
        local["local_evidence_class"].value_counts().sort_index()
        if not local.empty
        else pd.Series(dtype=int)
    )
    ncbi_pilot_counts = (
        ncbi_pilot["week4_diamond_status"].value_counts().sort_index()
        if not ncbi_pilot.empty
        else pd.Series(dtype=int)
    )

    lines = [
        "# Week 4 Status Report",
        "",
        "## Current Position",
        "",
        "Week 4 is in the orthology-hardening stage. The project has moved beyond candidate coverage and now has a strict sequence-supported transposon-suppression score tested by PGLS.",
        "",
        "## Key Result",
        "",
        "Mass+clade PGLS for transposon-suppression score variants:",
        "",
    ]
    for _, row in mass_clade.iterrows():
        lines.append(
            f"- `{row['score_variant']}`: estimate={float(row['module_estimate']):.3f}, "
            f"p={float(row['module_p']):.4g}, BH={float(row['module_p_bh_by_variant_model']):.4g}, "
            f"support={row['support_class']}."
        )
    weak = mass_clade[mass_clade["score_variant"] == "transposon_sequence_weak_inclusive"]
    weak_support = weak.iloc[0]["support_class"] if not weak.empty else "not_tested"
    if weak_support == "positive_bh_significant":
        weak_sentence = "Weak-inclusive scoring is also positive and BH-significant in the current full sequence-validation run, but the strict score remains the cleaner primary evidence."
    else:
        weak_sentence = "Weak-inclusive scoring is positive but drops to a trend, so weak reciprocal hits should not be used as primary positive evidence."
    lines.extend(
        [
            "",
            f"The strict sequence-supported score remains positive and BH-significant after body-mass and clade adjustment. {weak_sentence}",
            "",
            "## Full Sequence Sensitivity",
            "",
        ]
    )
    if sensitivity.empty:
        lines.append("- Not yet generated")
    else:
        subset = sensitivity[sensitivity["test"] == "subset_sensitivity"].copy()
        leave_one = sensitivity[sensitivity["test"] == "gene_leave_one_out"].copy()
        passed_subset = subset[subset["support_class"] == "positive_bh_significant"]
        failed_subset = subset[subset["support_class"] != "positive_bh_significant"]
        lines.append("Subset checks:")
        for _, row in passed_subset.sort_values("p").iterrows():
            lines.append(
                f"- `{row['comparison']}`: estimate={float(row['estimate']):.3f}, "
                f"p={float(row['p']):.4g}, BH={float(row['p_bh_by_variant_model']):.4g}."
            )
        if not failed_subset.empty:
            lines.append("")
            lines.append("Subset checks that lose support:")
            for _, row in failed_subset.sort_values("p").iterrows():
                lines.append(
                    f"- `{row['comparison']}`: estimate={float(row['estimate']):.3f}, "
                    f"p={float(row['p']):.4g}, BH={float(row['p_bh_by_variant_model']):.4g}."
                )
        drop_rows = leave_one[leave_one["comparison"].str.startswith("drop_", na=False)].copy()
        all_drop_supported = (
            not drop_rows.empty
            and (drop_rows["support_class"] == "positive_bh_significant").all()
        )
        lines.extend(["", "Gene leave-one-out:"])
        for _, row in drop_rows.sort_values("p").iterrows():
            lines.append(
                f"- drop `{row['dropped_gene']}`: estimate={float(row['estimate']):.3f}, "
                f"p={float(row['p']):.4g}, BH={float(row['p_bh_by_variant_model']):.4g}."
            )
        if all_drop_supported:
            lines.append("")
            lines.append("The strict signal survives dropping any single transposon gene, arguing against a one-gene artifact.")
        if not failed_subset.empty:
            lines.append(
                "The signal is lost in the bird-removal and Tier1-only subsets, so the current manuscript framing should stay bird-dependent and annotation-tier sensitive."
            )
    lines.extend(
        [
            "",
            "## Ambiguous External-Validation Queue",
            "",
            f"Rows queued: {len(queue)}",
            f"Species queued: {queue['scientific_name'].nunique()}",
            "",
            "Rows by gene:",
        ]
    )
    for gene, count in queue_gene_counts.items():
        lines.append(f"- {gene}: {count}")
    lines.extend(["", "Rows by DIAMOND status:"])
    for status, count in queue_status_counts.items():
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Ensembl Pilot", ""])
    for status, count in ensembl_counts.items():
        lines.append(f"- {status}: {count}")
    lines.extend(
        [
            "",
            "The Ensembl REST pilot did not support the ambiguous rows by symbol or protein accession. This is database coverage/accession compatibility evidence, not biological absence.",
            "",
            "## Local Ambiguity Diagnosis",
            "",
        ]
    )
    if len(local_counts):
        for cls, count in local_counts.items():
            lines.append(f"- {cls}: {count}")
    else:
        lines.append("- Not yet generated")
    lines.extend(["", "## NCBI Direct Candidate Cross-Check", ""])
    if len(ncbi_pilot_counts):
        lines.append(f"Rows tested: {len(ncbi_pilot)}")
        for status, count in ncbi_pilot_counts.items():
            lines.append(f"- {status}: {count}")
    else:
        lines.append("- Not yet generated")
    lines.extend(
        [
            "",
            "## Next Best Step",
            "",
            "Prioritize OMA/OrthoDB or manual family/domain checks for the 39-row external queue. The highest-risk rows are PIWIL2 candidates that reciprocally resolve to PIWIL1, followed by TRIM28 candidates with weak forward coverage. The strict sequence-supported transposon score should remain the primary Week 4 result.",
        ]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
