"""Build Phase 3 evidence ladder tables and validation waterfall figures."""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.pyplot as plt
import pandas as pd


def read(path: pathlib.Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def first_nonempty(*values: object) -> str:
    for value in values:
        text = str(value)
        if text and text != "nan":
            return text
    return ""


def classify(row: pd.Series) -> tuple[str, str, str, str]:
    gff_decision = row.get("phase3_gff_sequence_decision", "")
    cds_decision = row.get("phase3_cds_translation_decision", "")
    uniprot_decision = row.get("phase3_uniprot_sequence_decision", "")
    partial_class = row.get("phase3_partial_family_resolution_class", "")
    gff_call = row.get("gff_rescue_call", "")
    gff_annotation = row.get("gff_can_count_as_annotation_rescue", "")

    if gff_decision == "gff_sequence_strict_rescue":
        return (
            "local_gff_protein_strict",
            "primary_local_sequence",
            "scoreable_strict",
            "Assembly GFF protein_id plus reciprocal same-gene sequence validation.",
        )
    if cds_decision == "cds_translation_strict_rescue":
        return (
            "local_cds_translation_strict",
            "primary_local_sequence",
            "scoreable_strict",
            "Assembly CDS translation plus reciprocal same-gene sequence validation.",
        )
    if uniprot_decision == "uniprot_full_length_strict_rescue":
        return (
            "external_uniprot_strict",
            "external_sensitivity_sequence",
            "scoreable_sensitivity",
            "Target-species UniProt sequence plus reciprocal same-gene validation.",
        )
    if partial_class:
        if "ambiguous" in partial_class:
            return (
                partial_class,
                "family_or_fragment_ambiguous",
                "not_scoreable_not_absence",
                "Short or family-ambiguous fragment; keep out of strict score and do not count as absence.",
            )
        return (
            partial_class,
            "partial_fragment",
            "not_scoreable_not_absence",
            "Partial same-gene evidence; keep out of strict score and do not count as absence.",
        )
    if gff_decision == "gff_sequence_forward_supported_manual_review":
        return (
            "local_gff_forward_supported_manual_review",
            "manual_review_sequence",
            "not_scoreable_review",
            "Forward same-gene sequence evidence but reciprocal support is incomplete.",
        )
    if gff_decision == "gff_sequence_weak_same_gene_not_strict":
        return (
            "local_gff_weak_same_gene_not_strict",
            "weak_sequence",
            "not_scoreable_not_absence",
            "Same-gene sequence evidence exists but coverage/reciprocal thresholds are insufficient.",
        )
    if gff_decision == "reject_gff_sequence_no_same_gene_reference":
        return (
            "local_gff_rejected_no_same_gene_reference",
            "rejected_sequence",
            "not_scoreable_review",
            "GFF-linked candidate did not hit same-gene human reference.",
        )
    if gff_annotation == "True":
        return (
            "gff_annotation_only_pending_sequence",
            "annotation_only",
            "not_scoreable_review",
            "Assembly GFF annotation supports presence but sequence validation is missing or unavailable.",
        )
    if gff_call == "gff_probable_wrong_gene_family":
        return (
            "gff_probable_wrong_gene_family",
            "rejected_or_family_conflict",
            "not_scoreable_review",
            "GFF attributes point to a related family member rather than the target gene.",
        )
    if gff_call == "not_found_in_gff":
        return (
            "not_found_in_local_gff",
            "unresolved_missing_local_annotation",
            "not_scoreable_unknown",
            "No target gene symbol, alias, or product phrase was found in local GFF.",
        )
    return (
        "unresolved_no_sequence_evidence",
        "unresolved",
        "not_scoreable_unknown",
        "No strict local or external sequence evidence was available in Phase 3.",
    )


def plot_waterfall(counts: pd.DataFrame, output_png: pathlib.Path, output_pdf: pathlib.Path) -> None:
    order = [
        "input_priority1_rows",
        "gff_annotation_rescue",
        "local_gff_protein_strict",
        "local_cds_translation_strict",
        "external_uniprot_strict",
        "partial_or_family_not_scoreable",
        "not_found_or_unresolved",
    ]
    plot_df = counts[counts["waterfall_step"].isin(order)].copy()
    plot_df["waterfall_step"] = pd.Categorical(plot_df["waterfall_step"], categories=order, ordered=True)
    plot_df = plot_df.sort_values("waterfall_step")
    labels = [
        "Input\nrows",
        "GFF\nannotation",
        "GFF protein\nstrict",
        "CDS\nstrict",
        "UniProt\nstrict",
        "Partial/family\nnot scoreable",
        "Unresolved",
    ]
    colors = ["#4C566A", "#5E81AC", "#2E7D5B", "#3A9D8F", "#B48EAD", "#D08770", "#BF616A"]
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    bars = ax.bar(labels[: len(plot_df)], plot_df["rows"].astype(int), color=colors[: len(plot_df)], width=0.72)
    for bar, value in zip(bars, plot_df["rows"].astype(int), strict=False):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1, str(value), ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("High-priority species-gene rows")
    ax.set_title("Sequence evidence ladder")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=9)
    ax.grid(axis="y", color="#D8DEE9", linewidth=0.8, alpha=0.8)
    ax.set_axisbelow(True)
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=220)
    fig.savefig(output_pdf)
    plt.close(fig)


def plot_tier_by_gene(ladder: pd.DataFrame, output_png: pathlib.Path, output_pdf: pathlib.Path) -> None:
    tier_order = [
        "local_gff_protein_strict",
        "local_cds_translation_strict",
        "external_uniprot_strict",
        "dnmt1_partial_fragment_not_absence",
        "mbd2_short_mbd_domain_ambiguous",
        "mbd3_partial_fragment_not_absence",
        "gff_annotation_only_pending_sequence",
        "not_found_in_local_gff",
        "unresolved_no_sequence_evidence",
    ]
    tab = ladder.groupby(["human_gene_symbol", "phase3_evidence_level"], as_index=False).size()
    pivot = tab.pivot(index="human_gene_symbol", columns="phase3_evidence_level", values="size").fillna(0)
    cols = [col for col in tier_order if col in pivot.columns] + [col for col in pivot.columns if col not in tier_order]
    pivot = pivot[cols].sort_index()
    colors = [
        "#2E7D5B",
        "#3A9D8F",
        "#B48EAD",
        "#D08770",
        "#D08770",
        "#D08770",
        "#EBCB8B",
        "#BF616A",
        "#8F5F5F",
    ][: len(cols)]
    display_labels = {
        "local_gff_protein_strict": "Strict local GFF protein",
        "local_cds_translation_strict": "Strict local CDS translation",
        "external_uniprot_strict": "External UniProt sensitivity",
        "dnmt1_partial_fragment_not_absence": "DNMT1 partial fragment",
        "dnmt1_longer_local_isoform_available": "DNMT1 longer isoform, not scoreable",
        "mbd2_short_mbd_domain_ambiguous": "MBD2 family/domain ambiguous",
        "mbd3_partial_fragment_not_absence": "MBD3 partial fragment",
        "gff_annotation_only_pending_sequence": "GFF annotation only",
        "gff_probable_wrong_gene_family": "Probable wrong gene family",
        "not_found_in_local_gff": "Unresolved or not found",
        "unresolved_no_sequence_evidence": "Unresolved sequence evidence",
    }
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    bottom = None
    x = range(len(pivot.index))
    for color, col in zip(colors, cols, strict=False):
        values = pivot[col].astype(int)
        ax.bar(x, values, bottom=bottom, label=display_labels.get(col, col.replace("_", " ")), color=color)
        bottom = values if bottom is None else bottom + values
    ax.set_xticks(list(x))
    ax.set_xticklabels(pivot.index, rotation=45, ha="right")
    ax.set_ylabel("Rows")
    ax.set_title("Sequence evidence classes by gene")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D8DEE9", linewidth=0.8, alpha=0.8)
    ax.set_axisbelow(True)
    ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=220)
    fig.savefig(output_pdf)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff-input", type=pathlib.Path, required=True)
    parser.add_argument("--gff-hits", type=pathlib.Path, required=True)
    parser.add_argument("--gff-decisions", type=pathlib.Path, required=True)
    parser.add_argument("--cds-decisions", type=pathlib.Path, required=True)
    parser.add_argument("--uniprot-decisions", type=pathlib.Path, required=True)
    parser.add_argument("--partial-audit", type=pathlib.Path, required=True)
    parser.add_argument("--model-impact", type=pathlib.Path, required=True)
    parser.add_argument("--ladder-output", type=pathlib.Path, required=True)
    parser.add_argument("--counts-output", type=pathlib.Path, required=True)
    parser.add_argument("--gene-summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--waterfall-png", type=pathlib.Path, required=True)
    parser.add_argument("--waterfall-pdf", type=pathlib.Path, required=True)
    parser.add_argument("--gene-figure-png", type=pathlib.Path, required=True)
    parser.add_argument("--gene-figure-pdf", type=pathlib.Path, required=True)
    args = parser.parse_args()

    base = read(args.gff_input)
    gff_hits = read(args.gff_hits)
    gff_decisions = read(args.gff_decisions)
    cds = read(args.cds_decisions)
    uniprot = read(args.uniprot_decisions)
    partial = read(args.partial_audit)
    impact = read(args.model_impact)

    keys = ["scientific_name", "human_gene_symbol"]
    ladder = base[keys + ["clade", "flight_status", "species_taxid", "best_assembly_accession", "maintenance_module", "gene_family_risk"]].copy()
    ladder = ladder.merge(
        gff_hits[
            keys
            + [
                "gff_rescue_call",
                "gff_rescue_interpretation",
                "gff_can_count_as_annotation_rescue",
                "matched_id",
                "matched_gene",
            ]
        ],
        on=keys,
        how="left",
    )
    ladder = ladder.merge(
        gff_decisions[
            keys
            + [
                "phase3_gff_sequence_decision",
                "protein_id",
                "protein_length",
                "partial",
                "validation_calls",
                "reciprocal_best_genes",
                "max_pident",
                "max_qcovhsp",
                "max_scovhsp",
            ]
        ].rename(
            columns={
                "protein_id": "gff_protein_id",
                "protein_length": "gff_protein_length",
                "partial": "gff_partial",
                "validation_calls": "gff_validation_calls",
                "reciprocal_best_genes": "gff_reciprocal_best_genes",
                "max_pident": "gff_max_pident",
                "max_qcovhsp": "gff_max_qcovhsp",
                "max_scovhsp": "gff_max_scovhsp",
            }
        ),
        on=keys,
        how="left",
    )
    ladder = ladder.merge(
        cds[
            keys
            + [
                "phase3_cds_translation_decision",
                "cds_protein_id",
                "protein_length",
                "validation_calls",
                "max_scovhsp",
            ]
        ].rename(
            columns={
                "protein_length": "cds_protein_length",
                "validation_calls": "cds_validation_calls",
                "max_scovhsp": "cds_max_scovhsp",
            }
        ),
        on=keys,
        how="left",
    )
    ladder = ladder.merge(
        uniprot[
            keys
            + [
                "phase3_uniprot_sequence_decision",
                "fetched_accession",
                "protein_length",
                "validation_calls",
                "max_scovhsp",
            ]
        ].rename(
            columns={
                "fetched_accession": "uniprot_accession",
                "protein_length": "uniprot_protein_length",
                "validation_calls": "uniprot_validation_calls",
                "max_scovhsp": "uniprot_max_scovhsp",
            }
        ),
        on=keys,
        how="left",
    )
    ladder = ladder.merge(
        partial[
            keys
            + [
                "phase3_partial_family_resolution_class",
                "phase3_partial_family_resolution_reason",
                "recommended_next_step",
            ]
        ],
        on=keys,
        how="left",
    )
    ladder = ladder.fillna("")

    classified = ladder.apply(classify, axis=1, result_type="expand")
    ladder["phase3_evidence_level"] = classified[0]
    ladder["phase3_evidence_source_class"] = classified[1]
    ladder["phase3_scoring_status"] = classified[2]
    ladder["phase3_evidence_interpretation"] = classified[3]
    ladder["best_sequence_accession"] = ladder.apply(
        lambda row: first_nonempty(row.get("gff_protein_id"), row.get("cds_protein_id"), row.get("uniprot_accession")),
        axis=1,
    )

    strict_local = ladder["phase3_evidence_level"].isin(["local_gff_protein_strict", "local_cds_translation_strict"])
    strict_any = strict_local | (ladder["phase3_evidence_level"] == "external_uniprot_strict")
    partial_or_family = ladder["phase3_evidence_source_class"].isin(["partial_fragment", "family_or_fragment_ambiguous"])
    unresolved = ladder["phase3_evidence_source_class"].isin(["unresolved_missing_local_annotation", "unresolved"])

    counts_rows = [
        {"waterfall_step": "input_priority1_rows", "rows": len(ladder), "interpretation": "Full priority-1 species-gene rows."},
        {
            "waterfall_step": "gff_annotation_rescue",
            "rows": int((ladder["gff_can_count_as_annotation_rescue"] == "True").sum()),
            "interpretation": "Rows with assembly GFF annotation support.",
        },
        {
            "waterfall_step": "local_gff_protein_strict",
            "rows": int((ladder["phase3_evidence_level"] == "local_gff_protein_strict").sum()),
            "interpretation": "Primary local GFF protein strict sequence support.",
        },
        {
            "waterfall_step": "local_cds_translation_strict",
            "rows": int((ladder["phase3_evidence_level"] == "local_cds_translation_strict").sum()),
            "interpretation": "Additional local CDS translation strict sequence support.",
        },
        {
            "waterfall_step": "external_uniprot_strict",
            "rows": int((ladder["phase3_evidence_level"] == "external_uniprot_strict").sum()),
            "interpretation": "Additional external UniProt strict sequence support.",
        },
        {
            "waterfall_step": "strict_local_sequence_total",
            "rows": int(strict_local.sum()),
            "interpretation": "Total local assembly-derived strict rows.",
        },
        {
            "waterfall_step": "strict_any_sequence_total",
            "rows": int(strict_any.sum()),
            "interpretation": "Local strict rows plus external UniProt sensitivity strict rows.",
        },
        {
            "waterfall_step": "partial_or_family_not_scoreable",
            "rows": int(partial_or_family.sum()),
            "interpretation": "Partial or family-ambiguous rows retained as not absence but not scoreable.",
        },
        {
            "waterfall_step": "not_found_or_unresolved",
            "rows": int(unresolved.sum()),
            "interpretation": "Rows still unresolved after Phase 3 evidence integration.",
        },
    ]
    counts = pd.DataFrame(counts_rows)
    gene_summary = (
        ladder.groupby(["human_gene_symbol", "phase3_evidence_level", "phase3_scoring_status"], as_index=False)
        .agg(rows=("scientific_name", "count"), species=("scientific_name", "nunique"))
        .sort_values(["human_gene_symbol", "phase3_evidence_level"])
    )

    args.ladder_output.parent.mkdir(parents=True, exist_ok=True)
    args.counts_output.parent.mkdir(parents=True, exist_ok=True)
    args.gene_summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    ladder.to_csv(args.ladder_output, sep="\t", index=False)
    counts.to_csv(args.counts_output, sep="\t", index=False)
    gene_summary.to_csv(args.gene_summary_output, sep="\t", index=False)
    plot_waterfall(counts, args.waterfall_png, args.waterfall_pdf)
    plot_tier_by_gene(ladder, args.gene_figure_png, args.gene_figure_pdf)

    trans_rows = impact[(impact["comparison_layer"] == "all_module_pgls") & (impact["model"] == "residual_module")]
    model_line = ""
    if len(trans_rows) >= 2:
        base = trans_rows[trans_rows["dataset"] == "baseline"].head(1)
        rescued = trans_rows[trans_rows["dataset"] == "rescued"].head(1)
        if not base.empty and not rescued.empty:
            model_line = (
                f"Residual transposon model after combined rescue: rank {rescued.iloc[0]['rank']}, "
                f"estimate {float(rescued.iloc[0]['estimate']):.4g}, P {float(rescued.iloc[0]['p']):.4g} "
                f"(baseline rank {base.iloc[0]['rank']}, P {float(base.iloc[0]['p']):.4g})."
            )
    lines = [
        "# Phase 3 Evidence Ladder Synthesis Report",
        "",
        f"Priority-1 rows integrated: {len(ladder)}",
        f"Rows with GFF annotation support: {int((ladder['gff_can_count_as_annotation_rescue'] == 'True').sum())}",
        f"Local GFF protein strict rows: {int((ladder['phase3_evidence_level'] == 'local_gff_protein_strict').sum())}",
        f"Additional local CDS strict rows: {int((ladder['phase3_evidence_level'] == 'local_cds_translation_strict').sum())}",
        f"Additional external UniProt strict rows: {int((ladder['phase3_evidence_level'] == 'external_uniprot_strict').sum())}",
        f"Total local strict rows: {int(strict_local.sum())}",
        f"Total strict rows including external sensitivity: {int(strict_any.sum())}",
        f"Partial/family not-scoreable rows: {int(partial_or_family.sum())}",
        f"Unresolved rows: {int(unresolved.sum())}",
        "",
        "## Model Context",
        "",
        model_line or "Model impact row was not available.",
        "",
        "## Interpretation",
        "Phase 3 now supports a clear evidence hierarchy: primary local assembly evidence is strongest, CDS translation rescue repairs a real GFF protein_id failure mode, UniProt contributes only a small external sensitivity layer, and DNMT1/MBD2/MBD3 partial-family rows should remain out of strict scoring. This is enough for a validation-focused Results section, but not enough to claim a robust bird-specific mechanism because high-coverage and birds-only sensitivity tests remain weak.",
        "",
        "## Outputs",
        f"- evidence ladder: `{args.ladder_output}`",
        f"- waterfall counts: `{args.counts_output}`",
        f"- gene summary: `{args.gene_summary_output}`",
        f"- waterfall figure: `{args.waterfall_png}`",
        f"- gene evidence figure: `{args.gene_figure_png}`",
    ]
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.ladder_output}")


if __name__ == "__main__":
    main()
