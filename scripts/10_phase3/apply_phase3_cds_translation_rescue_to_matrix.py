"""Apply strict assembly CDS translation rescue decisions to a matrix."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


STRICT_STATUS = "cds_translation_strict_rescue"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, required=True)
    parser.add_argument("--gff-decisions", type=pathlib.Path, required=True)
    parser.add_argument("--cds-decisions", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t", dtype=str).fillna("")
    gff_decisions = pd.read_csv(args.gff_decisions, sep="\t", dtype=str).fillna("")
    cds_decisions = pd.read_csv(args.cds_decisions, sep="\t", dtype=str).fillna("")
    gff_strict_keys = set(
        zip(
            gff_decisions.loc[gff_decisions["phase3_gff_sequence_decision"] == "gff_sequence_strict_rescue", "scientific_name"],
            gff_decisions.loc[gff_decisions["phase3_gff_sequence_decision"] == "gff_sequence_strict_rescue", "human_gene_symbol"],
            strict=False,
        )
    )
    strict = cds_decisions[cds_decisions["phase3_cds_translation_decision"] == STRICT_STATUS].copy()

    phase3_cols = [
        "phase3_cds_translation_decision",
        "phase3_cds_translation_source",
        "phase3_cds_translation_protein_id",
        "phase3_cds_translation_validation_tools",
        "phase3_cds_translation_validation_calls",
        "phase3_cds_translation_top_reference_genes",
        "phase3_cds_translation_reciprocal_best_genes",
        "phase3_cds_translation_max_pident",
        "phase3_cds_translation_max_qcovhsp",
        "phase3_cds_translation_max_scovhsp",
    ]
    for col in phase3_cols:
        if col not in matrix.columns:
            matrix[col] = ""

    updated = []
    skipped_existing_gff = []
    missing = []
    for _, row in strict.iterrows():
        key = (row["scientific_name"], row["human_gene_symbol"])
        if key in gff_strict_keys:
            skipped_existing_gff.append({"scientific_name": key[0], "human_gene_symbol": key[1], "update_status": "already_gff_strict"})
            continue
        mask = (matrix["scientific_name"] == key[0]) & (matrix["human_gene_symbol"] == key[1])
        if not mask.any():
            missing.append({"scientific_name": key[0], "human_gene_symbol": key[1], "update_status": "missing_matrix_row"})
            continue
        idx = matrix.index[mask]
        matrix.loc[idx, "final_candidate_status"] = "week4_sequence_supported_candidate"
        matrix.loc[idx, "final_candidate_source"] = "Phase3 assembly CDS translation reciprocal DIAMOND/BLASTP"
        matrix.loc[idx, "final_candidate_confidence"] = "high"
        matrix.loc[idx, "week4_candidate_status"] = "week4_sequence_supported_candidate"
        matrix.loc[idx, "week4_candidate_source"] = "Phase3 assembly CDS translation reciprocal DIAMOND/BLASTP"
        matrix.loc[idx, "week4_candidate_confidence"] = "high"
        matrix.loc[idx, "ortholog_query_status"] = "candidate_found"
        matrix.loc[idx, "ortholog_status"] = "phase3_cds_translation_sequence_supported"
        matrix.loc[idx, "ortholog_source_database"] = "NCBI Assembly CDS translation"
        matrix.loc[idx, "ortholog_confidence"] = "high"
        matrix.loc[idx, "phase3_cds_translation_decision"] = row["phase3_cds_translation_decision"]
        matrix.loc[idx, "phase3_cds_translation_source"] = "assembly_cds_translation_reciprocal_validation"
        matrix.loc[idx, "phase3_cds_translation_protein_id"] = row["cds_protein_id"]
        matrix.loc[idx, "phase3_cds_translation_validation_tools"] = row["validation_tools"]
        matrix.loc[idx, "phase3_cds_translation_validation_calls"] = row["validation_calls"]
        matrix.loc[idx, "phase3_cds_translation_top_reference_genes"] = row["top_reference_genes"]
        matrix.loc[idx, "phase3_cds_translation_reciprocal_best_genes"] = row["reciprocal_best_genes"]
        matrix.loc[idx, "phase3_cds_translation_max_pident"] = row["max_pident"]
        matrix.loc[idx, "phase3_cds_translation_max_qcovhsp"] = row["max_qcovhsp"]
        matrix.loc[idx, "phase3_cds_translation_max_scovhsp"] = row["max_scovhsp"]
        updated.append({"scientific_name": key[0], "human_gene_symbol": key[1], "update_status": "updated"})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(args.output, sep="\t", index=False)
    summary = pd.DataFrame(updated + skipped_existing_gff + missing)
    summary.to_csv(args.summary_output, sep="\t", index=False)

    by_species = summary[summary["update_status"] == "updated"].groupby("scientific_name")["human_gene_symbol"].nunique()
    lines = [
        "# Phase 3 CDS Translation Rescue Matrix Overlay Report",
        "",
        f"Strict CDS translation decisions supplied: {len(strict)}",
        f"Matrix rows updated by local CDS translation rescue: {len(updated)}",
        f"Strict CDS rows skipped because already GFF strict: {len(skipped_existing_gff)}",
        f"Strict CDS decisions missing from matrix: {len(missing)}",
        "",
        "## Updated Rows by Species",
    ]
    for species, count in by_species.sort_index().items():
        lines.append(f"- {species}: {count}")
    if by_species.empty:
        lines.append("- none: 0")
    lines.extend(
        [
            "",
            "## Interpretation",
            "This matrix starts from the local assembly/GFF strict rescue matrix and overlays only additional reciprocal same-gene assembly CDS translation rows. It remains a separate Phase 3 sensitivity matrix, but the evidence source is local assembly-derived rather than external-only.",
            "",
            "## Outputs",
            f"- rescued matrix: `{args.output}`",
            f"- update summary: `{args.summary_output}`",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
