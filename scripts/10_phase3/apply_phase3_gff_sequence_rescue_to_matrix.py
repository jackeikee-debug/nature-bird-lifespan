"""Apply strict Phase 3 GFF sequence rescue decisions to an ortholog matrix."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


STRICT_STATUS = "gff_sequence_strict_rescue"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, required=True)
    parser.add_argument("--decisions", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t", dtype=str).fillna("")
    decisions = pd.read_csv(args.decisions, sep="\t", dtype=str).fillna("")
    strict = decisions[decisions["phase3_gff_sequence_decision"] == STRICT_STATUS].copy()

    phase3_cols = [
        "phase3_gff_sequence_decision",
        "phase3_gff_sequence_source",
        "phase3_gff_sequence_protein_id",
        "phase3_gff_sequence_validation_tools",
        "phase3_gff_sequence_validation_calls",
        "phase3_gff_sequence_top_reference_genes",
        "phase3_gff_sequence_reciprocal_best_genes",
        "phase3_gff_sequence_max_pident",
        "phase3_gff_sequence_max_qcovhsp",
        "phase3_gff_sequence_max_scovhsp",
    ]
    for col in phase3_cols:
        if col not in matrix.columns:
            matrix[col] = ""

    key_to_decision = {
        (row["scientific_name"], row["human_gene_symbol"]): row
        for _, row in strict.iterrows()
    }
    updated = []
    missing = []
    for key, row in key_to_decision.items():
        mask = (matrix["scientific_name"] == key[0]) & (matrix["human_gene_symbol"] == key[1])
        if not mask.any():
            missing.append({"scientific_name": key[0], "human_gene_symbol": key[1], "update_status": "missing_matrix_row"})
            continue
        idx = matrix.index[mask]
        matrix.loc[idx, "final_candidate_status"] = "week4_sequence_supported_candidate"
        matrix.loc[idx, "final_candidate_source"] = "Phase3 assembly GFF protein reciprocal DIAMOND/BLASTP"
        matrix.loc[idx, "final_candidate_confidence"] = "high"
        matrix.loc[idx, "week4_candidate_status"] = "week4_sequence_supported_candidate"
        matrix.loc[idx, "week4_candidate_source"] = "Phase3 assembly GFF protein reciprocal DIAMOND/BLASTP"
        matrix.loc[idx, "week4_candidate_confidence"] = "high"
        matrix.loc[idx, "ortholog_query_status"] = "candidate_found"
        matrix.loc[idx, "ortholog_status"] = "phase3_gff_sequence_supported"
        matrix.loc[idx, "ortholog_source_database"] = "NCBI Assembly GFF + NCBI Protein"
        matrix.loc[idx, "ortholog_confidence"] = "high"
        matrix.loc[idx, "phase3_gff_sequence_decision"] = row["phase3_gff_sequence_decision"]
        matrix.loc[idx, "phase3_gff_sequence_source"] = "assembly_gff_cds_protein_id_reciprocal_validation"
        matrix.loc[idx, "phase3_gff_sequence_protein_id"] = row["protein_id"]
        matrix.loc[idx, "phase3_gff_sequence_validation_tools"] = row["validation_tools"]
        matrix.loc[idx, "phase3_gff_sequence_validation_calls"] = row["validation_calls"]
        matrix.loc[idx, "phase3_gff_sequence_top_reference_genes"] = row["top_reference_genes"]
        matrix.loc[idx, "phase3_gff_sequence_reciprocal_best_genes"] = row["reciprocal_best_genes"]
        matrix.loc[idx, "phase3_gff_sequence_max_pident"] = row["max_pident"]
        matrix.loc[idx, "phase3_gff_sequence_max_qcovhsp"] = row["max_qcovhsp"]
        matrix.loc[idx, "phase3_gff_sequence_max_scovhsp"] = row["max_scovhsp"]
        updated.append({"scientific_name": key[0], "human_gene_symbol": key[1], "update_status": "updated"})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(args.output, sep="\t", index=False)

    summary = pd.DataFrame(updated + missing)
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary_output, sep="\t", index=False)

    by_species = summary[summary["update_status"] == "updated"].groupby("scientific_name")["human_gene_symbol"].nunique()
    lines = [
        "# Phase 3 GFF Sequence Rescue Matrix Overlay Report",
        "",
        f"Strict GFF sequence decisions supplied: {len(strict)}",
        f"Matrix rows updated: {len(updated)}",
        f"Strict decisions missing from matrix: {len(missing)}",
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
            "This overlay creates a Phase 3 rescued matrix without mutating the Phase 2 source matrix. Only reciprocal same-gene GFF-linked protein decisions are promoted to high-confidence sequence-supported candidates.",
            "",
            "## Outputs",
            f"- rescued matrix: `{args.output}`",
            f"- update summary: `{args.summary_output}`",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
