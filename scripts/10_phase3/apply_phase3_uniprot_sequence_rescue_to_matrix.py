"""Apply external UniProt sequence rescue decisions to a sensitivity matrix."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


STRICT_STATUS = "uniprot_full_length_strict_rescue"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, required=True)
    parser.add_argument("--gff-decisions", type=pathlib.Path, required=True)
    parser.add_argument("--uniprot-decisions", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t", dtype=str).fillna("")
    gff_decisions = pd.read_csv(args.gff_decisions, sep="\t", dtype=str).fillna("")
    uniprot_decisions = pd.read_csv(args.uniprot_decisions, sep="\t", dtype=str).fillna("")
    gff_strict_keys = set(
        zip(
            gff_decisions.loc[gff_decisions["phase3_gff_sequence_decision"] == "gff_sequence_strict_rescue", "scientific_name"],
            gff_decisions.loc[gff_decisions["phase3_gff_sequence_decision"] == "gff_sequence_strict_rescue", "human_gene_symbol"],
            strict=False,
        )
    )
    strict = uniprot_decisions[uniprot_decisions["phase3_uniprot_sequence_decision"] == STRICT_STATUS].copy()

    phase3_cols = [
        "phase3_uniprot_sequence_decision",
        "phase3_uniprot_sequence_source",
        "phase3_uniprot_sequence_accession",
        "phase3_uniprot_sequence_validation_tools",
        "phase3_uniprot_sequence_validation_calls",
        "phase3_uniprot_sequence_top_reference_genes",
        "phase3_uniprot_sequence_reciprocal_best_genes",
        "phase3_uniprot_sequence_max_pident",
        "phase3_uniprot_sequence_max_qcovhsp",
        "phase3_uniprot_sequence_max_scovhsp",
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
        matrix.loc[idx, "final_candidate_source"] = "Phase3 UniProt external full-length reciprocal DIAMOND/BLASTP"
        matrix.loc[idx, "final_candidate_confidence"] = "medium_high_external"
        matrix.loc[idx, "week4_candidate_status"] = "week4_sequence_supported_candidate"
        matrix.loc[idx, "week4_candidate_source"] = "Phase3 UniProt external full-length reciprocal DIAMOND/BLASTP"
        matrix.loc[idx, "week4_candidate_confidence"] = "medium_high_external"
        matrix.loc[idx, "ortholog_query_status"] = "candidate_found"
        matrix.loc[idx, "ortholog_status"] = "phase3_uniprot_sequence_supported"
        matrix.loc[idx, "ortholog_source_database"] = "UniProt target-species sequence"
        matrix.loc[idx, "ortholog_confidence"] = "medium_high_external"
        matrix.loc[idx, "phase3_uniprot_sequence_decision"] = row["phase3_uniprot_sequence_decision"]
        matrix.loc[idx, "phase3_uniprot_sequence_source"] = "uniprot_target_species_reciprocal_validation"
        matrix.loc[idx, "phase3_uniprot_sequence_accession"] = row["fetched_accession"]
        matrix.loc[idx, "phase3_uniprot_sequence_validation_tools"] = row["validation_tools"]
        matrix.loc[idx, "phase3_uniprot_sequence_validation_calls"] = row["validation_calls"]
        matrix.loc[idx, "phase3_uniprot_sequence_top_reference_genes"] = row["top_reference_genes"]
        matrix.loc[idx, "phase3_uniprot_sequence_reciprocal_best_genes"] = row["reciprocal_best_genes"]
        matrix.loc[idx, "phase3_uniprot_sequence_max_pident"] = row["max_pident"]
        matrix.loc[idx, "phase3_uniprot_sequence_max_qcovhsp"] = row["max_qcovhsp"]
        matrix.loc[idx, "phase3_uniprot_sequence_max_scovhsp"] = row["max_scovhsp"]
        updated.append({"scientific_name": key[0], "human_gene_symbol": key[1], "update_status": "updated"})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(args.output, sep="\t", index=False)
    summary = pd.DataFrame(updated + skipped_existing_gff + missing)
    summary.to_csv(args.summary_output, sep="\t", index=False)

    by_gene = summary[summary["update_status"] == "updated"].groupby("human_gene_symbol")["scientific_name"].nunique()
    lines = [
        "# Phase 3 UniProt Sequence Rescue Matrix Overlay Report",
        "",
        f"Strict UniProt decisions supplied: {len(strict)}",
        f"Matrix rows updated by UniProt external sensitivity rescue: {len(updated)}",
        f"Strict UniProt rows skipped because already GFF strict: {len(skipped_existing_gff)}",
        f"Strict UniProt decisions missing from matrix: {len(missing)}",
        "",
        "## Updated Rows by Gene",
    ]
    for gene, count in by_gene.sort_index().items():
        lines.append(f"- {gene}: {count}")
    if by_gene.empty:
        lines.append("- none: 0")
    lines.extend(
        [
            "",
            "## Interpretation",
            "This matrix starts from the local assembly/GFF strict rescue matrix and overlays only additional reciprocal same-gene UniProt full-length rows. The output is a sensitivity matrix and should be reported separately from the primary GFF sequence-supported matrix.",
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
