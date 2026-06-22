"""Merge priority-1 NCBI Gene expansion rows into the species-level matrix."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-matrix", type=pathlib.Path, required=True)
    parser.add_argument("--ncbi-expansion", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    base = pd.read_csv(args.base_matrix, sep="\t")
    expansion = pd.read_csv(args.ncbi_expansion, sep="\t")
    base_cols = list(base.columns)

    existing_pairs = set(zip(base["human_gene_symbol"], base["scientific_name"]))
    rows = []
    for _, row in expansion.iterrows():
        pair = (row["human_gene_symbol"], row["scientific_name"])
        if pair in existing_pairs:
            continue
        found = row["ncbi_pilot_status"] == "candidate_found"
        confidence = row["ncbi_symbol_confidence"] if found else ""
        candidate_status = "ncbi_gene_candidate" if found else "priority1_expansion_no_ncbi_gene_candidate"
        rec = {col: "" for col in base_cols}
        rec.update(
            {
                "genome_panel_version": "primary",
                "scientific_name": row["scientific_name"],
                "clade": row["clade"],
                "flight_status": row["flight_status"],
                "maintenance_module": "transposon_suppression",
                "human_gene_symbol": row["human_gene_symbol"],
                "source_set": "phase2_priority1_sequence_domain_expansion",
                "species_taxid": row["species_taxid"],
                "best_assembly_accession": row["best_assembly_accession"],
                "genome_analysis_tier": row["genome_analysis_tier"],
                "ortholog_query_status": row["ncbi_pilot_status"],
                "ortholog_gene_id": row["ncbi_gene_id"] if found else "",
                "ortholog_gene_symbol": row["ncbi_gene_symbol"] if found else "",
                "ortholog_status": row["ncbi_pilot_reason"] if found else "no_symbol_taxid_candidate",
                "ortholog_source_database": "NCBI Gene" if found else "",
                "ortholog_source_url": row["ncbi_gene_url"] if found else "",
                "ortholog_confidence": confidence,
                "copy_number_estimate": "not_estimated",
                "notes": row["ncbi_gene_description"] if found else "no NCBI Gene symbol/taxid candidate in P2.4 expansion",
                "combined_candidate_status": candidate_status,
                "combined_candidate_source": "NCBI Gene P2.4 expansion" if found else "",
                "combined_candidate_confidence": confidence,
                "final_candidate_status": candidate_status,
                "final_candidate_source": "NCBI Gene P2.4 expansion" if found else "",
                "final_candidate_confidence": confidence,
                "week4_validation_batch_source": "phase2_priority1_expansion",
                "week4_sequence_status": "phase2_gene_level_sequence_domain_supported",
                "week4_candidate_status": candidate_status,
                "week4_candidate_confidence": confidence,
                "week4_candidate_source": "NCBI Gene P2.4 expansion" if found else "",
            }
        )
        rows.append(rec)

    expanded = pd.DataFrame(rows, columns=base_cols)
    merged = pd.concat([base, expanded], ignore_index=True)
    merged = merged.sort_values(["scientific_name", "maintenance_module", "human_gene_symbol"])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)

    found_rows = expanded[expanded["final_candidate_status"] == "ncbi_gene_candidate"]
    lines = [
        "# Phase 2 Priority-1 Expanded Ortholog Matrix Report",
        "",
        f"Base matrix rows: {len(base)}",
        f"Priority-1 rows added: {len(expanded)}",
        f"Priority-1 NCBI Gene candidate rows added: {len(found_rows)}",
        f"Merged matrix rows: {len(merged)}",
        "",
        "## Interpretation",
        "",
        "These added rows are NCBI Gene symbol/taxid candidates for the 16 priority-1 sequence/domain-supported genes. They support interim priority1-expanded scoring, but absence rows remain protected from biological absence claims.",
    ]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
