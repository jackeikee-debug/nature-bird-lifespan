"""Merge full200 NCBI Gene expansion rows into an expanded ortholog matrix."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-matrix", type=pathlib.Path, required=True)
    parser.add_argument("--w2-expansion", type=pathlib.Path, required=True)
    parser.add_argument("--expansion-label", default="W2")
    parser.add_argument("--source-set", default="phase2_W2_crossdb_confirm_expansion")
    parser.add_argument("--sequence-status", default="phase2_W2_gene_level_crossdb_confirm")
    parser.add_argument("--no-candidate-status", default="w2_expansion_no_ncbi_gene_candidate")
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    base = pd.read_csv(args.base_matrix, sep="\t")
    expansion = pd.read_csv(args.w2_expansion, sep="\t")
    base_cols = list(base.columns)

    existing_pairs = set(zip(base["human_gene_symbol"], base["scientific_name"]))
    rows = []
    for _, row in expansion.iterrows():
        pair = (row["human_gene_symbol"], row["scientific_name"])
        if pair in existing_pairs:
            continue
        found = row["ncbi_pilot_status"] == "candidate_found"
        confidence = row["ncbi_symbol_confidence"] if found else ""
        candidate_status = "ncbi_gene_candidate" if found else args.no_candidate_status
        rec = {col: "" for col in base_cols}
        rec.update(
            {
                "genome_panel_version": "primary",
                "scientific_name": row["scientific_name"],
                "clade": row["clade"],
                "flight_status": row["flight_status"],
                "maintenance_module": row["maintenance_module_v2"],
                "human_gene_symbol": row["human_gene_symbol"],
                "source_set": args.source_set,
                "species_taxid": row["species_taxid"],
                "best_assembly_accession": row["best_assembly_accession"],
                "genome_analysis_tier": row["genome_analysis_tier"],
                "ortholog_query_status": row["ncbi_pilot_status"],
                "ortholog_gene_id": str(row["ncbi_gene_id"]).replace(".0", "") if found else "",
                "ortholog_gene_symbol": row["ncbi_gene_symbol"] if found else "",
                "ortholog_status": row["ncbi_pilot_reason"] if found else "no_symbol_taxid_candidate",
                "ortholog_source_database": "NCBI Gene" if found else "",
                "ortholog_source_url": row["ncbi_gene_url"] if found else "",
                "ortholog_confidence": confidence,
                "copy_number_estimate": "not_estimated",
                "notes": row["ncbi_gene_description"] if found else f"no NCBI Gene symbol/taxid candidate in {args.expansion_label} expansion",
                "combined_candidate_status": candidate_status,
                "combined_candidate_source": f"NCBI Gene {args.expansion_label} expansion" if found else "",
                "combined_candidate_confidence": confidence,
                "final_candidate_status": candidate_status,
                "final_candidate_source": f"NCBI Gene {args.expansion_label} expansion" if found else "",
                "final_candidate_confidence": confidence,
                "week4_validation_batch_source": args.source_set,
                "week4_sequence_status": args.sequence_status,
                "week4_candidate_status": candidate_status,
                "week4_candidate_confidence": confidence,
                "week4_candidate_source": f"NCBI Gene {args.expansion_label} expansion" if found else "",
            }
        )
        rows.append(rec)

    added = pd.DataFrame(rows, columns=base_cols)
    merged = pd.concat([base, added], ignore_index=True)
    merged = merged.sort_values(["scientific_name", "maintenance_module", "human_gene_symbol"])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)

    found_rows = added[added["final_candidate_status"] == "ncbi_gene_candidate"]
    module_counts = added.groupby("maintenance_module", as_index=False).size()
    lines = [
        f"# Phase 2 {args.expansion_label} Expanded Ortholog Matrix Report",
        "",
        f"Base matrix rows: {len(base)}",
        f"{args.expansion_label} rows added: {len(added)}",
        f"{args.expansion_label} NCBI Gene candidate rows added: {len(found_rows)}",
        f"Merged matrix rows: {len(merged)}",
        "",
        "## Added Rows by Module",
        "",
    ]
    for _, row in module_counts.iterrows():
        lines.append(f"- {row['maintenance_module']}: {row['size']}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"{args.expansion_label} rows are first-pass NCBI Gene symbol/taxid candidates. They are suitable for sensitivity scoring, but no-candidate rows remain protected from absence claims.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
