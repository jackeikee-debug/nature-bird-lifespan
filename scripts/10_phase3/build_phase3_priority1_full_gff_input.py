"""Build pseudo-decision input for full Phase 3 priority-1 GFF rescue."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    queue = pd.read_csv(args.queue, sep="\t", dtype=str).fillna("")
    priority1 = queue[queue["rescue_priority"].astype(str) == "1"].copy()
    priority1["phase3_batch_id"] = "phase3_priority1_full_transposon_bird_gff"
    priority1["phase3_batch_goal"] = "full_priority1_assembly_gff_sequence_rescue"
    priority1["phase3_rescue_decision"] = "not_rescued_no_protein_hit"
    priority1["phase3_rescue_reason"] = "Full priority-1 queue routed directly to assembly/GFF rescue."
    priority1["ncbi_protein_review_status"] = "not_run_full_queue_gff_route"
    priority1["can_count_as_rescued_for_coverage"] = "False"
    priority1["can_count_as_strict_sequence"] = "False"
    priority1["species_rank"] = priority1.groupby("scientific_name").ngroup()
    priority1["gene_rank"] = priority1.groupby("scientific_name").cumcount()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    priority1.to_csv(args.output, sep="\t", index=False)

    by_species = priority1.groupby("scientific_name")["human_gene_symbol"].nunique().sort_index()
    lines = [
        "# Phase 3 Priority-1 Full GFF Input Report",
        "",
        f"Priority-1 rows: {len(priority1)}",
        f"Species: {priority1['scientific_name'].nunique()}",
        f"Genes: {priority1['human_gene_symbol'].nunique()}",
        "",
        "## Rows by Species",
    ]
    for species, count in by_species.items():
        lines.append(f"- {species}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "This table is a routing adapter. It does not claim protein or sequence support; it marks all priority-1 unresolved rows as ready for assembly/GFF rescue.",
            "",
            "## Outputs",
            f"- input table: `{args.output}`",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
