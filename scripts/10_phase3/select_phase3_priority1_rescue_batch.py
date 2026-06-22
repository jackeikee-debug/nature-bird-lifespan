"""Select a reproducible Phase 3 priority-1 orthology rescue batch."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=pathlib.Path, required=True)
    parser.add_argument("--low-coverage", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--batch-id", default="phase3_priority1_transposon_bird_batch01")
    parser.add_argument("--skip-species-count", type=int, default=0)
    parser.add_argument("--species-count", type=int, default=4)
    parser.add_argument("--genes-per-species", type=int, default=10)
    args = parser.parse_args()

    queue = pd.read_csv(args.queue, sep="\t")
    low = pd.read_csv(args.low_coverage, sep="\t")

    low_trans_birds = low[
        low["maintenance_module"].eq("transposon_repeat_suppression")
        & low["clade"].eq("Aves")
    ].copy()
    low_trans_birds = low_trans_birds.sort_values(
        ["coverage_fraction", "scientific_name"],
        ascending=[True, True],
    )
    selected_species = (
        low_trans_birds["scientific_name"]
        .drop_duplicates()
        .iloc[args.skip_species_count : args.skip_species_count + args.species_count]
        .tolist()
    )

    priority_gene_order = [
        "DNMT1",
        "DNMT3A",
        "DNMT3B",
        "HELLS",
        "UHRF1",
        "SETDB2",
        "MBD2",
        "MBD3",
        "MORC3",
        "SAMHD1",
        "PIWIL1",
        "PIWIL2",
        "MOV10L1",
        "TDRD9",
        "TDRD12",
    ]
    gene_rank = {gene: idx for idx, gene in enumerate(priority_gene_order)}

    batch_rows = []
    q = queue[
        queue["rescue_priority"].eq(1)
        & queue["scientific_name"].isin(selected_species)
        & queue["maintenance_module"].eq("transposon_repeat_suppression")
    ].copy()
    q["species_rank"] = q["scientific_name"].map({name: idx for idx, name in enumerate(selected_species)})
    q["gene_rank"] = q["human_gene_symbol"].map(gene_rank).fillna(999).astype(int)
    q = q.sort_values(["species_rank", "gene_rank", "human_gene_symbol"])

    for species in selected_species:
        sub = q[q["scientific_name"].eq(species)].head(args.genes_per_species)
        batch_rows.append(sub)

    batch = pd.concat(batch_rows, ignore_index=True)
    batch.insert(0, "phase3_batch_id", args.batch_id)
    batch.insert(1, "phase3_batch_goal", "NCBI_Protein_symbol_species_rescue_then_sequence_confirmation")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    batch.to_csv(args.output, sep="\t", index=False)

    title = args.batch_id.replace("phase3_priority1_transposon_bird_", "").upper()
    lines = [
        f"# Phase 3 Priority-1 Rescue {title}",
        "",
        "## Selection Rule",
        "",
        f"- Batch id: {args.batch_id}",
        f"- Skipped low-coverage species: {args.skip_species_count}",
        f"- Species count: {args.species_count}",
        f"- Genes per species: {args.genes_per_species}",
        "- Species selected from the lowest transposon/repeat coverage birds.",
        "- Genes selected from high-priority repeat/chromatin rescue genes first.",
        "",
        "## Selected Species",
        "",
    ]
    for species in selected_species:
        cov = low_trans_birds.loc[low_trans_birds["scientific_name"].eq(species), "coverage_fraction"].iloc[0]
        lines.append(f"- {species}: transposon coverage={cov:.3f}")
    lines.extend(
        [
            "",
            "## Batch Summary",
            "",
            f"Rows: {len(batch)}",
            f"Genes: {batch['human_gene_symbol'].nunique()}",
            f"Species: {batch['scientific_name'].nunique()}",
            "",
            "## Interpretation",
            "",
            "This Phase 3 rescue batch tests whether NCBI Protein symbol/species searches can recover candidate evidence for low-coverage bird species that drive the current coverage vulnerability.",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
