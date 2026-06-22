"""Build targeted secondary-domain rescue FASTA records."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--candidate-proteins", type=pathlib.Path, required=True)
    parser.add_argument("--gene", default="TDRD12")
    parser.add_argument("--species", default="Agelaius phoeniceus")
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t")
    proteins = pd.read_csv(args.candidate_proteins, sep="\t")
    selected = manifest[
        (manifest["human_gene_symbol"] == args.gene) & (manifest["scientific_name"] == args.species)
    ].copy()
    selected = selected.sort_values(["protein_rank_for_gene_species", "protein_accession"]).head(1)
    selected = selected.merge(
        proteins[["human_gene_symbol", "scientific_name", "protein_accession", "sequence"]],
        on=["human_gene_symbol", "scientific_name", "protein_accession"],
        how="left",
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        for _, row in selected.iterrows():
            handle.write(f">{row['domain_batch_record_id']}\n{wrap(str(row['sequence']))}\n")

    args.report.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Phase 2 Priority 1 Domain Rescue Batch Report",
        "",
        f"Target gene: {args.gene}",
        f"Target species: {args.species}",
        f"Records written: {len(selected)}",
    ]
    if len(selected):
        row = selected.iloc[0]
        lines.append(f"Selected protein: {row['protein_accession']}")
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
