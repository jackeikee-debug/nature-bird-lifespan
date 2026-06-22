#!/usr/bin/env python3
"""Prepare the 10 reviewed human reference proteins for InterProScan/Pfam."""

from pathlib import Path

import pandas as pd


GENES = ["DNMT1", "DNMT3A", "DNMT3B", "HELLS", "MBD2", "MBD3", "MORC3", "SAMHD1", "SETDB2", "UHRF1"]
WORK = Path("data/interim/protein_conservation")
OUTPUT = WORK / "human_reference_10genes.faa"
MANIFEST = WORK / "human_reference_10genes_manifest.tsv"


def read_first_fasta(path: Path) -> tuple[str, str]:
    header = ""
    sequence = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if header:
                break
            header = line[1:]
        elif header and line:
            sequence.append(line)
    return header, "".join(sequence)


def main() -> None:
    rows = pd.read_csv("results/tables/targeted_protein_conservation_rows.tsv", sep="\t")
    records = []
    manifest = []
    for gene in GENES:
        header, sequence = read_first_fasta(WORK / f"{gene}.raw.faa")
        if header != "REF_Homo_sapiens" or not sequence:
            raise RuntimeError(f"Missing human reference at the start of {gene}.raw.faa")
        gene_rows = rows.loc[rows["human_gene_symbol"].eq(gene)]
        accession = gene_rows["human_reference_accession"].dropna().iloc[0]
        expected_length = int(gene_rows["human_reference_length"].dropna().iloc[0])
        if len(sequence) != expected_length:
            raise RuntimeError(f"Reference length mismatch for {gene}: {len(sequence)} != {expected_length}")
        record_id = f"REF_{gene}"
        records.append((record_id, sequence))
        manifest.append(
            {
                "record_id": record_id,
                "human_gene_symbol": gene,
                "human_reference_accession": accession,
                "protein_length": len(sequence),
            }
        )
    with OUTPUT.open("w", encoding="ascii", newline="\n") as handle:
        for record_id, sequence in records:
            handle.write(f">{record_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    pd.DataFrame(manifest).to_csv(MANIFEST, sep="\t", index=False)
    print(f"Wrote {OUTPUT} and {MANIFEST}")


if __name__ == "__main__":
    main()
