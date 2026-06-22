#!/usr/bin/env python
"""Trim gappy columns from targeted family alignments for sensitivity trees."""

from __future__ import annotations

import argparse
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = None
    seq_parts: list[str] = []
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts)))
            header = line[1:]
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def wrap(seq: str, width: int = 70) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for header, seq in records:
            handle.write(f">{header}\n{wrap(seq)}\n")


def trim_alignment(records: list[tuple[str, str]], max_col_gap: float, min_non_gap_aa: int):
    if not records:
        raise ValueError("empty alignment")
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise ValueError(f"alignment sequences have nonuniform lengths: {sorted(lengths)}")
    n_seq = len(records)
    aln_len = len(records[0][1])
    keep_cols = []
    for idx in range(aln_len):
        col = [seq[idx] for _, seq in records]
        gap_fraction = sum(aa in "-?." for aa in col) / n_seq
        if gap_fraction <= max_col_gap:
            keep_cols.append(idx)
    trimmed = []
    dropped_low_info = 0
    for header, seq in records:
        new_seq = "".join(seq[i] for i in keep_cols)
        non_gap = sum(aa not in "-?." for aa in new_seq)
        if non_gap < min_non_gap_aa:
            dropped_low_info += 1
            continue
        trimmed.append((header, new_seq))
    return trimmed, keep_cols, dropped_low_info


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--families", nargs="+", default=["DNMT_family", "MBD_family"])
    parser.add_argument("--input-dir", default="data/interim/phase3/gene_family_trees")
    parser.add_argument("--max-col-gap", type=float, default=0.70)
    parser.add_argument("--min-non-gap-aa", type=int, default=80)
    parser.add_argument("--summary-output", default="results/tables/gene_family_trimmed_alignment_qc.tsv")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    rows = []
    for family in args.families:
        in_path = input_dir / f"{family}_mafft_aligned.faa"
        out_path = input_dir / f"{family}_mafft_trimmed_gappy70.faa"
        records = read_fasta(in_path)
        trimmed, keep_cols, dropped = trim_alignment(records, args.max_col_gap, args.min_non_gap_aa)
        write_fasta(trimmed, out_path)
        rows.append(
            {
                "family": family,
                "input_alignment": str(in_path).replace("\\", "/"),
                "trimmed_alignment": str(out_path).replace("\\", "/"),
                "input_sequences": len(records),
                "trimmed_sequences": len(trimmed),
                "input_columns": len(records[0][1]) if records else 0,
                "trimmed_columns": len(keep_cols),
                "columns_retained_fraction": len(keep_cols) / len(records[0][1]) if records else 0,
                "dropped_low_information_sequences": dropped,
                "max_column_gap_fraction": args.max_col_gap,
                "min_non_gap_aa": args.min_non_gap_aa,
            }
        )
        print(f"{family}: kept {len(trimmed)}/{len(records)} sequences and {len(keep_cols)}/{len(records[0][1])} columns")

    out = Path(args.summary_output)
    out.parent.mkdir(parents=True, exist_ok=True)
    headers = list(rows[0].keys())
    with out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(headers) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[h]) for h in headers) + "\n")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
