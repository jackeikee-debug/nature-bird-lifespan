#!/usr/bin/env python
"""Run MAFFT and IQ-TREE for prepared targeted gene-family FASTA files."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def run_mafft(raw_fasta: Path, aligned_fasta: Path) -> None:
    aligned_fasta.parent.mkdir(parents=True, exist_ok=True)
    with aligned_fasta.open("w", encoding="ascii", newline="\n") as handle:
        subprocess.run(["mafft.bat", "--auto", str(raw_fasta)], check=True, stdout=handle)


def run_iqtree(aligned_fasta: Path, prefix: Path, bootstrap: int, threads: str) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "iqtree2.cmd",
            "-s",
            str(aligned_fasta),
            "-m",
            "MFP",
            "-B",
            str(bootstrap),
            "-alrt",
            str(bootstrap),
            "-T",
            threads,
            "-pre",
            str(prefix),
            "-redo",
        ],
        check=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--families", nargs="+", default=["DNMT_family", "MBD_family"])
    parser.add_argument("--input-dir", default="data/interim/phase3/gene_family_trees")
    parser.add_argument("--tree-dir", default="results/trees")
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--threads", default="AUTO")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    tree_dir = Path(args.tree_dir)
    for family in args.families:
        raw_fasta = input_dir / f"{family}_raw.faa"
        aligned_fasta = input_dir / f"{family}_mafft_aligned.faa"
        prefix = tree_dir / f"{family}_iqtree"
        print(f"Running MAFFT for {family}")
        run_mafft(raw_fasta, aligned_fasta)
        print(f"Running IQ-TREE for {family}")
        run_iqtree(aligned_fasta, prefix, args.bootstrap, args.threads)


if __name__ == "__main__":
    main()
