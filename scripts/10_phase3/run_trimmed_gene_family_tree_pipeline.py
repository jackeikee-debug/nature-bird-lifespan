#!/usr/bin/env python
"""Run IQ-TREE for trimmed targeted gene-family alignments."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def run_iqtree(alignment: Path, prefix: Path, bootstrap: int, threads: str) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "iqtree2.cmd",
            "-s",
            str(alignment),
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
    parser.add_argument("--suffix", default="trimmed_gappy70")
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--threads", default="AUTO")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    tree_dir = Path(args.tree_dir)
    for family in args.families:
        alignment = input_dir / f"{family}_mafft_{args.suffix}.faa"
        prefix = tree_dir / f"{family}_iqtree_{args.suffix}"
        print(f"Running trimmed IQ-TREE for {family}: {alignment}")
        run_iqtree(alignment, prefix, args.bootstrap, args.threads)


if __name__ == "__main__":
    main()
