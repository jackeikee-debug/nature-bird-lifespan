"""Select one batch from the full200 sensitivity species queue."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=pathlib.Path, required=True)
    parser.add_argument("--batch-id", required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    queue = pd.read_csv(args.queue, sep="\t")
    batch = queue[queue["batch_id"] == args.batch_id].copy()
    if batch.empty:
        raise SystemExit(f"No rows found for batch_id={args.batch_id}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    batch.to_csv(args.output, sep="\t", index=False)
    lines = [
        "# Full200 Sensitivity Selected Batch Report",
        "",
        f"Batch: {args.batch_id}",
        f"Rows: {len(batch)}",
        f"Genes: {batch['human_gene_symbol'].nunique()}",
        "",
        "## Genes",
        "",
    ]
    for gene, sub in batch.groupby("human_gene_symbol"):
        lines.append(f"- {gene}: {len(sub)} rows")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
