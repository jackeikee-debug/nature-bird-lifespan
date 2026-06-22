"""Combine Phase 3 matrix overlay summaries for multi-layer sensitivity reports."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", type=pathlib.Path, required=True)
    parser.add_argument("--labels", nargs="+", required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()
    if len(args.inputs) != len(args.labels):
        raise ValueError("--inputs and --labels must have the same length")

    parts = []
    for path, label in zip(args.inputs, args.labels, strict=True):
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        df["overlay_layer"] = label
        parts.append(df)
    combined = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output, sep="\t", index=False)

    updated = combined[combined["update_status"] == "updated"] if not combined.empty else combined
    counts = updated.groupby("overlay_layer")["human_gene_symbol"].count() if not updated.empty else {}
    lines = [
        "# Phase 3 Combined Overlay Summary Report",
        "",
        f"Overlay summaries combined: {len(args.inputs)}",
        f"Updated rows across combined overlays: {len(updated)}",
        "",
        "## Updated Rows by Overlay Layer",
    ]
    for layer, count in getattr(counts, "items", lambda: [])():
        lines.append(f"- {layer}: {count}")
    if len(updated) == 0:
        lines.append("- none: 0")
    lines.extend(["", "## Outputs", f"- combined summary: `{args.output}`"])
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
