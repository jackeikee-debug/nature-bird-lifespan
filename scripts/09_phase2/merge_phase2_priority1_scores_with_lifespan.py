"""Merge Phase 2 priority1-expanded scores with lifespan/PGLS traits."""

from __future__ import annotations

import argparse
import pathlib

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores-wide", type=pathlib.Path, required=True)
    parser.add_argument("--traits", type=pathlib.Path, required=True)
    parser.add_argument("--pgls-residuals", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    scores = pd.read_csv(args.scores_wide, sep="\t")
    traits = pd.read_csv(args.traits, sep="\t")
    pgls = pd.read_csv(args.pgls_residuals, sep="\t")
    keep_traits = [
        "scientific_name",
        "opentree_tip_label",
        "body_mass_g",
        "max_lifespan_years",
        "lifespan_residual_log10",
        "lifespan_residual_ratio",
        "data_quality",
        "sample_size",
        "specimen_origin",
    ]
    merged = scores.merge(traits[keep_traits], on="scientific_name", how="left")
    merged["log10_body_mass_g"] = np.log10(pd.to_numeric(merged["body_mass_g"], errors="coerce"))
    merged = merged.merge(pgls, on="opentree_tip_label", how="left")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)

    lines = [
        "# Phase 2 Priority1-Expanded Lifespan Merge Report",
        "",
        f"Score rows: {len(scores)}",
        f"Merged rows: {len(merged)}",
        f"Species: {merged['scientific_name'].nunique()}",
        f"Variants: {merged['score_variant'].nunique()}",
        f"Rows missing lifespan residual: {int(merged['lifespan_residual_log10'].isna().sum())}",
        f"Rows missing PGLS residual: {int(merged['pgls_model_c_mass_clade_residual'].isna().sum())}",
    ]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
