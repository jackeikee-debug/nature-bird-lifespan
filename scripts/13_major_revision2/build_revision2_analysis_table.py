"""Merge lifespan, module, genome quality, and study effort covariates."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lifespan", type=pathlib.Path, default=pathlib.Path("data/processed/maintenance_lifespan_phase2_W3_full_background_expanded.tsv"))
    parser.add_argument("--quality", type=pathlib.Path, default=pathlib.Path("data/processed/phase2_annotation_bias_covariates.tsv"))
    parser.add_argument("--genome-audit", type=pathlib.Path, default=pathlib.Path("data/processed/genome_availability_audit.tsv"))
    parser.add_argument("--effort", type=pathlib.Path, default=pathlib.Path("data/processed/jme_revision2_study_effort_covariates.tsv"))
    parser.add_argument("--output", type=pathlib.Path, default=pathlib.Path("data/processed/jme_revision2_analysis_table.tsv"))
    args = parser.parse_args()

    lifespan = pd.read_csv(args.lifespan, sep="\t")
    if "score_variant" in lifespan.columns:
        preferred = "phase2_W3_full_background_sensitivity"
        if preferred in set(lifespan["score_variant"]):
            lifespan = lifespan[lifespan["score_variant"] == preferred].copy()
    quality = pd.read_csv(args.quality, sep="\t")
    genome_audit = pd.read_csv(args.genome_audit, sep="\t", low_memory=False)
    effort = pd.read_csv(args.effort, sep="\t")
    merged = lifespan.merge(quality, on=["scientific_name", "clade", "flight_status", "genome_analysis_tier"], how="left")
    lineage = genome_audit[["scientific_name", "busco_lineage"]].drop_duplicates("scientific_name")
    merged = merged.merge(lineage, on="scientific_name", how="left")
    merged = merged.merge(effort, on="scientific_name", how="left")
    if len(merged) != 68 or merged["scientific_name"].nunique() != 68:
        raise RuntimeError(f"Expected 68 unique species, observed {len(merged)} rows and {merged['scientific_name'].nunique()} names")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {args.output}: {len(merged)} species, {len(merged.columns)} columns")


if __name__ == "__main__":
    main()
