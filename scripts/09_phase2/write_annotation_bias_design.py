"""Write Phase 2 annotation-bias model design before expanded PGLS."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", type=pathlib.Path, required=True)
    parser.add_argument("--orthology-audit", type=pathlib.Path, required=True)
    parser.add_argument("--covariates-output", type=pathlib.Path, required=True)
    parser.add_argument("--model-spec-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    species = pd.read_csv(args.species, sep="\t")
    audit = pd.read_csv(args.orthology_audit, sep="\t")

    covariates = species[
        [
            "scientific_name",
            "clade",
            "flight_status",
            "genome_analysis_tier",
            "genome_quality_risk",
            "has_annotation_report",
            "scaffold_n50",
            "contig_n50",
            "coverage",
            "busco_complete",
            "busco_singlecopy",
            "busco_duplicated",
            "busco_fragmented",
            "busco_missing",
        ]
    ].copy()
    covariates["has_annotation_report"] = covariates["has_annotation_report"].astype(str)
    covariates["tier_numeric"] = covariates["genome_analysis_tier"].map(
        {
            "tier1_refseq_annotated_chromosome": 1,
            "tier2_annotated": 2,
            "tier3_assembly_only": 3,
            "tier4_low_contiguity_or_unclear": 4,
            "none": 5,
        }
    )
    args.covariates_output.parent.mkdir(parents=True, exist_ok=True)
    covariates.to_csv(args.covariates_output, sep="\t", index=False)

    module_rows = []
    for module, group in audit.groupby("maintenance_module_v2"):
        n_genes = group["human_gene_symbol"].nunique()
        n_high = int(
            (
                group["orthology_feasibility_class"]
                == "new_high_priority_validation_required"
            ).sum()
        )
        module_rows.append(
            {
                "model_family": "phase2_expanded_pgls",
                "maintenance_module_v2": module,
                "primary_formula": "lifespan_residual ~ module_score + log_body_mass + clade",
                "annotation_adjusted_formula": "lifespan_residual ~ module_score + log_body_mass + clade + tier_numeric + module_missingness",
                "coverage_adjusted_formula": "lifespan_residual ~ module_score + log_body_mass + clade + module_coverage_fraction",
                "strict_sensitivity": "run_on_strict_orthology_only",
                "tier_sensitivity": "tier1_only_and_tier1_tier2_primary",
                "outlier_sensitivity": "exclude_top_lifespan_residual_outliers",
                "clade_sensitivity": "bird_only_and_leave_out_birds",
                "submodule_sensitivity": "leave_one_submodule_out",
                "n_genes": n_genes,
                "n_high_priority_orthology": n_high,
            }
        )
    model_spec = pd.DataFrame(module_rows)
    args.model_spec_output.parent.mkdir(parents=True, exist_ok=True)
    model_spec.to_csv(args.model_spec_output, sep="\t", index=False)

    tier_counts = covariates["genome_analysis_tier"].value_counts().to_dict()
    risk_counts = covariates["genome_quality_risk"].value_counts().to_dict()
    lines = [
        "# Phase 2 P2.3 Annotation-Bias Model Design",
        "",
        "## Purpose",
        "",
        "This design is frozen before expanded module association testing. Its role is to prevent the repeat/chromatin signal from being interpreted as biology if it is better explained by genome annotation quality, genome tier, or module missingness.",
        "",
        "## Species Covariates",
        "",
        f"Primary species rows: {len(covariates)}",
        f"Genome tier counts: {tier_counts}",
        f"Genome quality risk counts: {risk_counts}",
        "",
        "Required covariates:",
        "",
        "- genome_analysis_tier / tier_numeric",
        "- genome_quality_risk",
        "- has_annotation_report",
        "- scaffold_n50 and contig_n50 where available",
        "- coverage where available",
        "- BUSCO completeness/missingness where available",
        "- module_missingness and module_coverage_fraction after v2 scoring",
        "",
        "## Predeclared Model Ladder",
        "",
        "1. Primary biology model: lifespan_residual ~ module_score + log_body_mass + clade.",
        "2. Annotation-adjusted model: add tier_numeric and module_missingness.",
        "3. Coverage-adjusted model: add module_coverage_fraction.",
        "4. Strict-only model: use strict orthology-supported genes only.",
        "5. Sensitivity models: Tier1-only, bird-only, leave-out-birds, outlier removal, leave-one-submodule-out.",
        "",
        "## Pass / Caution / Fail Logic",
        "",
        "- Pass: repeat/chromatin scores remain positive and high-ranked after tier and missingness adjustment.",
        "- Caution: effect weakens but remains directionally positive and biologically specific.",
        "- Fail: effect disappears, reverses, or is matched by many unrelated modules after annotation adjustment.",
        "",
        "## Decision Gate",
        "",
        "**pass_design_frozen**",
        "",
        "This is a design gate, not a result gate. Expanded PGLS should not be interpreted unless these annotation-bias covariates are carried forward.",
    ]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.report}")


if __name__ == "__main__":
    main()
