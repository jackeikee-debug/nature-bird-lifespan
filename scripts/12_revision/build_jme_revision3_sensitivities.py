#!/usr/bin/env python
"""Build JME revision-3 matching, score-weight, and SAMHD1 alignment audits."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
REPORTS = ROOT / "results" / "reports"

FOUND_STATUSES = {
    "ncbi_gene_candidate",
    "gff_rescue_candidate",
    "diamond_validated_protein_candidate",
    "week4_sequence_supported_candidate",
}
SCHEMES = {
    "primary_1.0_0.8_0.5": {"high": 1.0, "medium": 0.8, "low": 0.5},
    "unweighted_local_presence": {"high": 1.0, "medium": 1.0, "low": 1.0},
    "conservative_1.0_0.7_0.25": {"high": 1.0, "medium": 0.7, "low": 0.25},
}
CLADE_ORDER = ["Aves", "Mammalia_Chiroptera", "Mammalia_nonChiroptera", "Reptilia"]
CLADE_LABELS = {
    "Aves": "Birds",
    "Mammalia_Chiroptera": "Bats",
    "Mammalia_nonChiroptera": "Other mammals",
    "Reptilia": "Reptiles",
}
CLADE_COLORS = {
    "Aves": "#D06B4F",
    "Mammalia_Chiroptera": "#7D65A8",
    "Mammalia_nonChiroptera": "#2F6F8F",
    "Reptilia": "#4D8B63",
}


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current is not None:
                records[current].append(line)
    return {name: "".join(parts) for name, parts in records.items()}


def select_confidence(row: pd.Series) -> str:
    for column in ("week4_candidate_confidence", "final_candidate_confidence", "ortholog_confidence"):
        value = row.get(column, "")
        if isinstance(value, str) and value and value != "nan":
            return value
    return ""


def select_status(row: pd.Series) -> str:
    for column in ("week4_candidate_status", "final_candidate_status", "combined_candidate_status"):
        value = row.get(column, "")
        if isinstance(value, str) and value and value != "nan":
            return value
    return ""


def build_matched_random_audit() -> pd.DataFrame:
    summary = pd.read_csv(TABLES / "phase2_W3_full_background_matched_random_set_tests.tsv", sep="\t")
    null = pd.read_csv(TABLES / "phase2_W3_full_background_matched_random_set_null.tsv", sep="\t")
    rows = []
    for _, result in summary.iterrows():
        module = result["target_module"]
        subset = null.loc[null["target_module"].eq(module)].copy()
        canonical = subset["sampled_genes"].str.split(",").apply(lambda values: ",".join(sorted(values)))
        rows.append(
            {
                "target_module": module,
                "target_gene_count": int(result["target_gene_count"]),
                "background_gene_pool_n": 139,
                "background_definition": "final 200-gene maintenance panel excluding chromatin and transposon/repeat target modules",
                "candidate_draws": int(result["attempts"]),
                "retained_nearest_sets": len(subset),
                "random_seed": 20260612,
                "sampling_without_replacement_within_set": True,
                "unique_retained_set_compositions": int(canonical.nunique()),
                "target_mean_gene_coverage": result["target_mean_gene_coverage"],
                "matched_mean_gene_coverage": result["matched_mean_coverage"],
                "mean_absolute_coverage_difference": subset["coverage_distance"].mean(),
                "p95_absolute_coverage_difference": subset["coverage_distance"].quantile(0.95),
                "maximum_absolute_coverage_difference": subset["coverage_distance"].max(),
                "gene_length_matched": False,
                "sequence_conservation_matched": False,
                "observed_residual_slope": result["observed_estimate"],
                "empirical_p_residual_slope_greater": result["empirical_p_estimate_greater"],
                "observed_residual_correlation": result["observed_r"],
                "empirical_p_residual_correlation_greater": result["empirical_p_r_greater"],
                "interpretation_boundary": "observability-matched maintenance-gene null; not a genome-wide or phylogenetic-null test",
            }
        )
    output = pd.DataFrame(rows)
    output.to_csv(TABLES / "jme_matched_random_gene_set_audit.tsv", sep="\t", index=False)
    return output


def build_weight_sensitivity() -> pd.DataFrame:
    matrix = pd.read_csv(
        ROOT / "data/processed/ortholog_matrix_primary_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv",
        sep="\t",
        low_memory=False,
    )
    eligibility = pd.read_csv(
        ROOT / "data/processed/phase2_strict_v2_scoring_eligibility_sequence_updated.tsv", sep="\t"
    )
    traits = pd.read_csv(
        ROOT / "data/processed/maintenance_lifespan_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv",
        sep="\t",
    )
    traits = traits.loc[traits["score_variant"].eq("phase2_W3_full_background_sensitivity")].drop_duplicates(
        "scientific_name"
    )
    gene_module = eligibility.set_index("human_gene_symbol")["maintenance_module_v2"].to_dict()
    matrix["maintenance_module_v2"] = matrix["human_gene_symbol"].map(gene_module)
    if matrix["maintenance_module_v2"].isna().any():
        raise ValueError("At least one final-matrix gene lacks a maintenance_module_v2 assignment")
    if matrix["human_gene_symbol"].nunique() != 200:
        raise ValueError("Weight sensitivity requires the final 200-gene matrix")

    matrix["selected_status"] = matrix.apply(select_status, axis=1)
    matrix["selected_confidence"] = matrix.apply(select_confidence, axis=1)
    matrix["accepted_local"] = matrix["selected_status"].isin(FOUND_STATUSES) & matrix[
        "selected_confidence"
    ].isin({"high", "medium", "low"})

    metadata = traits.set_index("scientific_name")
    rows = []
    for scheme, weights in SCHEMES.items():
        weighted = matrix["selected_confidence"].map(weights).fillna(0.0)
        work = matrix.assign(row_weight=np.where(matrix["accepted_local"], weighted, 0.0))
        for (species, module), subset in work.groupby(["scientific_name", "maintenance_module_v2"], sort=True):
            meta = metadata.loc[species]
            total = subset["human_gene_symbol"].nunique()
            observed = int(subset["accepted_local"].sum())
            rows.append(
                {
                    "weight_scheme": scheme,
                    "high_weight": weights["high"],
                    "medium_weight": weights["medium"],
                    "low_weight": weights["low"],
                    "scientific_name": species,
                    "opentree_tip_label": meta["opentree_tip_label"],
                    "clade": meta["clade"],
                    "log10_body_mass_g": meta["log10_body_mass_g"],
                    "log10_max_lifespan_years": np.log10(meta["max_lifespan_years"]),
                    "pgls_model_c_mass_clade_residual": meta["pgls_model_c_mass_clade_residual"],
                    "maintenance_module": module,
                    "genes_total": total,
                    "genes_observed_local": observed,
                    "coverage_fraction": observed / total,
                    "module_score": subset["row_weight"].sum() / total,
                    "external_sensitivity_rows_scored": False,
                }
            )
    output = pd.DataFrame(rows).sort_values(["weight_scheme", "maintenance_module", "scientific_name"])
    output.to_csv(TABLES / "module_weight_sensitivity_species.tsv", sep="\t", index=False)
    return output


def sequence_metrics(sequence: str, reference: str, columns: list[int]) -> dict[str, float | int]:
    if not columns:
        return {"reference_residues": 0, "paired_residues": 0, "identical_residues": 0, "coverage": np.nan, "identity": np.nan, "identity_coverage_product": np.nan}
    paired = [column for column in columns if sequence[column] != "-"]
    identical = sum(sequence[column] == reference[column] for column in paired)
    coverage = len(paired) / len(columns)
    identity = identical / len(paired) if paired else np.nan
    return {
        "reference_residues": len(columns),
        "paired_residues": len(paired),
        "identical_residues": identical,
        "coverage": coverage,
        "identity": identity,
        "identity_coverage_product": identity * coverage if np.isfinite(identity) else np.nan,
    }


def build_samhd1_alignment_audit() -> tuple[pd.DataFrame, pd.DataFrame]:
    alignment_path = ROOT / "data/interim/protein_conservation/SAMHD1.aligned.faa"
    sequences = read_fasta(alignment_path)
    reference = sequences["REF_Homo_sapiens"]
    domain_rows = pd.read_csv(TABLES / "targeted_domain_conservation_rows.tsv", sep="\t")
    domain_rows = domain_rows.loc[
        domain_rows["human_gene_symbol"].eq("SAMHD1") & domain_rows["sequence_available"].astype(bool)
    ].copy()
    traits = pd.read_csv(
        ROOT / "data/processed/maintenance_lifespan_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv",
        sep="\t",
    )
    traits = traits.loc[traits["score_variant"].eq("phase2_W3_full_background_sensitivity")].drop_duplicates(
        "scientific_name"
    )
    trait_columns = ["scientific_name", "opentree_tip_label", "pgls_model_c_mass_clade_residual"]
    domain_rows = domain_rows.drop(columns=["opentree_tip_label"], errors="ignore").merge(
        traits[trait_columns], on="scientific_name", how="left"
    )

    reference_columns = [column for column, residue in enumerate(reference) if residue != "-"]
    if len(reference_columns) != 626:
        raise ValueError(f"Expected 626 SAMHD1 reference residues, found {len(reference_columns)}")
    target_sequences = [sequences[record_id] for record_id in domain_rows["alignment_record_id"]]
    gap_fraction = {
        column: np.mean([sequence[column] == "-" for sequence in target_sequences]) for column in reference_columns
    }
    position_rows = []
    for residue_index, column in enumerate(reference_columns, start=1):
        domain_class = "SAM" if 42 <= residue_index <= 107 else "HD" if 164 <= residue_index <= 227 else "non_domain"
        position_rows.append(
            {
                "human_reference_residue": residue_index,
                "alignment_column_1based": column + 1,
                "reference_amino_acid": reference[column],
                "domain_class": domain_class,
                "target_sequence_gap_fraction": gap_fraction[column],
                "retain_original": True,
                "retain_max_gap_0.70": gap_fraction[column] <= 0.70,
                "retain_max_gap_0.50": gap_fraction[column] <= 0.50,
            }
        )
    position_qc = pd.DataFrame(position_rows)
    position_qc.to_csv(TABLES / "samhd1_alignment_position_qc.tsv", sep="\t", index=False)

    variants = {
        "original_reference_positions": lambda row: True,
        "trim_columns_gap_gt_0.70": lambda row: row["target_sequence_gap_fraction"] <= 0.70,
        "trim_columns_gap_gt_0.50": lambda row: row["target_sequence_gap_fraction"] <= 0.50,
    }
    species_rows = []
    for variant, keep in variants.items():
        retained = position_qc.loc[position_qc.apply(keep, axis=1)]
        class_columns = {
            name: (retained.loc[retained["domain_class"].isin(classes), "alignment_column_1based"] - 1).astype(int).tolist()
            for name, classes in {
                "domain": ["SAM", "HD"],
                "sam": ["SAM"],
                "hd": ["HD"],
                "nondomain": ["non_domain"],
                "whole": ["SAM", "HD", "non_domain"],
            }.items()
        }
        for _, row in domain_rows.iterrows():
            sequence = sequences[row["alignment_record_id"]]
            metrics = {name: sequence_metrics(sequence, reference, columns) for name, columns in class_columns.items()}
            species_rows.append(
                {
                    "alignment_variant": variant,
                    "scientific_name": row["scientific_name"],
                    "opentree_tip_label": row["opentree_tip_label"],
                    "clade": row["clade"],
                    "alignment_record_id": row["alignment_record_id"],
                    "selected_accession": row["selected_accession"],
                    "qualified_original_domain_coverage_ge_0.5": row["domain_reference_coverage"] >= 0.5,
                    "original_domain_reference_coverage": row["domain_reference_coverage"],
                    "pgls_model_c_mass_clade_residual": row["pgls_model_c_mass_clade_residual"],
                    "whole_reference_residues": metrics["whole"]["reference_residues"],
                    "whole_reference_coverage": metrics["whole"]["coverage"],
                    "sam_reference_residues": metrics["sam"]["reference_residues"],
                    "sam_reference_coverage": metrics["sam"]["coverage"],
                    "hd_reference_residues": metrics["hd"]["reference_residues"],
                    "hd_reference_coverage": metrics["hd"]["coverage"],
                    "domain_reference_residues": metrics["domain"]["reference_residues"],
                    "domain_reference_coverage": metrics["domain"]["coverage"],
                    "domain_aligned_identity": metrics["domain"]["identity"],
                    "domain_identity_coverage_product": metrics["domain"]["identity_coverage_product"],
                    "nondomain_reference_residues": metrics["nondomain"]["reference_residues"],
                    "nondomain_reference_coverage": metrics["nondomain"]["coverage"],
                    "nondomain_aligned_identity": metrics["nondomain"]["identity"],
                    "nondomain_identity_coverage_product": metrics["nondomain"]["identity_coverage_product"],
                    "domain_minus_nondomain_product": metrics["domain"]["identity_coverage_product"]
                    - metrics["nondomain"]["identity_coverage_product"],
                }
            )
    species_qc = pd.DataFrame(species_rows).sort_values(["alignment_variant", "clade", "scientific_name"])
    species_qc.to_csv(TABLES / "samhd1_alignment_species_qc.tsv", sep="\t", index=False)
    return position_qc, species_qc


def build_clade_faceted_samhd1_heatmap(species_qc: pd.DataFrame) -> None:
    original = species_qc.loc[
        species_qc["alignment_variant"].eq("original_reference_positions")
        & species_qc["qualified_original_domain_coverage_ge_0.5"]
    ].copy()
    sequences = read_fasta(ROOT / "data/interim/protein_conservation/SAMHD1.aligned.faa")
    reference = sequences["REF_Homo_sapiens"]
    reference_columns = [column for column, residue in enumerate(reference) if residue != "-"]
    heights = [max(2, int((original["clade"] == clade).sum())) for clade in CLADE_ORDER]
    fig, axes = plt.subplots(
        4,
        1,
        figsize=(13.5, 13.5),
        gridspec_kw={"height_ratios": heights, "hspace": 0.36},
        sharex=True,
    )
    cmap = ListedColormap(["#E7EAED", "#E6B8A9", "#277D78"])
    for axis, clade in zip(axes, CLADE_ORDER):
        subset = original.loc[original["clade"].eq(clade)].sort_values(
            ["domain_minus_nondomain_product", "scientific_name"]
        )
        matrix = []
        for record_id in subset["alignment_record_id"]:
            sequence = sequences[record_id]
            matrix.append(
                [0 if sequence[col] == "-" else 2 if sequence[col] == reference[col] else 1 for col in reference_columns]
            )
        matrix = np.asarray(matrix)
        axis.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=2, extent=[1, 626, len(subset), 0])
        axis.axvspan(42, 107, facecolor="#6AA49D", alpha=0.12, edgecolor="#4B837D", linewidth=0.9)
        axis.axvspan(164, 227, facecolor="#D06B4F", alpha=0.11, edgecolor="#D06B4F", linewidth=0.9)
        axis.set_yticks(np.arange(len(subset)) + 0.5, subset["scientific_name"].str.replace("_", " "), fontsize=6.8)
        axis.set_title(
            f"{CLADE_LABELS[clade]} (n = {len(subset)})",
            loc="left",
            fontsize=10.5,
            weight="bold",
            color=CLADE_COLORS[clade],
        )
        axis.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xlabel("Human SAMHD1 reference position")
    legend = [
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#277D78", markeredgecolor="none", markersize=8, label="Identical"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#E6B8A9", markeredgecolor="none", markersize=8, label="Substitution"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#E7EAED", markeredgecolor="none", markersize=8, label="Gap"),
    ]
    axes[0].legend(handles=legend, frameon=False, ncol=3, loc="upper right", bbox_to_anchor=(1, 1.34), fontsize=8.5)
    fig.suptitle(
        "Clade-faceted SAMHD1 residue identity and alignment coverage",
        x=0.12,
        y=0.995,
        ha="left",
        fontsize=15,
        weight="bold",
    )
    fig.subplots_adjust(left=0.24, right=0.98, top=0.95, bottom=0.06)
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / "jme_figureS7_samhd1_clade_heatmap.png", dpi=500, bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / "jme_figureS7_samhd1_clade_heatmap.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_report(matched: pd.DataFrame, weights: pd.DataFrame, positions: pd.DataFrame, species_qc: pd.DataFrame) -> None:
    qualified = species_qc.loc[
        species_qc["alignment_variant"].eq("original_reference_positions")
        & species_qc["qualified_original_domain_coverage_ge_0.5"]
    ]
    lines = [
        "# JME Revision 3 Sensitivity Inputs",
        "",
        "## Matched random sets",
        "",
    ]
    for _, row in matched.iterrows():
        lines.append(
            f"- {row['target_module']}: {int(row['retained_nearest_sets'])} unique sets retained from "
            f"{int(row['candidate_draws'])} draws; mean absolute coverage difference "
            f"{row['mean_absolute_coverage_difference']:.6f}; maximum {row['maximum_absolute_coverage_difference']:.6f}."
        )
    lines.extend(
        [
            "",
            "Sets were matched on gene count and mean observability only. Gene length and sequence conservation were not matching variables.",
            "",
            "## Score-weight sensitivity",
            "",
            f"- {weights['weight_scheme'].nunique()} schemes x {weights['maintenance_module'].nunique()} modules x {weights['scientific_name'].nunique()} species.",
            "- External sensitivity rows remain excluded from primary scoring in every weight scheme.",
            "",
            "## SAMHD1 alignment quality",
            "",
            f"- Human-reference positions audited: {len(positions)}.",
            f"- Qualified species retained under the original domain-coverage threshold: {qualified['scientific_name'].nunique()}.",
            f"- Reference positions removed at gap fractions >0.70: {(~positions['retain_max_gap_0.70']).sum()}.",
            f"- Reference positions removed at gap fractions >0.50: {(~positions['retain_max_gap_0.50']).sum()}.",
        ]
    )
    REPORTS.mkdir(parents=True, exist_ok=True)
    (REPORTS / "jme_revision3_sensitivity_inputs.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    TABLES.mkdir(parents=True, exist_ok=True)
    matched = build_matched_random_audit()
    weights = build_weight_sensitivity()
    positions, species_qc = build_samhd1_alignment_audit()
    build_clade_faceted_samhd1_heatmap(species_qc)
    write_report(matched, weights, positions, species_qc)
    print("Built JME revision-3 sensitivity inputs and Supplementary Figure 7.")


if __name__ == "__main__":
    main()
