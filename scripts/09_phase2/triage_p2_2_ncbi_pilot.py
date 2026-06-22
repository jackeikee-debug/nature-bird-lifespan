"""Triage Phase 2 P2.2 NCBI pilot results into validation actions."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def triage_class(row: pd.Series) -> str:
    frac = row["candidate_fraction"]
    submodule = row["submodule_v2"]
    gene = row["human_gene_symbol"]
    if frac == 1:
        return "pilot_symbol_stable"
    if frac >= 0.8:
        return "pilot_mostly_stable"
    if frac >= 0.4:
        return "pilot_mixed_requires_crossdb"
    if frac == 0 and (
        gene.startswith("APOBEC3") or gene in {"PIWIL3", "PIWIL4", "DNMT3L", "TREX1", "ZCCHC3"}
    ):
        return "no_bird_symbol_hit_high_biological_ambiguity"
    if submodule == "somatic_retroelement_restriction" and frac == 0:
        return "no_bird_symbol_hit_possible_lineage_specificity"
    return "pilot_low_support_requires_external_validation"


def next_action(row: pd.Series) -> str:
    klass = row["triage_class"]
    if klass == "pilot_symbol_stable":
        return "include_in_batch1_sequence_or_crossdb_confirmation"
    if klass == "pilot_mostly_stable":
        return "include_but_check_missing_species_and_crossdb"
    if klass == "pilot_mixed_requires_crossdb":
        return "prioritize_OMA_OrthoDB_Ensembl_Compara_before_scoring"
    if klass == "no_bird_symbol_hit_high_biological_ambiguity":
        return "do_not_score_as_absence_until_domain_or_crossdb_validation"
    if klass == "no_bird_symbol_hit_possible_lineage_specificity":
        return "test_lineage_specific_loss_or_annotation_gap_with_external_sources"
    return "manual_review_before_strict_panel"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot-summary", type=pathlib.Path, required=True)
    parser.add_argument("--pilot-results", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--review-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    summary = pd.read_csv(args.pilot_summary, sep="\t")
    results = pd.read_csv(args.pilot_results, sep="\t")
    summary["triage_class"] = summary.apply(triage_class, axis=1)
    summary["recommended_next_action"] = summary.apply(next_action, axis=1)
    summary["strict_panel_ready"] = summary["triage_class"].isin(
        {"pilot_symbol_stable", "pilot_mostly_stable"}
    )

    detail = (
        results.groupby("human_gene_symbol", as_index=False)
        .agg(
            no_candidate_species=(
                "scientific_name",
                lambda x: ";".join(
                    results.loc[x.index][
                        results.loc[x.index, "ncbi_pilot_status"] == "no_candidate"
                    ]["scientific_name"].astype(str)
                ),
            ),
            candidate_species=(
                "scientific_name",
                lambda x: ";".join(
                    results.loc[x.index][
                        results.loc[x.index, "ncbi_pilot_status"] == "candidate_found"
                    ]["scientific_name"].astype(str)
                ),
            ),
        )
    )
    triage = summary.merge(detail, on="human_gene_symbol", how="left")
    triage = triage.sort_values(["strict_panel_ready", "candidate_fraction", "human_gene_symbol"])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    triage.to_csv(args.output, sep="\t", index=False)
    review = triage[~triage["strict_panel_ready"]].copy()
    args.review_output.parent.mkdir(parents=True, exist_ok=True)
    review.to_csv(args.review_output, sep="\t", index=False)

    class_counts = triage["triage_class"].value_counts().sort_index()
    strict_ready = int(triage["strict_panel_ready"].sum())
    lines = [
        "# Phase 2 P2.2 NCBI Pilot Triage Report",
        "",
        f"Pilot genes triaged: {len(triage)}",
        f"Strict-panel ready from pilot: {strict_ready}",
        f"Requires review/external validation: {len(triage) - strict_ready}",
        "",
        "## Triage Classes",
        "",
    ]
    for klass, count in class_counts.items():
        lines.append(f"- {klass}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "A stable NCBI symbol hit is not final orthology evidence, but it is enough to prioritize sequence or cross-database confirmation. Genes with zero bird symbol hits must not be interpreted as biological absence until OMA, OrthoDB, Ensembl Compara, or domain-level checks are complete.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
