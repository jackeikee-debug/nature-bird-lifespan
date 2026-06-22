"""Summarize unresolved routes after full priority-1 GFF sequence rescue."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def route(row: pd.Series) -> tuple[str, str]:
    decision = row.get("phase3_gff_sequence_decision", "")
    gff_call = row.get("gff_rescue_call", "")
    interpretation = row.get("gff_rescue_interpretation", "")
    fetch = row.get("sequence_fetch_status", "")
    gene = row.get("human_gene_symbol", "")

    if row.get("can_count_as_strict_sequence_after_gff_sequence", "") == "True":
        return ("resolved_strict_sequence", "Already rescued by GFF-linked protein reciprocal validation.")
    if interpretation == "not_found_in_local_gff":
        return ("external_db_or_assembly_reannotation", "No local GFF symbol/product hit; check UniProt, OMA, OrthoDB, Ensembl Compara, or updated assembly annotation.")
    if interpretation == "gff_conflicting_family_signal" or gff_call == "gff_probable_wrong_gene_family":
        return ("reject_or_family_level_domain_review", "GFF hit points to a related family member; do not score unless domain/orthology review resolves it.")
    if fetch == "no_protein_id_in_gff_cds":
        return ("cds_or_genome_sequence_extraction", "GFF annotation exists but no protein_id was available; extract CDS/genomic sequence or query external protein databases.")
    if decision == "gff_sequence_weak_same_gene_not_strict":
        if gene in {"MBD2", "MBD3", "DNMT1", "DNMT3A", "DNMT3B"}:
            return ("family_specific_domain_review", "Same-gene signal is weak in a paralog-prone family; needs domain architecture and lineage-aware orthology.")
        return ("manual_sequence_threshold_review", "Same-gene hit exists but coverage/reciprocal criteria are weak; inspect length, partial flag, and reference isoforms.")
    if decision == "reject_gff_sequence_no_same_gene_reference":
        return ("probable_paralog_or_reference_gap", "Best hits do not include same-gene human reference; usually reject, but check reference completeness for difficult families.")
    if decision == "gff_sequence_forward_supported_manual_review":
        return ("manual_reciprocal_review", "Forward support exists but reciprocal evidence is incomplete; inspect hit tables and alternative references.")
    if decision == "gff_sequence_not_validated":
        return ("reference_or_validation_gap", "Protein was not validated; often missing human reference, no hit, or absent candidate FASTA record.")
    return ("unclassified_manual_review", "No route rule matched; inspect row manually.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff-hits", type=pathlib.Path, required=True)
    parser.add_argument("--sequence-decisions", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    gff = pd.read_csv(args.gff_hits, sep="\t", dtype=str).fillna("")
    decisions = pd.read_csv(args.sequence_decisions, sep="\t", dtype=str).fillna("")
    key_cols = ["scientific_name", "human_gene_symbol"]
    merged = gff.merge(
        decisions,
        on=key_cols,
        how="left",
        suffixes=("_gff", ""),
    )
    merged = merged.fillna("")
    # Prefer decision columns when present, otherwise keep GFF-only status.
    for col in [
        "clade",
        "flight_status",
        "species_taxid",
        "best_assembly_accession",
        "maintenance_module",
        "matched_gene",
        "gff_rescue_call",
        "gff_rescue_interpretation",
        "sequence_fetch_status",
        "phase3_gff_sequence_decision",
        "can_count_as_strict_sequence_after_gff_sequence",
    ]:
        gff_col = f"{col}_gff"
        if col not in merged.columns and gff_col in merged.columns:
            merged[col] = merged[gff_col]
        elif gff_col in merged.columns:
            merged[col] = merged[col].where(merged[col] != "", merged[gff_col])

    routes = merged.apply(route, axis=1, result_type="expand")
    merged["next_rescue_route"] = routes[0]
    merged["next_rescue_reason"] = routes[1]

    keep = [
        "scientific_name",
        "human_gene_symbol",
        "clade",
        "flight_status",
        "species_taxid",
        "best_assembly_accession",
        "maintenance_module",
        "gff_rescue_call",
        "gff_rescue_interpretation",
        "matched_gene",
        "protein_id",
        "sequence_fetch_status",
        "validation_calls",
        "top_reference_genes",
        "reciprocal_best_genes",
        "phase3_gff_sequence_decision",
        "can_count_as_strict_sequence_after_gff_sequence",
        "next_rescue_route",
        "next_rescue_reason",
    ]
    for col in keep:
        if col not in merged.columns:
            merged[col] = ""
    out = merged[keep].sort_values(["next_rescue_route", "scientific_name", "human_gene_symbol"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    summary = (
        out.groupby(["next_rescue_route"], as_index=False)
        .agg(rows=("human_gene_symbol", "count"), species=("scientific_name", "nunique"), genes=("human_gene_symbol", "nunique"))
        .sort_values(["next_rescue_route"])
    )
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary_output, sep="\t", index=False)

    unresolved = out[out["next_rescue_route"] != "resolved_strict_sequence"]
    by_route = unresolved["next_rescue_route"].value_counts().sort_index()
    lines = [
        "# Phase 3 Priority-1 Unresolved Route Audit",
        "",
        f"Rows audited: {len(out)}",
        f"Strict sequence-resolved rows: {int((out['next_rescue_route'] == 'resolved_strict_sequence').sum())}",
        f"Rows still needing review/rescue: {len(unresolved)}",
        "",
        "## Unresolved Routes",
    ]
    for route_name, count in by_route.items():
        lines.append(f"- {route_name}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "This audit separates biological non-support from missing evidence. Rows not found in local GFF should be routed to external orthology databases or updated annotations; weak/paralog-prone rows should be handled by family/domain review rather than scored as absence.",
            "",
            "## Outputs",
            f"- row audit: `{args.output}`",
            f"- summary: `{args.summary_output}`",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
