"""Summarize InterProScan domain confirmation for Phase 2 priority-1 proteins."""

from __future__ import annotations

import argparse
import pathlib
import re

import pandas as pd


IPR_COLUMNS = [
    "protein_accession_query",
    "sequence_md5",
    "sequence_length",
    "analysis",
    "signature_accession",
    "signature_description",
    "start",
    "end",
    "score",
    "status",
    "date",
    "interpro_accession",
    "interpro_description",
    "go_terms",
    "pathway_terms",
]

DOMAIN_PATTERNS = {
    "tudor": re.compile(r"\btudor\b|PF00567|IPR002999", re.IGNORECASE),
    "helicase": re.compile(
        r"DEAD|DEXD|helicase|P-loop containing nucleoside triphosphate hydrolase|PF00270|PF00271|IPR011545|IPR014001|IPR027417",
        re.IGNORECASE,
    ),
    "zinc_knuckle": re.compile(r"zinc|knuckle|CCHC|CHHC|PF00098|IPR001878", re.IGNORECASE),
    "kh": re.compile(r"\bKH\b|K homology|PF00013|IPR004087", re.IGNORECASE),
}


def split_record_id(record_id: str) -> dict[str, str]:
    parts = str(record_id).split("|")
    out = {
        "record_gene": parts[0] if len(parts) > 0 else "",
        "record_species": parts[1] if len(parts) > 1 else "",
        "record_protein_accession": parts[2] if len(parts) > 2 else "",
        "record_rank": parts[3].replace("rank:", "") if len(parts) > 3 else "",
    }
    return out


def has_pattern(row: pd.Series, key: str) -> bool:
    text = " ".join(
        str(row.get(col, ""))
        for col in ["signature_accession", "signature_description", "interpro_accession", "interpro_description"]
    )
    return bool(DOMAIN_PATTERNS[key].search(text))


def evaluate_rule(rule: str, flags: dict[str, bool]) -> tuple[bool, str]:
    if rule == "tudor_domain_required":
        return flags["has_tudor"], "has_tudor"
    if rule == "tudor_plus_helicase_expected":
        return flags["has_tudor"] and flags["has_helicase"], "has_tudor_and_helicase"
    if rule == "dead_box_helicase_required":
        return flags["has_helicase"], "has_helicase"
    if rule == "zinc_knuckle_required":
        return flags["has_zinc_knuckle"], "has_zinc_knuckle"
    if rule == "tudor_or_kh_domain_required":
        return flags["has_tudor"] or flags["has_kh"], "has_tudor_or_kh"
    return False, "manual_rule_not_encoded"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--interproscan-tsv", type=pathlib.Path, required=True)
    parser.add_argument("--secondary-interproscan-tsv", type=pathlib.Path, default=None)
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--row-output", type=pathlib.Path, required=True)
    parser.add_argument("--gene-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t")
    ipr_tables = []
    if args.interproscan_tsv.exists() and args.interproscan_tsv.stat().st_size:
        primary = pd.read_csv(args.interproscan_tsv, sep="\t", names=IPR_COLUMNS, dtype=str)
        primary["domain_evidence_source"] = "primary_pfam"
        ipr_tables.append(primary)
    if (
        args.secondary_interproscan_tsv is not None
        and args.secondary_interproscan_tsv.exists()
        and args.secondary_interproscan_tsv.stat().st_size
    ):
        secondary = pd.read_csv(args.secondary_interproscan_tsv, sep="\t", names=IPR_COLUMNS, dtype=str)
        secondary["domain_evidence_source"] = "secondary_targeted_interproscan"
        ipr_tables.append(secondary)
    if ipr_tables:
        ipr = pd.concat(ipr_tables, ignore_index=True)
    else:
        ipr = pd.DataFrame(columns=IPR_COLUMNS + ["domain_evidence_source"])

    if not ipr.empty:
        parsed = ipr["protein_accession_query"].apply(split_record_id).apply(pd.Series)
        ipr = pd.concat([ipr, parsed], axis=1)
        for key in DOMAIN_PATTERNS:
            ipr[f"match_{key}"] = ipr.apply(has_pattern, axis=1, key=key)

    grouped_rows = []
    for _, row in manifest.iterrows():
        record_id = row["domain_batch_record_id"]
        hits = ipr[ipr["protein_accession_query"] == record_id] if not ipr.empty else ipr
        flags = {
            "has_tudor": bool(hits.get("match_tudor", pd.Series(dtype=bool)).any()),
            "has_helicase": bool(hits.get("match_helicase", pd.Series(dtype=bool)).any()),
            "has_zinc_knuckle": bool(hits.get("match_zinc_knuckle", pd.Series(dtype=bool)).any()),
            "has_kh": bool(hits.get("match_kh", pd.Series(dtype=bool)).any()),
        }
        rule_passed, rule_basis = evaluate_rule(str(row["domain_rule"]), flags)
        hit_desc = "; ".join(
            sorted(
                set(
                    hits["signature_accession"].fillna("").astype(str)
                    + ":"
                    + hits["signature_description"].fillna("").astype(str)
                )
            )
        )
        grouped_rows.append(
            {
                **row.to_dict(),
                **flags,
                "interproscan_hit_count": len(hits),
                "primary_pfam_hit_count": int((hits.get("domain_evidence_source", "") == "primary_pfam").sum())
                if len(hits)
                else 0,
                "secondary_hit_count": int(
                    (hits.get("domain_evidence_source", "") == "secondary_targeted_interproscan").sum()
                )
                if len(hits)
                else 0,
                "domain_rule_basis": rule_basis,
                "domain_rule_passed": rule_passed,
                "domain_match_summary": hit_desc[:1000],
            }
        )

    row_summary = pd.DataFrame(grouped_rows)

    gene_rows = []
    for gene, sub in row_summary.groupby("human_gene_symbol", dropna=False):
        species = sub.groupby("scientific_name")["domain_rule_passed"].any()
        proteins_with_rule = int(sub["domain_rule_passed"].sum())
        species_with_rule = int(species.sum())
        species_tested = int(species.shape[0])
        all_species_pass = species_tested > 0 and species_with_rule == species_tested
        rule = ";".join(sorted(set(sub["domain_rule"].astype(str))))

        if str(gene).startswith("TDRD"):
            recommendation = (
                "domain_supported_keep_paralog_guard"
                if all_species_pass
                else "domain_incomplete_keep_manual_review"
            )
            strict_upgrade_allowed = False
            sensitivity_upgrade_allowed = bool(all_species_pass)
            absence_claim_allowed = False
        elif all_species_pass and species_with_rule >= 3:
            recommendation = "domain_supported_manual_queue_upgrade_candidate"
            strict_upgrade_allowed = False
            sensitivity_upgrade_allowed = True
            absence_claim_allowed = False
        else:
            recommendation = "domain_incomplete_keep_manual_review"
            strict_upgrade_allowed = False
            sensitivity_upgrade_allowed = False
            absence_claim_allowed = False

        gene_rows.append(
            {
                "human_gene_symbol": gene,
                "domain_rule": rule,
                "protein_records_tested": int(len(sub)),
                "protein_records_rule_passed": proteins_with_rule,
                "species_tested": species_tested,
                "species_rule_passed": species_with_rule,
                "all_species_domain_rule_passed": all_species_pass,
                "domain_confirmation_recommendation": recommendation,
                "strict_upgrade_allowed": strict_upgrade_allowed,
                "sensitivity_upgrade_allowed": sensitivity_upgrade_allowed,
                "absence_claim_allowed": absence_claim_allowed,
            }
        )

    gene_summary = pd.DataFrame(gene_rows).sort_values("human_gene_symbol")

    args.row_output.parent.mkdir(parents=True, exist_ok=True)
    args.gene_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    row_summary.to_csv(args.row_output, sep="\t", index=False)
    gene_summary.to_csv(args.gene_output, sep="\t", index=False)

    lines = [
        "# Phase 2 Priority 1 Domain Confirmation Report",
        "",
        "## Summary",
        "",
        f"Protein records tested: {len(row_summary)}",
        f"Genes tested: {gene_summary['human_gene_symbol'].nunique()}",
        f"Genes with all species passing encoded domain rule: {int(gene_summary['all_species_domain_rule_passed'].sum())}",
        "",
        "## Gene-Level Recommendations",
        "",
    ]
    for _, row in gene_summary.iterrows():
        lines.append(
            f"- {row['human_gene_symbol']}: {row['domain_confirmation_recommendation']} "
            f"({row['species_rule_passed']}/{row['species_tested']} species pass)"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "TDRD genes remain protected from automatic strict upgrade because shared Tudor-domain architecture does not by itself resolve TDRD paralogs. Passing TDRD rows are suitable as domain-supported sensitivity evidence and as a queue for follow-up tree/HMM discrimination.",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.gene_output}")


if __name__ == "__main__":
    main()
