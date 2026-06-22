"""Audit DNMT1 partial fragments and MBD2/MBD3 family ambiguity."""

from __future__ import annotations

import argparse
import gzip
import pathlib
import re

import pandas as pd


TARGET_CLASSES = {
    "dnmt1_partial_fragment_weak_support",
    "mbd2_mbd3_short_fragment_ambiguity",
    "mbd3_partial_fragment_weak_support",
}

PRODUCT_HINTS = {
    "DNMT1": ["dnmt1", "dna methyltransferase 1", "dnmt1 methyltransferase"],
    "MBD2": ["mbd2", "methyl-cpg-binding domain protein 2", "methyl-cpg binding domain protein 2"],
    "MBD3": ["mbd3", "methyl-cpg-binding domain protein 3", "methyl-cpg binding domain protein 3"],
}


def parse_fasta_lengths(path: pathlib.Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    records = []
    header = ""
    seq_parts: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if header:
                    records.append({"header": header, "length": len("".join(seq_parts))})
                header = line[1:].strip()
                seq_parts = []
            elif line.strip():
                seq_parts.append(line.strip())
    if header:
        records.append({"header": header, "length": len("".join(seq_parts))})
    return records


def find_asset(accession: str, asset_dir: pathlib.Path, file_type: str) -> pathlib.Path | None:
    suffix = "_protein.faa.gz" if file_type == "protein" else "_cds_from_genomic.fna.gz"
    candidates = sorted((asset_dir / accession).glob(f"*{suffix}"))
    return candidates[0] if candidates else None


def header_matches_gene(header: str, gene: str) -> bool:
    lower = header.lower()
    symbol_pattern = re.compile(rf"(?<![a-z0-9]){re.escape(gene.lower())}(?![a-z0-9])")
    hints = PRODUCT_HINTS.get(gene, [gene.lower()])
    if gene in {"MBD2", "MBD3"}:
        other = "mbd3" if gene == "MBD2" else "mbd2"
        other_phrase = "domain protein 3" if gene == "MBD2" else "domain protein 2"
        other_pattern = re.compile(rf"(?<![a-z0-9]){re.escape(other)}(?![a-z0-9])")
        if other_pattern.search(lower) or other_phrase in lower:
            return False
    return any((symbol_pattern.search(lower) if hint == gene.lower() else hint in lower) for hint in hints)


def candidate_summary(accession: str, gene: str, asset_dir: pathlib.Path, file_type: str) -> dict[str, object]:
    path = find_asset(accession, asset_dir, file_type)
    records = parse_fasta_lengths(path) if path else []
    matches = [record for record in records if header_matches_gene(str(record["header"]), gene)]
    matches = sorted(matches, key=lambda record: int(record["length"]), reverse=True)
    return {
        f"{file_type}_asset_path": str(path or ""),
        f"{file_type}_matching_records": len(matches),
        f"{file_type}_best_length": int(matches[0]["length"]) if matches else "",
        f"{file_type}_best_header": str(matches[0]["header"]) if matches else "",
    }


def classify(row: pd.Series) -> tuple[str, str]:
    gene = row["human_gene_symbol"]
    length = pd.to_numeric(row.get("protein_length", ""), errors="coerce")
    scov = pd.to_numeric(row.get("max_scovhsp", ""), errors="coerce")
    reciprocal = row.get("reciprocal_best_genes", "")
    partial = row.get("partial", "")
    protein_matches = pd.to_numeric(row.get("protein_matching_records", ""), errors="coerce")
    protein_best = pd.to_numeric(row.get("protein_best_length", ""), errors="coerce")

    if gene == "DNMT1":
        if pd.notna(scov) and scov >= 50:
            return "dnmt1_borderline_manual_review", "DNMT1 has same-gene support near useful coverage; inspect alignment/domain boundaries before scoring."
        if pd.notna(protein_best) and pd.notna(length) and protein_best > length * 1.5:
            return "dnmt1_longer_local_isoform_available", "Local assembly protein FASTA contains a longer DNMT1-like record than the GFF-selected protein."
        return "dnmt1_partial_fragment_not_absence", "DNMT1 same-gene evidence exists but subject coverage is too low for strict scoring; treat as assembly fragmentation/partial annotation, not absence."

    if gene == "MBD2":
        if reciprocal == "MBD2" and pd.notna(scov) and scov >= 50:
            return "mbd2_local_sequence_strict_candidate", "MBD2 reciprocal best is same gene with adequate subject coverage; should already be strict or reviewed for threshold consistency."
        return "mbd2_short_mbd_domain_ambiguous", "MBD2 candidate is a short MBD-domain fragment; same-gene forward evidence exists but low subject coverage and MBD2/MBD3 similarity prevent strict scoring."

    if gene == "MBD3":
        if reciprocal and reciprocal != "MBD3":
            return "mbd3_reciprocal_paralog_ambiguous", "MBD3 candidate reciprocally prefers another MBD-family member; do not score strict."
        if pd.notna(scov) and scov >= 50:
            return "mbd3_borderline_manual_review", "MBD3 same-gene support is stronger than typical short fragments; inspect why strict threshold failed."
        return "mbd3_partial_fragment_not_absence", "MBD3 evidence is partial and below strict coverage thresholds; treat as presence-like but not scoreable."

    if "true" in str(partial).lower() or "5'" in str(partial):
        return "partial_fragment_not_absence", "Candidate is explicitly marked partial and should not be treated as gene absence."
    if pd.notna(protein_matches) and protein_matches == 0:
        return "no_local_family_candidate_found", "No same-family record was found in local protein FASTA."
    return "manual_review", "Manual review required."


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family-review", type=pathlib.Path, required=True)
    parser.add_argument("--gff-decisions", type=pathlib.Path, required=True)
    parser.add_argument("--uniprot-decisions", type=pathlib.Path, required=True)
    parser.add_argument("--asset-dir", type=pathlib.Path, default=pathlib.Path("data/raw/annotation_rescue"))
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    family = pd.read_csv(args.family_review, sep="\t", dtype=str).fillna("")
    gff = pd.read_csv(args.gff_decisions, sep="\t", dtype=str).fillna("")
    uniprot = pd.read_csv(args.uniprot_decisions, sep="\t", dtype=str).fillna("")
    uniprot = uniprot[
        [
            "scientific_name",
            "human_gene_symbol",
            "fetched_accession",
            "protein_length",
            "phase3_uniprot_sequence_decision",
            "validation_calls",
            "max_scovhsp",
        ]
    ].rename(
        columns={
            "fetched_accession": "uniprot_fetched_accession",
            "protein_length": "uniprot_protein_length",
            "validation_calls": "uniprot_validation_calls",
            "max_scovhsp": "uniprot_max_scovhsp",
        }
    )

    targets = family[family["family_domain_review_class"].isin(TARGET_CLASSES)].copy()
    details = targets.merge(
        gff[
            [
                "scientific_name",
                "human_gene_symbol",
                "best_assembly_accession",
                "gff_rescue_call",
                "matched_id",
                "matched_gene",
                "protein_id",
                "partial",
                "protein_length",
                "validation_calls",
                "top_reference_genes",
                "reciprocal_best_genes",
                "max_pident",
                "max_qcovhsp",
                "max_scovhsp",
                "phase3_gff_sequence_decision",
            ]
        ],
        on=["scientific_name", "human_gene_symbol"],
        how="left",
        suffixes=("", "_gff"),
    ).merge(uniprot, on=["scientific_name", "human_gene_symbol"], how="left")

    rows = []
    for _, row in details.iterrows():
        record = row.to_dict()
        accession = str(row.get("best_assembly_accession", ""))
        gene = str(row.get("human_gene_symbol", ""))
        record.update(candidate_summary(accession, gene, args.asset_dir, "protein"))
        record.update(candidate_summary(accession, gene, args.asset_dir, "cds"))
        cls, reason = classify(pd.Series(record))
        record["phase3_partial_family_resolution_class"] = cls
        record["phase3_partial_family_resolution_reason"] = reason
        record["recommended_next_step"] = {
            "dnmt1_partial_fragment_not_absence": "Do not score as absence; prioritize genomic/CDS boundary inspection or long-read/updated annotation if this species matters for coverage.",
            "dnmt1_longer_local_isoform_available": "Extract and validate the longer local protein candidate.",
            "dnmt1_borderline_manual_review": "Inspect alignment/domain boundaries before considering sensitivity-only scoring.",
            "mbd2_short_mbd_domain_ambiguous": "Keep out of strict score; require full-length MBD2 candidate or MBD-family gene tree/domain architecture.",
            "mbd2_local_sequence_strict_candidate": "Review threshold consistency and possible strict upgrade.",
            "mbd3_reciprocal_paralog_ambiguous": "Keep out of strict score; require MBD-family gene tree/domain architecture.",
            "mbd3_borderline_manual_review": "Inspect alignment/domain boundaries before considering sensitivity-only scoring.",
            "mbd3_partial_fragment_not_absence": "Do not score as absence; require longer candidate for strict scoring.",
        }.get(cls, "Manual review.")
        rows.append(record)

    output = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)

    summary = (
        output.groupby(["human_gene_symbol", "phase3_partial_family_resolution_class"], as_index=False)
        .agg(rows=("scientific_name", "count"), species=("scientific_name", "nunique"))
        .sort_values(["human_gene_symbol", "phase3_partial_family_resolution_class"])
    )
    summary.to_csv(args.summary_output, sep="\t", index=False)

    lines = [
        "# Phase 3 Partial/Family Resolution Audit Report",
        "",
        f"Rows audited: {len(output)}",
        "",
        "## Resolution Classes",
    ]
    for _, row in summary.iterrows():
        lines.append(f"- {row['human_gene_symbol']} / {row['phase3_partial_family_resolution_class']}: {row['rows']}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "This audit separates assembly fragmentation from probable family ambiguity. Rows classified as partial fragments should not be treated as absences; rows classified as MBD-family ambiguous should stay out of strict scoring until full-length candidates or gene-family/domain evidence resolve them.",
            "",
            "## Outputs",
            f"- row audit: `{args.output}`",
            f"- summary: `{args.summary_output}`",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
