"""Phase 3 assembly/GFF rescue for priority transposon bird rows.

This pass searches locally downloaded NCBI assembly GFF files for genes that
were not rescued by NCBI Protein. GFF evidence is annotation support only; it
does not replace reciprocal protein validation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import pathlib
import re
import urllib.parse
from collections import Counter, defaultdict


DEFAULT_DECISIONS = [
    pathlib.Path("data/processed/phase3_priority1_transposon_bird_batch01_rescue_decisions.tsv"),
    pathlib.Path("data/processed/phase3_priority1_transposon_bird_batch02_rescue_decisions.tsv"),
]
DEFAULT_ELIGIBILITY = pathlib.Path("data/processed/phase2_strict_v2_scoring_eligibility_sequence_updated.tsv")
DEFAULT_ANNOTATION_DIR = pathlib.Path("data/raw/annotation_rescue")
DEFAULT_OUTPUT = pathlib.Path("data/processed/phase3_assembly_gff_rescue_hits.tsv")
DEFAULT_SUMMARY = pathlib.Path("results/tables/phase3_assembly_gff_rescue_summary.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/phase3_assembly_gff_rescue_report.md")

TARGET_DECISIONS = {
    "not_rescued_no_protein_hit",
    "protein_broad_hit_unresolved",
    "partial_same_gene_support_not_strict",
}
FEATURE_TYPES = {"gene", "mRNA", "transcript", "CDS"}
ATTRIBUTE_FIELDS = ["gene", "Name", "gene_synonym", "standard_name", "product", "description", "Dbxref"]

PRODUCT_HINTS = {
    "DNMT1": ["dna methyltransferase 1", "dnmt1"],
    "DNMT3A": ["dna methyltransferase 3a", "dna methyltransferase 3 alpha", "dnmt3a"],
    "DNMT3B": ["dna methyltransferase 3b", "dna methyltransferase 3 beta", "dnmt3b"],
    "HELLS": ["lymphoid-specific helicase", "helicase lymphoid specific", "hells"],
    "UHRF1": ["ubiquitin like with phd and ring finger 1", "ubiquitin-like with phd and ring finger 1", "uhrf1"],
    "SETDB2": ["set domain bifurcated 2", "histone-lysine n-methyltransferase setdb2", "setdb2"],
    "MBD2": ["methyl-cpg binding domain protein 2", "methyl-cpg-binding domain protein 2", "mbd2"],
    "MBD3": ["methyl-cpg binding domain protein 3", "methyl-cpg-binding domain protein 3", "mbd3"],
    "MORC3": ["morc family cw-type zinc finger 3", "morc3"],
    "SAMHD1": ["sam and hd domain-containing protein 1", "samhd1"],
}


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_attrs(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in text.split(";"):
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        attrs[key] = urllib.parse.unquote(value)
    return attrs


def norm_token(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def word_pattern(symbol: str) -> re.Pattern[str]:
    return re.compile(rf"(?<![A-Za-z0-9]){re.escape(symbol)}(?![A-Za-z0-9])", re.IGNORECASE)


def split_aliases(text: str) -> list[str]:
    if not text:
        return []
    aliases = []
    for token in re.split(r"[;,|]", text):
        token = token.strip()
        if len(token) >= 3:
            aliases.append(token)
    return aliases


def eligibility_lookup(path: pathlib.Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(path)
    return {row["human_gene_symbol"]: row for row in rows}


def find_gff(annotation_dir: pathlib.Path, accession: str) -> pathlib.Path | None:
    acc_dir = annotation_dir / accession
    if not acc_dir.exists():
        return None
    candidates = sorted(acc_dir.glob("*genomic.gff.gz"))
    if candidates:
        return candidates[0]
    candidates = sorted(acc_dir.glob("*.gff.gz"))
    return candidates[0] if candidates else None


def wrong_family_hit(target_symbol: str, text: str) -> str:
    normalized = norm_token(text)
    family = {
        "DNMT1": ["dnmt3a", "dnmt3b"],
        "DNMT3A": ["dnmt1", "dnmt3b"],
        "DNMT3B": ["dnmt1", "dnmt3a"],
        "MBD2": ["mbd3"],
        "MBD3": ["mbd2"],
        "SETDB2": ["setdb1"],
    }
    for other in family.get(target_symbol, []):
        if other in normalized:
            return other.upper()
    return ""


def classify_match(symbol: str, attrs: dict[str, str], eligibility: dict[str, str]) -> tuple[str, str, str]:
    fields = {field: attrs.get(field, "") for field in ATTRIBUTE_FIELDS if attrs.get(field)}
    haystack = "; ".join(f"{key}={value}" for key, value in fields.items())
    wrong = wrong_family_hit(symbol, haystack)

    symbol_norm = norm_token(symbol)
    exact_fields = ["gene", "Name", "standard_name"]
    for field in exact_fields:
        value = attrs.get(field, "")
        if value and norm_token(value) == symbol_norm:
            return "gff_symbol_exact", f"{field}={value}", wrong

    for field in exact_fields:
        value = attrs.get(field, "")
        if value and word_pattern(symbol).search(value):
            return "gff_symbol_like", f"{field}={value}", wrong

    for alias in split_aliases(eligibility.get("mygene_aliases", "")):
        if alias.upper() == symbol:
            continue
        if len(alias) < 4:
            continue
        pattern = word_pattern(alias)
        for field in ["gene", "Name", "gene_synonym", "standard_name", "product", "description"]:
            value = attrs.get(field, "")
            if value and pattern.search(value):
                return "gff_alias_hit", f"{field}={value};alias={alias}", wrong

    product = " ".join([attrs.get("product", ""), attrs.get("description", "")]).lower()
    for hint in PRODUCT_HINTS.get(symbol, []):
        if re.fullmatch(r"[a-z0-9]+", hint):
            if not word_pattern(hint).search(product):
                continue
        elif hint not in product:
            continue
        if hint in product or re.fullmatch(r"[a-z0-9]+", hint):
            return "gff_product_hit", f"product_hint={hint};product={attrs.get('product', '')}", wrong

    if wrong:
        return "gff_probable_wrong_gene_family", haystack, wrong
    return "not_found_in_gff", "", ""


def better_call(current: dict[str, str] | None, candidate: dict[str, str]) -> bool:
    if current is None:
        return True
    rank = {
        "gff_symbol_exact": 5,
        "gff_symbol_like": 4,
        "gff_alias_hit": 3,
        "gff_product_hit": 2,
        "gff_probable_wrong_gene_family": 1,
        "not_found_in_gff": 0,
    }
    return rank.get(candidate["gff_rescue_call"], 0) > rank.get(current["gff_rescue_call"], 0)


def search_gff(gff_path: pathlib.Path, target_rows: list[dict[str, str]], eligibility: dict[str, dict[str, str]]) -> dict[str, dict[str, str]]:
    wanted = {row["human_gene_symbol"] for row in target_rows}
    quick_terms = {}
    for symbol in wanted:
        terms = {symbol.lower()}
        terms.update(term.lower() for term in PRODUCT_HINTS.get(symbol, []))
        for alias in split_aliases(eligibility.get(symbol, {}).get("mygene_aliases", "")):
            if len(alias) >= 4:
                terms.add(alias.lower())
        quick_terms[symbol] = terms
    best: dict[str, dict[str, str]] = {}
    seen_features: set[tuple[str, str]] = set()
    with gzip.open(gff_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] not in FEATURE_TYPES:
                continue
            lower_line = parts[8].lower()
            candidate_symbols = [
                symbol for symbol, terms in quick_terms.items()
                if any(term and term in lower_line for term in terms)
            ]
            if not candidate_symbols:
                continue
            attrs = parse_attrs(parts[8])
            feature_key = (parts[2], attrs.get("ID", parts[8]))
            if feature_key in seen_features:
                continue
            seen_features.add(feature_key)
            for symbol in candidate_symbols:
                call, match_text, wrong = classify_match(symbol, attrs, eligibility.get(symbol, {}))
                if call == "not_found_in_gff":
                    continue
                candidate = {
                    "gff_rescue_call": call,
                    "matched_feature_type": parts[2],
                    "matched_seqid": parts[0],
                    "matched_start": parts[3],
                    "matched_end": parts[4],
                    "matched_strand": parts[6],
                    "matched_id": attrs.get("ID", ""),
                    "matched_parent": attrs.get("Parent", ""),
                    "matched_gene": attrs.get("gene", ""),
                    "matched_name": attrs.get("Name", ""),
                    "matched_product": attrs.get("product", ""),
                    "matched_dbxref": attrs.get("Dbxref", ""),
                    "matched_text": match_text,
                    "wrong_family_signal": wrong,
                }
                if better_call(best.get(symbol), candidate):
                    best[symbol] = candidate
    return best


def load_decision_rows(paths: list[pathlib.Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        if not path.exists():
            continue
        rows.extend(read_tsv(path))
    return rows


def as_bool(text: str) -> bool:
    return str(text).strip().lower() == "true"


def rescue_interpretation(call: str, wrong_signal: str) -> tuple[str, str, str]:
    if wrong_signal and call != "gff_symbol_exact":
        return (
            "gff_conflicting_family_signal",
            "False",
            "GFF attributes contain a related-family symbol that may indicate paralogy or wrong-gene annotation.",
        )
    if call in {"gff_symbol_exact", "gff_symbol_like", "gff_alias_hit", "gff_product_hit"}:
        return (
            "annotation_presence_supported_pending_sequence",
            "True",
            "Assembly GFF supports gene presence, but protein/domain validation is still required for strict scoring.",
        )
    if call == "gff_probable_wrong_gene_family":
        return (
            "probable_wrong_gene_or_paralog",
            "False",
            "The best GFF match points to a related gene family member rather than the requested symbol.",
        )
    return ("not_found_in_local_gff", "False", "No matching gene symbol, alias, or product phrase was found in the local assembly GFF.")


def build_report(
    report_path: pathlib.Path,
    output: pathlib.Path,
    summary_output: pathlib.Path,
    rows: list[dict[str, str]],
    input_rows: list[dict[str, str]],
) -> None:
    total = len(rows)
    supported = sum(1 for row in rows if row["gff_can_count_as_annotation_rescue"] == "True")
    strict = sum(1 for row in rows if row["can_count_as_strict_sequence_after_gff"] == "True")
    by_call = Counter(row["gff_rescue_call"] for row in rows)
    by_species_supported = Counter(row["scientific_name"] for row in rows if row["gff_can_count_as_annotation_rescue"] == "True")
    no_gff_species = sorted({row["scientific_name"] for row in rows if row["gff_file_status"] != "found"})

    lines = [
        "# Phase 3 Assembly/GFF Rescue Report",
        "",
        "## Scope",
        f"- Input decision rows considered: {len(input_rows)}",
        f"- GFF-assessed rows: {total}",
        f"- Annotation-level rescue rows: {supported}",
        f"- Strict sequence rescue rows from this pass: {strict}",
        "",
        "GFF evidence is treated as assembly annotation support only. It can improve coverage triage and nominate rows for sequence extraction, but it does not by itself upgrade a row into the strict sequence-supported score.",
        "",
        "## GFF Calls",
    ]
    lines.extend([f"- {call}: {count}" for call, count in sorted(by_call.items())] or ["- none: 0"])
    lines.extend(["", "## Annotation Rescue by Species"])
    lines.extend([f"- {species}: {count}" for species, count in sorted(by_species_supported.items())] or ["- none: 0"])
    lines.extend(["", "## Missing Local GFF"])
    lines.extend([f"- {species}" for species in no_gff_species] or ["- none"])
    lines.extend(
        [
            "",
            "## Interpretation",
            "- Rows with `annotation_presence_supported_pending_sequence` are the next candidates for protein/CDS extraction from the same assembly or for UniProt/OMA/OrthoDB/Ensembl Compara confirmation.",
            "- Rows with related-family conflict should not be rescued unless reciprocal sequence/domain checks resolve the ambiguity.",
            "- If this GFF pass finds little support, the project should prioritize external orthology databases rather than continuing blind NCBI Protein searches.",
            "",
            "## Outputs",
            f"- Row-level table: `{output}`",
            f"- Summary table: `{summary_output}`",
        ]
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--decisions", nargs="+", type=pathlib.Path, default=DEFAULT_DECISIONS)
    parser.add_argument("--eligibility", type=pathlib.Path, default=DEFAULT_ELIGIBILITY)
    parser.add_argument("--annotation-dir", type=pathlib.Path, default=DEFAULT_ANNOTATION_DIR)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary-output", type=pathlib.Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    decisions = load_decision_rows(args.decisions)
    input_rows = [
        row for row in decisions
        if row.get("phase3_rescue_decision", "") in TARGET_DECISIONS
        and not as_bool(row.get("can_count_as_strict_sequence", "False"))
    ]
    eligibility = eligibility_lookup(args.eligibility)

    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in input_rows:
        grouped[(row["scientific_name"], row["best_assembly_accession"])].append(row)

    rows: list[dict[str, str]] = []
    for (species, accession), species_rows in sorted(grouped.items()):
        gff_path = find_gff(args.annotation_dir, accession)
        gff_status = "found" if gff_path else "missing"
        hits = search_gff(gff_path, species_rows, eligibility) if gff_path else {}
        for source_row in species_rows:
            symbol = source_row["human_gene_symbol"]
            hit = hits.get(symbol, {"gff_rescue_call": "not_found_in_gff"})
            interpretation, count_as_annotation, reason = rescue_interpretation(
                hit.get("gff_rescue_call", "not_found_in_gff"),
                hit.get("wrong_family_signal", ""),
            )
            rows.append(
                {
                    "phase3_batch_id": source_row.get("phase3_batch_id", ""),
                    "scientific_name": species,
                    "clade": source_row.get("clade", ""),
                    "flight_status": source_row.get("flight_status", ""),
                    "species_taxid": source_row.get("species_taxid", ""),
                    "best_assembly_accession": accession,
                    "genome_analysis_tier": source_row.get("genome_analysis_tier", ""),
                    "maintenance_module": source_row.get("maintenance_module", ""),
                    "human_gene_symbol": symbol,
                    "gene_family_risk": source_row.get("gene_family_risk", ""),
                    "previous_phase3_rescue_decision": source_row.get("phase3_rescue_decision", ""),
                    "previous_ncbi_protein_review_status": source_row.get("ncbi_protein_review_status", ""),
                    "gff_file_status": gff_status,
                    "gff_local_path": str(gff_path or ""),
                    "gff_rescue_call": hit.get("gff_rescue_call", "not_found_in_gff"),
                    "gff_rescue_interpretation": interpretation,
                    "gff_can_count_as_annotation_rescue": count_as_annotation,
                    "can_count_as_strict_sequence_after_gff": "False",
                    "gff_rescue_reason": reason,
                    "matched_feature_type": hit.get("matched_feature_type", ""),
                    "matched_seqid": hit.get("matched_seqid", ""),
                    "matched_start": hit.get("matched_start", ""),
                    "matched_end": hit.get("matched_end", ""),
                    "matched_strand": hit.get("matched_strand", ""),
                    "matched_id": hit.get("matched_id", ""),
                    "matched_parent": hit.get("matched_parent", ""),
                    "matched_gene": hit.get("matched_gene", ""),
                    "matched_name": hit.get("matched_name", ""),
                    "matched_product": hit.get("matched_product", ""),
                    "matched_dbxref": hit.get("matched_dbxref", ""),
                    "matched_text": hit.get("matched_text", ""),
                    "wrong_family_signal": hit.get("wrong_family_signal", ""),
                }
            )

    fields = [
        "phase3_batch_id",
        "scientific_name",
        "clade",
        "flight_status",
        "species_taxid",
        "best_assembly_accession",
        "genome_analysis_tier",
        "maintenance_module",
        "human_gene_symbol",
        "gene_family_risk",
        "previous_phase3_rescue_decision",
        "previous_ncbi_protein_review_status",
        "gff_file_status",
        "gff_local_path",
        "gff_rescue_call",
        "gff_rescue_interpretation",
        "gff_can_count_as_annotation_rescue",
        "can_count_as_strict_sequence_after_gff",
        "gff_rescue_reason",
        "matched_feature_type",
        "matched_seqid",
        "matched_start",
        "matched_end",
        "matched_strand",
        "matched_id",
        "matched_parent",
        "matched_gene",
        "matched_name",
        "matched_product",
        "matched_dbxref",
        "matched_text",
        "wrong_family_signal",
    ]
    write_tsv(args.output, rows, fields)

    summary_rows = []
    species_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        species_counts[row["scientific_name"]][row["gff_rescue_interpretation"]] += 1
        species_counts[row["scientific_name"]]["total_rows"] += 1
        if row["gff_can_count_as_annotation_rescue"] == "True":
            species_counts[row["scientific_name"]]["annotation_rescue_rows"] += 1
    for species, counts in sorted(species_counts.items()):
        summary_rows.append(
            {
                "scientific_name": species,
                "total_rows": str(counts["total_rows"]),
                "annotation_rescue_rows": str(counts["annotation_rescue_rows"]),
                "not_found_in_local_gff": str(counts["not_found_in_local_gff"]),
                "annotation_presence_supported_pending_sequence": str(counts["annotation_presence_supported_pending_sequence"]),
                "probable_wrong_gene_or_paralog": str(counts["probable_wrong_gene_or_paralog"]),
                "gff_conflicting_family_signal": str(counts["gff_conflicting_family_signal"]),
            }
        )
    summary_fields = [
        "scientific_name",
        "total_rows",
        "annotation_rescue_rows",
        "not_found_in_local_gff",
        "annotation_presence_supported_pending_sequence",
        "probable_wrong_gene_or_paralog",
        "gff_conflicting_family_signal",
    ]
    write_tsv(args.summary_output, summary_rows, summary_fields)
    build_report(args.report, args.output, args.summary_output, rows, input_rows)
    print(f"Wrote {args.output}, {args.summary_output}, and {args.report}")


if __name__ == "__main__":
    main()
