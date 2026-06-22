"""Run a UniProt-backed cross-database confirmation pilot for v2 genes."""

from __future__ import annotations

import argparse
import csv
import io
import json
import pathlib
import re
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd


UNIPROT_FIELDS = ",".join(
    [
        "accession",
        "reviewed",
        "gene_names",
        "organism_name",
        "protein_name",
        "length",
        "xref_refseq",
        "xref_orthodb",
        "xref_oma",
    ]
)


def symbol_tokens(gene_names: str) -> set[str]:
    tokens = set()
    for token in re.split(r"[\s,;]+", str(gene_names).upper()):
        clean = re.sub(r"[^A-Z0-9_.-]", "", token)
        if clean:
            tokens.add(clean)
    return tokens


def cache_name(gene: str, taxid: str) -> str:
    safe_gene = re.sub(r"[^A-Za-z0-9_.-]", "_", gene)
    safe_taxid = re.sub(r"[^0-9]", "_", str(taxid))
    return f"{safe_gene}__taxid_{safe_taxid}.tsv"


def fetch_uniprot(gene: str, taxid: str, cache_dir: pathlib.Path, sleep: float) -> str:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / cache_name(gene, taxid)
    if cache_path.exists():
        return cache_path.read_text(encoding="utf-8")

    query = f"(gene:{gene}) AND (taxonomy_id:{taxid})"
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {
            "query": query,
            "fields": UNIPROT_FIELDS,
            "format": "tsv",
            "size": "10",
        }
    )
    last_error = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                text = response.read().decode("utf-8", errors="replace")
            cache_path.write_text(text, encoding="utf-8")
            if sleep:
                time.sleep(sleep)
            return text
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            time.sleep(2 * (attempt + 1))
    error_text = f"ERROR\t{type(last_error).__name__}: {last_error}\n"
    cache_path.write_text(error_text, encoding="utf-8")
    return error_text


def parse_uniprot_tsv(text: str) -> list[dict[str, str]]:
    if text.startswith("ERROR\t"):
        return []
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    return [dict(row) for row in reader]


def summarize_hits(gene: str, hits: list[dict[str, str]]) -> dict[str, object]:
    gene_upper = gene.upper()
    symbol_like = []
    broad = []
    for hit in hits:
        if gene_upper in symbol_tokens(hit.get("Gene Names", "")):
            symbol_like.append(hit)
        else:
            broad.append(hit)

    selected = symbol_like or broad
    top = selected[:3]
    oma_xrefs = [hit.get("OMA", "") for hit in top if str(hit.get("OMA", "")).strip()]
    orthodb_xrefs = [hit.get("OrthoDB", "") for hit in top if str(hit.get("OrthoDB", "")).strip()]
    refseq_xrefs = [hit.get("RefSeq", "") for hit in top if str(hit.get("RefSeq", "")).strip()]
    reviewed = [hit.get("Reviewed", "") for hit in top if str(hit.get("Reviewed", "")).strip()]

    if symbol_like:
        decision = "uniprot_symbol_like_hit"
    elif hits:
        decision = "uniprot_broad_hit_no_exact_symbol"
    else:
        decision = "no_uniprot_hit"

    return {
        "uniprot_hit_count": len(hits),
        "uniprot_symbol_like_hit_count": len(symbol_like),
        "uniprot_broad_hit_count": len(broad),
        "uniprot_decision": decision,
        "top_uniprot_accessions": ";".join(hit.get("Entry", "") for hit in top),
        "top_uniprot_reviewed": ";".join(reviewed),
        "top_uniprot_gene_names": " | ".join(hit.get("Gene Names", "") for hit in top),
        "top_uniprot_protein_names": " | ".join(hit.get("Protein names", "") for hit in top),
        "top_uniprot_lengths": ";".join(str(hit.get("Length", "")) for hit in top),
        "uniprot_refseq_xrefs": " | ".join(refseq_xrefs),
        "uniprot_orthodb_xrefs": " | ".join(orthodb_xrefs),
        "uniprot_oma_xrefs": " | ".join(oma_xrefs),
        "uniprot_oma_xref_supported": bool(oma_xrefs),
        "uniprot_orthodb_xref_supported": bool(orthodb_xrefs),
    }


def row_confirmation_status(row: pd.Series) -> str:
    if row["uniprot_decision"] == "uniprot_symbol_like_hit":
        if row.get("gene_family_risk", "") == "high_paralog_family" or "domain_check" in str(row["crossdb_route"]):
            return "protein_supported_sequence_or_domain_check_required"
        if row["uniprot_oma_xref_supported"] or row["uniprot_orthodb_xref_supported"]:
            return "protein_and_crossdb_xref_supported"
        return "protein_supported_crossdb_xref_missing"
    if row["uniprot_decision"] == "uniprot_broad_hit_no_exact_symbol":
        return "ambiguous_broad_uniprot_hit"
    return "no_uniprot_support_not_absence"


def gene_upgrade_decision(row: pd.Series) -> str:
    symbol_hits = int(row["species_with_symbol_like_uniprot"])
    rows = int(row["species_rows"])
    high_paralog = bool(row["any_high_paralog_family"])
    route = str(row["crossdb_route_examples"])
    if symbol_hits >= max(4, rows - 1):
        if high_paralog or "domain_check" in route:
            return "sequence_or_domain_confirmation_required_before_strict_upgrade"
        return "upgrade_candidate_pending_matrix_merge"
    if symbol_hits >= 3:
        return "mixed_support_expand_species_or_sequence_check"
    if symbol_hits > 0:
        return "weak_support_manual_review"
    return "hold_no_uniprot_support_not_absence"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-queue", type=pathlib.Path, required=True)
    parser.add_argument("--priority", type=int, default=1)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--gene-summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--sequence-queue-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--cache-dir", type=pathlib.Path, required=True)
    parser.add_argument("--sleep", type=float, default=0.1)
    args = parser.parse_args()

    queue = pd.read_csv(args.species_queue, sep="\t")
    if "gene_family_risk" not in queue.columns:
        queue["gene_family_risk"] = "not_available_in_species_queue"
    pilot = queue[queue["crossdb_priority"].astype(int) <= args.priority].copy()

    rows = []
    for _, row in pilot.iterrows():
        text = fetch_uniprot(
            str(row["human_gene_symbol"]),
            str(row["species_taxid"]),
            args.cache_dir,
            args.sleep,
        )
        hits = parse_uniprot_tsv(text)
        summary = summarize_hits(str(row["human_gene_symbol"]), hits)
        record = row.to_dict()
        record.update(summary)
        rows.append(record)

    result = pd.DataFrame(rows)
    result["confirmation_status"] = result.apply(row_confirmation_status, axis=1)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)

    gene_summary = (
        result.groupby(
            [
                "crossdb_priority",
                "human_gene_symbol",
                "maintenance_module_v2",
                "submodule_v2",
            ],
            as_index=False,
        )
        .agg(
            species_rows=("scientific_name", "nunique"),
            species_with_symbol_like_uniprot=("uniprot_decision", lambda s: int((s == "uniprot_symbol_like_hit").sum())),
            species_with_oma_xref=("uniprot_oma_xref_supported", "sum"),
            species_with_orthodb_xref=("uniprot_orthodb_xref_supported", "sum"),
            species_without_uniprot_support=("uniprot_decision", lambda s: int((s == "no_uniprot_hit").sum())),
            any_high_paralog_family=("gene_family_risk", lambda s: bool((s == "high_paralog_family").any())),
            crossdb_route_examples=("crossdb_route", lambda s: ";".join(sorted(set(map(str, s))))),
            confirmation_status_examples=("confirmation_status", lambda s: ";".join(sorted(set(map(str, s))))),
            top_uniprot_accessions=("top_uniprot_accessions", lambda s: " | ".join(x for x in s if str(x).strip())),
        )
    )
    gene_summary["uniprot_symbol_like_fraction"] = (
        gene_summary["species_with_symbol_like_uniprot"] / gene_summary["species_rows"]
    ).round(3)
    gene_summary["strict_upgrade_decision"] = gene_summary.apply(gene_upgrade_decision, axis=1)
    args.gene_summary_output.parent.mkdir(parents=True, exist_ok=True)
    gene_summary.to_csv(args.gene_summary_output, sep="\t", index=False)

    sequence_queue = result[
        result["confirmation_status"].isin(
            {
                "protein_supported_sequence_or_domain_check_required",
                "ambiguous_broad_uniprot_hit",
                "no_uniprot_support_not_absence",
            }
        )
    ].copy()
    args.sequence_queue_output.parent.mkdir(parents=True, exist_ok=True)
    sequence_queue.to_csv(args.sequence_queue_output, sep="\t", index=False)

    status_counts = result["confirmation_status"].value_counts().sort_index()
    decision_counts = gene_summary["strict_upgrade_decision"].value_counts().sort_index()
    lines = [
        "# Phase 2 Cross-Database UniProt Pilot Report",
        "",
        "## Summary",
        "",
        f"Priority threshold: {args.priority}",
        f"Gene-species rows queried: {len(result)}",
        f"Genes queried: {result['human_gene_symbol'].nunique()}",
        f"Species queried: {result['scientific_name'].nunique()}",
        "",
        "## Row-Level Confirmation Status",
        "",
    ]
    for status, count in status_counts.items():
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Gene-Level Upgrade Decisions", ""])
    for decision, count in decision_counts.items():
        lines.append(f"- {decision}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "A UniProt symbol-like hit with OMA or OrthoDB cross-reference is strong cross-database support, but high-paralog or domain-check routes are still held for domain or sequence confirmation before strict v2 upgrade. No UniProt hit is not treated as biological absence.",
            "",
            "## Outputs",
            "",
            f"- row results: {args.output.as_posix()}",
            f"- gene summary: {args.gene_summary_output.as_posix()}",
            f"- sequence/domain queue: {args.sequence_queue_output.as_posix()}",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    meta = {
        "uniprot_rest_endpoint": "https://rest.uniprot.org/uniprotkb/search",
        "fields": UNIPROT_FIELDS,
        "priority": args.priority,
    }
    (args.cache_dir / "run_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
