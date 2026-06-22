"""UniProtKB cross-check for Phase 3 priority-1 rescue rows."""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd


UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

PRODUCT_HINTS = {
    "DNMT1": ["dnmt1", "dna (cytosine-5)-methyltransferase", "methyltransferase 1"],
    "DNMT3A": ["dnmt3a", "methyltransferase 3a", "methyltransferase 3 alpha"],
    "DNMT3B": ["dnmt3b", "methyltransferase 3b", "methyltransferase 3 beta"],
    "HELLS": ["hells", "lymphoid-specific helicase"],
    "MBD2": ["mbd2", "methyl-cpg-binding domain protein 2"],
    "MBD3": ["mbd3", "methyl-cpg-binding domain protein 3"],
    "MORC3": ["morc3", "morc family cw-type zinc finger 3"],
    "SAMHD1": ["samhd1", "sam and hd domain-containing protein 1", "deoxynucleoside triphosphate triphosphohydrolase"],
    "SETDB2": ["setdb2", "histone-lysine n-methyltransferase setdb2"],
    "UHRF1": ["uhrf1", "e3 ubiquitin-protein ligase uhrf1"],
}

FAMILY_CONFLICTS = {
    "DNMT1": ["dnmt3a", "dnmt3b"],
    "DNMT3A": ["dnmt1", "dnmt3b"],
    "DNMT3B": ["dnmt1", "dnmt3a"],
    "MBD2": ["mbd3", "lmbd2"],
    "MBD3": ["mbd2", "lmbd3"],
    "SETDB2": ["setdb1"],
}


def safe_key(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(value))


def urlopen_text(url: str, timeout: int = 45) -> str:
    last_error: Exception | None = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code in {429, 500, 502, 503, 504}:
                time.sleep(4 * (attempt + 1))
                continue
            raise
        except urllib.error.URLError as exc:
            last_error = exc
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"Failed URL after retries: {url}; last_error={last_error}")


def query_uniprot(taxid: str, symbol: str, cache_dir: pathlib.Path, sleep: float, size: int = 5) -> str:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"uniprot_{safe_key(taxid)}_{safe_key(symbol)}.tsv"
    if path.exists():
        return path.read_text(encoding="utf-8")
    query = f"(organism_id:{taxid}) AND (gene:{symbol} OR gene_exact:{symbol} OR protein_name:{symbol})"
    params = {
        "query": query,
        "fields": "accession,id,gene_names,protein_name,organism_id,organism_name,xref_refseq",
        "format": "tsv",
        "size": str(size),
    }
    url = UNIPROT_SEARCH + "?" + urllib.parse.urlencode(params)
    text = urlopen_text(url)
    path.write_text(text, encoding="utf-8")
    if sleep:
        time.sleep(sleep)
    return text


def parse_tsv(text: str) -> list[dict[str, str]]:
    lines = [line for line in text.splitlines() if line.strip() and not line.startswith("% ")]
    if len(lines) < 2:
        return []
    reader = csv.DictReader(lines, delimiter="\t")
    return [dict(row) for row in reader]


def norm(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(text).lower())


def token_match(symbol: str, text: str) -> bool:
    return re.search(rf"(?<![A-Za-z0-9]){re.escape(symbol)}(?![A-Za-z0-9])", text, re.IGNORECASE) is not None


def classify(symbol: str, hits: list[dict[str, str]], local_protein_id: str) -> tuple[str, str, dict[str, str]]:
    if not hits:
        return "no_uniprot_hit", "No UniProtKB hit returned for taxid+gene query.", {}
    symbol_norm = norm(symbol)
    best = hits[0]
    best_text = " ".join(str(v) for v in best.values())
    lower = best_text.lower()
    conflicts = [item for item in FAMILY_CONFLICTS.get(symbol, []) if item in norm(best_text)]
    gene_names = best.get("Gene Names", "")
    protein_name = best.get("Protein names", "")
    refseq = best.get("RefSeq", "")
    protein_lower = protein_name.lower()
    if symbol.lower() in protein_lower and re.search(rf"{re.escape(symbol.lower())}[-\s]?interacting", protein_lower):
        return "uniprot_interactor_not_target", "UniProt top protein is an interacting protein mentioning the requested symbol, not the target gene itself.", best
    if local_protein_id and local_protein_id in refseq:
        return "uniprot_refseq_accession_match", "Local GFF protein_id appears in UniProt RefSeq/xref field.", best
    if any(norm(token) == symbol_norm for token in re.split(r"[\s;,]+", gene_names) if token):
        return "uniprot_gene_exact_support", "UniProt gene_names contain exact requested symbol.", best
    if token_match(symbol, gene_names):
        return "uniprot_symbol_like_support", "UniProt gene/protein fields contain requested symbol-like token.", best
    if token_match(symbol, protein_name) and "interacting" not in protein_lower:
        return "uniprot_symbol_like_support", "UniProt protein field contains requested symbol-like token.", best
    for hint in PRODUCT_HINTS.get(symbol, []):
        if hint in lower:
            if conflicts:
                return "uniprot_family_conflict", f"Product hint matched but family conflict detected: {','.join(conflicts)}.", best
            return "uniprot_product_support", f"UniProt protein name matches product hint: {hint}.", best
    if conflicts:
        return "uniprot_family_conflict", f"UniProt top hit contains related-family signal: {','.join(conflicts)}.", best
    return "uniprot_broad_hit_manual_review", "UniProt returned hits but no exact symbol/product rule matched.", best


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--cache-dir", type=pathlib.Path, required=True)
    parser.add_argument("--sleep", type=float, default=0.05)
    args = parser.parse_args()

    audit = pd.read_csv(args.audit, sep="\t", dtype=str).fillna("")
    rows = []
    for _, row in audit.iterrows():
        text = query_uniprot(row["species_taxid"], row["human_gene_symbol"], args.cache_dir, args.sleep)
        hits = parse_tsv(text)
        call, reason, best = classify(row["human_gene_symbol"], hits, row.get("protein_id", ""))
        rows.append(
            {
                "scientific_name": row["scientific_name"],
                "human_gene_symbol": row["human_gene_symbol"],
                "species_taxid": row["species_taxid"],
                "next_rescue_route": row["next_rescue_route"],
                "phase3_gff_sequence_decision": row["phase3_gff_sequence_decision"],
                "can_count_as_strict_sequence_after_gff_sequence": row["can_count_as_strict_sequence_after_gff_sequence"],
                "local_gff_protein_id": row.get("protein_id", ""),
                "uniprot_hit_count": len(hits),
                "uniprot_crosscheck_call": call,
                "uniprot_crosscheck_reason": reason,
                "uniprot_accession_top": best.get("Entry", ""),
                "uniprot_entry_name_top": best.get("Entry Name", ""),
                "uniprot_gene_names_top": best.get("Gene Names", ""),
                "uniprot_protein_name_top": best.get("Protein names", ""),
                "uniprot_organism_id_top": best.get("Organism (ID)", ""),
                "uniprot_organism_top": best.get("Organism", ""),
                "uniprot_refseq_top": best.get("RefSeq", ""),
            }
        )
    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    summary = (
        out.groupby(["next_rescue_route", "uniprot_crosscheck_call"], as_index=False)
        .agg(rows=("human_gene_symbol", "count"), species=("scientific_name", "nunique"), genes=("human_gene_symbol", "nunique"))
        .sort_values(["next_rescue_route", "uniprot_crosscheck_call"])
    )
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary_output, sep="\t", index=False)

    call_counts = out["uniprot_crosscheck_call"].value_counts().sort_index()
    strict = out[out["next_rescue_route"] == "resolved_strict_sequence"]
    strict_support = strict[strict["uniprot_crosscheck_call"].isin(["uniprot_refseq_accession_match", "uniprot_gene_exact_support", "uniprot_symbol_like_support", "uniprot_product_support"])]
    external_gap = out[out["next_rescue_route"] == "external_db_or_assembly_reannotation"]
    external_support = external_gap[external_gap["uniprot_crosscheck_call"].isin(["uniprot_gene_exact_support", "uniprot_symbol_like_support", "uniprot_product_support"])]
    lines = [
        "# Phase 3 Priority-1 UniProt Crosscheck Report",
        "",
        f"Rows queried: {len(out)}",
        f"Strict sequence rows with UniProt support: {len(strict_support)} / {len(strict)}",
        f"Local-GFF-missing rows with UniProt support: {len(external_support)} / {len(external_gap)}",
        "",
        "## UniProt Calls",
    ]
    for call, count in call_counts.items():
        lines.append(f"- {call}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "UniProt support is treated as independent annotation/cross-database evidence, not as a replacement for reciprocal sequence validation. Rows with UniProt support but no local GFF hit should be queued for updated assembly annotation or sequence-level confirmation before strict scoring.",
            "",
            "## Outputs",
            f"- row table: `{args.output}`",
            f"- summary: `{args.summary_output}`",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
