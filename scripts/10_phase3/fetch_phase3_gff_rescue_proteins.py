"""Fetch protein sequences for Phase 3 assembly/GFF rescue hits."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import pathlib
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict

import pandas as pd


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def safe_key(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(value))


def parse_attrs(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in text.split(";"):
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        attrs[key] = urllib.parse.unquote(value)
    return attrs


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def urlopen_text(url: str, timeout: int = 60) -> str:
    last_error: Exception | None = None
    for attempt in range(6):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code in {429, 500, 502, 503, 504}:
                time.sleep(5 * (attempt + 1))
                continue
            raise
        except urllib.error.URLError as exc:
            last_error = exc
            time.sleep(4 * (attempt + 1))
    raise RuntimeError(f"Failed URL after retries: {url}; last_error={last_error}")


def cached_efetch_proteins(accessions: list[str], cache_dir: pathlib.Path, prefix: str, sleep: float) -> str:
    accessions = [acc.strip() for acc in accessions if str(acc).strip()]
    if not accessions:
        return ""
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{prefix}_{safe_key('_'.join(accessions[:25]))}_fasta.txt"
    if path.exists():
        return path.read_text(encoding="utf-8")
    url = f"{EUTILS}/efetch.fcgi?" + urllib.parse.urlencode(
        {
            "db": "protein",
            "id": ",".join(accessions),
            "rettype": "fasta",
            "retmode": "text",
            "tool": "bird_lifespan_phase3",
        }
    )
    text = urlopen_text(url)
    path.write_text(text, encoding="utf-8")
    if sleep:
        time.sleep(sleep)
    return text


def cached_elink_gene_to_protein(gene_id: str, cache_dir: pathlib.Path, sleep: float) -> list[str]:
    gene_id = str(gene_id).strip()
    if not gene_id:
        return []
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"elink_gene_protein_{safe_key(gene_id)}.json"
    if path.exists():
        text = path.read_text(encoding="utf-8")
    else:
        url = f"{EUTILS}/elink.fcgi?" + urllib.parse.urlencode(
            {
                "dbfrom": "gene",
                "db": "protein",
                "id": gene_id,
                "retmode": "json",
                "tool": "bird_lifespan_phase3",
            }
        )
        text = urlopen_text(url)
        path.write_text(text, encoding="utf-8")
        if sleep:
            time.sleep(sleep)
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        return []
    links: list[str] = []
    for linkname in ("gene_protein_refseq", "gene_protein"):
        for linkset in data.get("linksets", []):
            for linkdb in linkset.get("linksetdbs", []):
                if linkdb.get("linkname") == linkname:
                    links.extend(map(str, linkdb.get("links", [])))
        if links:
            break
    return list(dict.fromkeys(links))


def parse_fasta(text: str) -> list[dict[str, object]]:
    records = []
    header = None
    seq_parts: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                seq = "".join(seq_parts)
                records.append({"fasta_header": header, "sequence": seq, "protein_length": len(seq)})
            header = line[1:].strip()
            seq_parts = []
        elif line.strip():
            seq_parts.append(line.strip())
    if header is not None:
        seq = "".join(seq_parts)
        records.append({"fasta_header": header, "sequence": seq, "protein_length": len(seq)})
    return records


def accession_from_header(header: str) -> str:
    return str(header).split()[0] if header else ""


def collect_gene_protein_ids(gff_path: pathlib.Path, gene_ids: set[str]) -> dict[str, list[dict[str, str]]]:
    mrna_to_gene: dict[str, str] = {}
    mrna_attrs: dict[str, dict[str, str]] = {}
    proteins: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    with gzip.open(gff_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            ftype = parts[2]
            attrs = parse_attrs(parts[8])
            if ftype in {"mRNA", "transcript"} and attrs.get("Parent") in gene_ids:
                mrna_id = attrs.get("ID", "")
                if mrna_id:
                    mrna_to_gene[mrna_id] = attrs["Parent"]
                    mrna_attrs[mrna_id] = attrs
            elif ftype == "CDS" and attrs.get("Parent") in mrna_to_gene:
                gene_id = mrna_to_gene[attrs["Parent"]]
                protein_id = attrs.get("protein_id") or attrs.get("Name", "")
                if not protein_id and attrs.get("Dbxref", "").startswith("NCBI_GP:"):
                    protein_id = attrs.get("Dbxref", "").split("NCBI_GP:", 1)[1].split(",", 1)[0]
                if not protein_id:
                    continue
                proteins[gene_id][protein_id] = {
                    "gff_gene_id": gene_id,
                    "gff_mrna_id": attrs.get("Parent", ""),
                    "gff_cds_id": attrs.get("ID", ""),
                    "protein_id": protein_id,
                    "cds_name": attrs.get("Name", ""),
                    "cds_dbxref": attrs.get("Dbxref", ""),
                    "cds_product": attrs.get("product", ""),
                    "mrna_product": mrna_attrs.get(attrs.get("Parent", ""), {}).get("product", ""),
                    "partial": attrs.get("partial", mrna_attrs.get(attrs.get("Parent", ""), {}).get("partial", "")),
                }
    return {gene_id: list(records.values()) for gene_id, records in proteins.items()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff-rescue-hits", type=pathlib.Path, required=True)
    parser.add_argument("--eligibility", type=pathlib.Path, required=True)
    parser.add_argument("--metadata-output", type=pathlib.Path, required=True)
    parser.add_argument("--candidate-fasta-output", type=pathlib.Path, required=True)
    parser.add_argument("--human-fasta-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--cache-dir", type=pathlib.Path, required=True)
    parser.add_argument("--max-human-proteins", type=int, default=5)
    parser.add_argument("--sleep", type=float, default=0.1)
    args = parser.parse_args()

    hits = pd.read_csv(args.gff_rescue_hits, sep="\t", dtype=str).fillna("")
    eligible_hits = hits[
        (hits["gff_can_count_as_annotation_rescue"] == "True")
        & (hits["gff_rescue_call"].isin(["gff_symbol_exact", "gff_symbol_like", "gff_alias_hit", "gff_product_hit"]))
        & (hits["matched_id"] != "")
    ].copy()

    grouped = eligible_hits.groupby("gff_local_path")
    extracted: dict[tuple[str, str], list[dict[str, str]]] = {}
    for gff_path_text, group in grouped:
        gff_path = pathlib.Path(gff_path_text)
        gene_ids = set(group["matched_id"])
        gene_proteins = collect_gene_protein_ids(gff_path, gene_ids) if gff_path.exists() else {}
        for _, row in group.iterrows():
            extracted[(row["scientific_name"], row["human_gene_symbol"])] = gene_proteins.get(row["matched_id"], [])

    candidate_parts: list[str] = []
    metadata_rows: list[dict[str, object]] = []
    for _, row in eligible_hits.iterrows():
        records = extracted.get((row["scientific_name"], row["human_gene_symbol"]), [])
        if not records:
            metadata_rows.append(
                {
                    **row.to_dict(),
                    "protein_id": "",
                    "gff_mrna_id": "",
                    "gff_cds_id": "",
                    "cds_product": "",
                    "mrna_product": "",
                    "partial": "",
                    "sequence": "",
                    "protein_length": 0,
                    "sequence_fetch_status": "no_protein_id_in_gff_cds",
                    "fasta_header": "",
                    "fetched_accession": "",
                }
            )
            continue
        accessions = [record["protein_id"] for record in records]
        fasta = cached_efetch_proteins(accessions, args.cache_dir, f"gff_candidate_{row['scientific_name']}_{row['human_gene_symbol']}", args.sleep)
        fetched = {accession_from_header(str(record["fasta_header"])): record for record in parse_fasta(fasta)}
        for rank, record in enumerate(records, start=1):
            protein_id = record["protein_id"]
            fasta_record = fetched.get(protein_id)
            if fasta_record is None:
                metadata_rows.append(
                    {
                        **row.to_dict(),
                        **record,
                        "sequence": "",
                        "protein_length": 0,
                        "sequence_fetch_status": "protein_id_efetch_failed",
                        "fasta_header": "",
                        "fetched_accession": "",
                    }
                )
                continue
            safe_species = row["scientific_name"].replace(" ", "_")
            header = (
                f">{row['human_gene_symbol']}|{safe_species}|assembly_gff_protein:{protein_id}|"
                f"rank:{rank}|{protein_id} {fasta_record['fasta_header']}"
            )
            candidate_parts.append(header + "\n" + str(fasta_record["sequence"]))
            metadata_rows.append(
                {
                    **row.to_dict(),
                    **record,
                    "sequence": fasta_record["sequence"],
                    "protein_length": fasta_record["protein_length"],
                    "sequence_fetch_status": "protein_sequence_found",
                    "fasta_header": fasta_record["fasta_header"],
                    "fetched_accession": protein_id,
                }
            )

    metadata = pd.DataFrame(metadata_rows)
    eligibility = pd.read_csv(args.eligibility, sep="\t", dtype=str).fillna("")
    genes = sorted(set(metadata.loc[metadata["sequence_fetch_status"] == "protein_sequence_found", "human_gene_symbol"])) if not metadata.empty else []
    human_parts: list[str] = []
    human_ref = eligibility[eligibility["human_gene_symbol"].isin(genes)][["human_gene_symbol", "entrezgene"]].drop_duplicates()
    for _, row in human_ref.iterrows():
        protein_ids = cached_elink_gene_to_protein(row["entrezgene"], args.cache_dir, args.sleep)[: args.max_human_proteins]
        fasta = cached_efetch_proteins(protein_ids, args.cache_dir, f"human_ref_{row['human_gene_symbol']}", args.sleep)
        for rank, record in enumerate(parse_fasta(fasta), start=1):
            header = (
                f">{row['human_gene_symbol']}|Homo_sapiens|gene:{row['entrezgene']}|"
                f"rank:{rank}|{accession_from_header(str(record['fasta_header']))} {record['fasta_header']}"
            )
            human_parts.append(header + "\n" + str(record["sequence"]))

    args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
    args.candidate_fasta_output.parent.mkdir(parents=True, exist_ok=True)
    args.human_fasta_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(args.metadata_output, sep="\t", index=False)
    args.candidate_fasta_output.write_text("\n".join(candidate_parts) + ("\n" if candidate_parts else ""), encoding="utf-8")
    args.human_fasta_output.write_text("\n".join(human_parts) + ("\n" if human_parts else ""), encoding="utf-8")

    status_counts = metadata["sequence_fetch_status"].value_counts().sort_index() if not metadata.empty else {}
    lines = [
        "# Phase 3 GFF Rescue Protein Fetch Report",
        "",
        f"GFF annotation rescue rows eligible for sequence fetch: {len(eligible_hits)}",
        f"Metadata rows: {len(metadata)}",
        f"Candidate FASTA records: {len(candidate_parts)}",
        f"Human reference FASTA records: {len(human_parts)}",
        "",
        "## Fetch Status",
    ]
    for status, count in getattr(status_counts, "items", lambda: [])():
        lines.append(f"- {status}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "Candidate proteins fetched from GFF-linked CDS `protein_id` values are ready for reciprocal DIAMOND/BLASTP validation. They should not be counted as strict sequence support until reciprocal same-gene validation passes.",
            "",
            "## Outputs",
            f"- metadata: `{args.metadata_output}`",
            f"- candidate FASTA: `{args.candidate_fasta_output}`",
            f"- human FASTA: `{args.human_fasta_output}`",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.metadata_output}")


if __name__ == "__main__":
    main()
