"""Fetch NCBI linked protein sequences for priority 1 confirmation candidates."""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from xml.etree import ElementTree

import pandas as pd


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def clean_id(value: object) -> str:
    text = str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]
    return text


def cache_path(cache_dir: pathlib.Path, prefix: str, key: str, suffix: str) -> pathlib.Path:
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", str(key))
    return cache_dir / f"{prefix}_{safe}.{suffix}"


def urlopen_text(url: str, timeout: int = 60) -> str:
    last_error: Exception | None = None
    for attempt in range(6):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code == 429:
                time.sleep(10 * (attempt + 1))
                continue
            raise
        except urllib.error.URLError as exc:
            last_error = exc
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"Failed URL after retries: {url}; last_error={last_error}")


def looks_like_ncbi_error(text: str) -> bool:
    return '"ERROR"' in text or "NCBI C++ Exception" in text or "Read failed" in text


def elink_gene_to_protein(gene_id: str, cache_dir: pathlib.Path, sleep: float) -> list[str]:
    path = cache_path(cache_dir, "elink_gene_protein", gene_id, "json")
    if path.exists() and not looks_like_ncbi_error(path.read_text(encoding="utf-8")):
        text = path.read_text(encoding="utf-8")
    else:
        url = f"{EUTILS}/elink.fcgi?" + urllib.parse.urlencode(
            {"dbfrom": "gene", "db": "protein", "id": gene_id, "retmode": "json"}
        )
        last_text = ""
        for attempt in range(4):
            text = urlopen_text(url)
            last_text = text
            if not looks_like_ncbi_error(text):
                break
            time.sleep(2 * (attempt + 1))
        text = last_text
        path.write_text(text, encoding="utf-8")
        if sleep:
            time.sleep(sleep)

    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        return []
    links: list[str] = []
    for linkset in data.get("linksets", []):
        for linkdb in linkset.get("linksetdbs", []):
            if linkdb.get("linkname") == "gene_protein_refseq":
                links.extend(linkdb.get("links", []))
    if not links:
        for linkset in data.get("linksets", []):
            for linkdb in linkset.get("linksetdbs", []):
                if linkdb.get("linkname") == "gene_protein":
                    links.extend(linkdb.get("links", []))
    return list(dict.fromkeys(map(str, links)))


def efetch_protein_fasta(protein_ids: list[str], cache_dir: pathlib.Path, sleep: float) -> str:
    if not protein_ids:
        return ""
    key = "_".join(protein_ids[:20])
    path = cache_path(cache_dir, "efetch_protein_fasta", key, "faa")
    if path.exists():
        return path.read_text(encoding="utf-8")
    url = f"{EUTILS}/efetch.fcgi?" + urllib.parse.urlencode(
        {
            "db": "protein",
            "id": ",".join(protein_ids),
            "rettype": "fasta",
            "retmode": "text",
        }
    )
    text = urlopen_text(url)
    path.write_text(text, encoding="utf-8")
    if sleep:
        time.sleep(sleep)
    return text


def efetch_protein_docsum(protein_ids: list[str], cache_dir: pathlib.Path, sleep: float) -> list[dict[str, str]]:
    if not protein_ids:
        return []
    key = "_".join(protein_ids[:20])
    path = cache_path(cache_dir, "efetch_protein_docsum", key, "xml")
    if path.exists():
        text = path.read_text(encoding="utf-8")
    else:
        url = f"{EUTILS}/efetch.fcgi?" + urllib.parse.urlencode(
            {
                "db": "protein",
                "id": ",".join(protein_ids),
                "rettype": "docsum",
                "retmode": "xml",
            }
        )
        text = urlopen_text(url)
        path.write_text(text, encoding="utf-8")
        if sleep:
            time.sleep(sleep)
    root = ElementTree.fromstring(text)
    rows = []
    for docsum in root.findall(".//DocSum"):
        record: dict[str, str] = {"protein_uid": docsum.findtext("Id", default="")}
        for item in docsum.findall("Item"):
            name = item.attrib.get("Name", "")
            if name:
                record[name] = item.text or ""
        rows.append(record)
    return rows


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


def fetch_for_gene(gene_id: str, cache_dir: pathlib.Path, sleep: float, max_proteins: int) -> tuple[list[dict[str, object]], str]:
    protein_ids = elink_gene_to_protein(gene_id, cache_dir, sleep)
    selected_ids = protein_ids[:max_proteins]
    fasta = efetch_protein_fasta(selected_ids, cache_dir, sleep)
    docsums = efetch_protein_docsum(selected_ids, cache_dir, sleep)
    docsum_by_uid = {row.get("protein_uid", ""): row for row in docsums}
    records = parse_fasta(fasta)
    for protein_id, record in zip(selected_ids, records):
        record["protein_uid"] = protein_id
        record["protein_accession"] = accession_from_header(str(record["fasta_header"]))
        record.update({f"docsum_{k}": v for k, v in docsum_by_uid.get(protein_id, {}).items()})
    return records, fasta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=pathlib.Path, required=True)
    parser.add_argument("--metadata-output", type=pathlib.Path, required=True)
    parser.add_argument("--candidate-fasta-output", type=pathlib.Path, required=True)
    parser.add_argument("--human-fasta-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--cache-dir", type=pathlib.Path, required=True)
    parser.add_argument("--max-proteins-per-gene", type=int, default=5)
    parser.add_argument("--sleep", type=float, default=0.1)
    args = parser.parse_args()

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    batch = pd.read_csv(args.batch, sep="\t", dtype=str)
    rows = []
    candidate_fasta_parts = []
    human_fasta_parts = []

    for _, row in batch.iterrows():
        gene_id = clean_id(row["ncbi_gene_id"])
        records, fasta = fetch_for_gene(gene_id, args.cache_dir, args.sleep, args.max_proteins_per_gene)
        for rank, record in enumerate(records, start=1):
            out = row.to_dict()
            out.update(record)
            out["protein_rank_for_gene_species"] = rank
            out["protein_fetch_status"] = "protein_sequence_found"
            rows.append(out)
        if records:
            tagged = []
            for rank, record in enumerate(records, start=1):
                header = (
                    f">{row['human_gene_symbol']}|{row['scientific_name'].replace(' ', '_')}|"
                    f"gene:{gene_id}|rank:{rank}|{record['fasta_header']}"
                )
                tagged.append(header + "\n" + str(record["sequence"]))
            candidate_fasta_parts.append("\n".join(tagged))
        else:
            out = row.to_dict()
            out.update(
                {
                    "protein_uid": "",
                    "protein_accession": "",
                    "fasta_header": "",
                    "sequence": "",
                    "protein_length": 0,
                    "protein_rank_for_gene_species": "",
                    "protein_fetch_status": "no_linked_protein_sequence",
                }
            )
            rows.append(out)

    human_ref = batch[["human_gene_symbol", "entrezgene"]].drop_duplicates()
    for _, row in human_ref.iterrows():
        gene_id = clean_id(row["entrezgene"])
        records, _ = fetch_for_gene(gene_id, args.cache_dir, args.sleep, args.max_proteins_per_gene)
        for rank, record in enumerate(records, start=1):
            header = f">{row['human_gene_symbol']}|Homo_sapiens|gene:{gene_id}|rank:{rank}|{record['fasta_header']}"
            human_fasta_parts.append(header + "\n" + str(record["sequence"]))

    metadata = pd.DataFrame(rows)
    metadata = metadata.sort_values(
        ["human_gene_symbol", "scientific_name", "protein_rank_for_gene_species"]
    )
    args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(args.metadata_output, sep="\t", index=False)
    args.candidate_fasta_output.parent.mkdir(parents=True, exist_ok=True)
    args.candidate_fasta_output.write_text("\n".join(candidate_fasta_parts) + "\n", encoding="utf-8")
    args.human_fasta_output.parent.mkdir(parents=True, exist_ok=True)
    args.human_fasta_output.write_text("\n".join(human_fasta_parts) + "\n", encoding="utf-8")

    status_counts = metadata["protein_fetch_status"].value_counts().sort_index()
    lines = [
        "# Phase 2 Priority 1 Protein Sequence Fetch Report",
        "",
        "## Summary",
        "",
        f"Input gene-species rows: {len(batch)}",
        f"Output protein rows: {len(metadata)}",
        f"Genes with candidate protein rows: {metadata.loc[metadata['protein_fetch_status'] == 'protein_sequence_found', 'human_gene_symbol'].nunique()}",
        f"Species with candidate protein rows: {metadata.loc[metadata['protein_fetch_status'] == 'protein_sequence_found', 'scientific_name'].nunique()}",
        "",
        "## Fetch Status",
        "",
    ]
    for status, count in status_counts.items():
        lines.append(f"- {status}: {count}")
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- metadata: {args.metadata_output.as_posix()}",
            f"- candidate FASTA: {args.candidate_fasta_output.as_posix()}",
            f"- human reference FASTA: {args.human_fasta_output.as_posix()}",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.metadata_output}")


if __name__ == "__main__":
    main()
