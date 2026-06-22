"""Fetch UniProt candidate sequences and human references for Phase 3 rescue."""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"


def safe_key(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(value))


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
    token = str(header).split()[0]
    if "|" in token:
        parts = token.split("|")
        if len(parts) >= 2:
            return parts[1]
    return token


def cached_uniprot_fasta(accession: str, cache_dir: pathlib.Path, sleep: float) -> str:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"uniprot_{safe_key(accession)}.faa"
    if path.exists():
        return path.read_text(encoding="utf-8")
    text = urlopen_text(UNIPROT_FASTA.format(accession=urllib.parse.quote(accession)))
    path.write_text(text, encoding="utf-8")
    if sleep:
        time.sleep(sleep)
    return text


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=pathlib.Path, required=True)
    parser.add_argument("--eligibility", type=pathlib.Path, required=True)
    parser.add_argument("--metadata-output", type=pathlib.Path, required=True)
    parser.add_argument("--candidate-fasta-output", type=pathlib.Path, required=True)
    parser.add_argument("--human-fasta-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--cache-dir", type=pathlib.Path, required=True)
    parser.add_argument("--max-human-proteins", type=int, default=5)
    parser.add_argument("--sleep", type=float, default=0.1)
    args = parser.parse_args()

    batch = pd.read_csv(args.batch, sep="\t", dtype=str).fillna("")
    metadata_rows: list[dict[str, object]] = []
    candidate_parts: list[str] = []
    for _, row in batch.iterrows():
        accession = row["uniprot_accession_top"]
        try:
            fasta = cached_uniprot_fasta(accession, args.cache_dir, args.sleep)
            records = parse_fasta(fasta)
        except Exception as exc:  # noqa: BLE001 - preserve row-level fetch failures in metadata.
            records = []
            fetch_error = str(exc)
        else:
            fetch_error = ""
        if not records:
            metadata_rows.append(
                {
                    **row.to_dict(),
                    "fetched_accession": accession,
                    "fasta_header": "",
                    "sequence": "",
                    "protein_length": 0,
                    "sequence_fetch_status": "uniprot_fetch_failed",
                    "sequence_fetch_error": fetch_error,
                }
            )
            continue
        record = records[0]
        safe_species = row["scientific_name"].replace(" ", "_")
        header = (
            f">{row['human_gene_symbol']}|{safe_species}|uniprot_accession:{accession}|"
            f"rank:1|{accession} {record['fasta_header']}"
        )
        candidate_parts.append(header + "\n" + str(record["sequence"]))
        metadata_rows.append(
            {
                **row.to_dict(),
                "fetched_accession": accession,
                "fasta_header": record["fasta_header"],
                "sequence": record["sequence"],
                "protein_length": record["protein_length"],
                "sequence_fetch_status": "uniprot_sequence_found",
                "sequence_fetch_error": "",
            }
        )

    metadata = pd.DataFrame(metadata_rows)
    eligibility = pd.read_csv(args.eligibility, sep="\t", dtype=str).fillna("")
    genes = sorted(set(metadata.loc[metadata["sequence_fetch_status"] == "uniprot_sequence_found", "human_gene_symbol"])) if not metadata.empty else []
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
    length_num = pd.to_numeric(metadata.get("protein_length", pd.Series(dtype=str)), errors="coerce")
    lines = [
        "# Phase 3 Full-Length UniProt Rescue Sequence Fetch Report",
        "",
        f"Rows in UniProt rescue batch: {len(batch)}",
        f"Metadata rows: {len(metadata)}",
        f"Candidate FASTA records: {len(candidate_parts)}",
        f"Human reference FASTA records: {len(human_parts)}",
        f"Median fetched candidate length: {length_num.median() if not length_num.empty else ''}",
        "",
        "## Fetch Status",
    ]
    for status, count in getattr(status_counts, "items", lambda: [])():
        lines.append(f"- {status}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "These candidate sequences come from UniProt target-species records and are intended as an external full-length sensitivity rescue. Reciprocal same-gene DIAMOND/BLASTP validation is required before any row is promoted in a sensitivity matrix.",
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
