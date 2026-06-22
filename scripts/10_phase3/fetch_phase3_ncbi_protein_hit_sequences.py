"""Fetch protein FASTA sequences for Phase 3 NCBI Protein rescue hits."""

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


def clean_id(value: object) -> str:
    text = str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]
    return text


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
            if exc.code == 429:
                time.sleep(8 * (attempt + 1))
                continue
            raise
        except urllib.error.URLError as exc:
            last_error = exc
            time.sleep(4 * (attempt + 1))
    raise RuntimeError(f"Failed URL after retries: {url}; last_error={last_error}")


def cached_efetch(ids: list[str], cache_dir: pathlib.Path, prefix: str, rettype: str, sleep: float) -> str:
    ids = [clean_id(x) for x in ids if clean_id(x)]
    if not ids:
        return ""
    path = cache_dir / f"{prefix}_{safe_key('_'.join(ids[:20]))}_{rettype}.txt"
    if path.exists():
        return path.read_text(encoding="utf-8")
    url = f"{EUTILS}/efetch.fcgi?" + urllib.parse.urlencode(
        {
            "db": "protein",
            "id": ",".join(ids),
            "rettype": rettype,
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
    gene_id = clean_id(gene_id)
    if not gene_id:
        return []
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
    preferred = ("gene_protein_refseq", "gene_protein")
    for linkname in preferred:
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


def split_accessions(value: object, max_accessions: int) -> list[str]:
    accessions = [x.strip() for x in str(value).split(";") if x.strip()]
    return accessions[:max_accessions]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-review", type=pathlib.Path, required=True)
    parser.add_argument("--eligibility", type=pathlib.Path, required=True)
    parser.add_argument("--metadata-output", type=pathlib.Path, required=True)
    parser.add_argument("--candidate-fasta-output", type=pathlib.Path, required=True)
    parser.add_argument("--human-fasta-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--cache-dir", type=pathlib.Path, required=True)
    parser.add_argument("--max-candidate-accessions", type=int, default=4)
    parser.add_argument("--max-human-proteins", type=int, default=5)
    parser.add_argument("--sleep", type=float, default=0.1)
    args = parser.parse_args()

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    review = pd.read_csv(args.protein_review, sep="\t", dtype=str).fillna("")
    eligibility = pd.read_csv(args.eligibility, sep="\t", dtype=str).fillna("")
    genes = review.loc[
        review["ncbi_protein_review_status"].isin(["protein_symbol_like_hit", "protein_broad_hit"]),
        ["human_gene_symbol", "scientific_name", "ncbi_protein_accessions_top", "ncbi_protein_review_status"],
    ].copy()

    metadata_rows = []
    candidate_parts = []
    for _, row in genes.iterrows():
        accessions = split_accessions(row["ncbi_protein_accessions_top"], args.max_candidate_accessions)
        fasta = cached_efetch(accessions, args.cache_dir, "candidate", "fasta", args.sleep)
        records = parse_fasta(fasta)
        if not records:
            metadata_rows.append(
                {
                    **row.to_dict(),
                    "protein_accession": "",
                    "fasta_header": "",
                    "sequence": "",
                    "protein_length": 0,
                    "protein_fetch_status": "no_sequence_fetched",
                }
            )
            continue
        for rank, record in enumerate(records, start=1):
            header = (
                f">{row['human_gene_symbol']}|{row['scientific_name'].replace(' ', '_')}|"
                f"protein_accession:{accession_from_header(str(record['fasta_header']))}|rank:{rank}|{record['fasta_header']}"
            )
            candidate_parts.append(header + "\n" + str(record["sequence"]))
            metadata_rows.append(
                {
                    **row.to_dict(),
                    "protein_accession": accession_from_header(str(record["fasta_header"])),
                    "fasta_header": record["fasta_header"],
                    "sequence": record["sequence"],
                    "protein_length": record["protein_length"],
                    "protein_fetch_status": "protein_sequence_found",
                }
            )

    human_parts = []
    human_ref = eligibility[
        eligibility["human_gene_symbol"].isin(sorted(set(genes["human_gene_symbol"])))
    ][["human_gene_symbol", "entrezgene"]].drop_duplicates()
    for _, row in human_ref.iterrows():
        protein_ids = cached_elink_gene_to_protein(row["entrezgene"], args.cache_dir, args.sleep)[: args.max_human_proteins]
        fasta = cached_efetch(protein_ids, args.cache_dir, "human_ref", "fasta", args.sleep)
        for rank, record in enumerate(parse_fasta(fasta), start=1):
            header = (
                f">{row['human_gene_symbol']}|Homo_sapiens|gene:{clean_id(row['entrezgene'])}|"
                f"rank:{rank}|{record['fasta_header']}"
            )
            human_parts.append(header + "\n" + str(record["sequence"]))

    metadata = pd.DataFrame(metadata_rows)
    args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
    args.candidate_fasta_output.parent.mkdir(parents=True, exist_ok=True)
    args.human_fasta_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(args.metadata_output, sep="\t", index=False)
    args.candidate_fasta_output.write_text("\n".join(candidate_parts) + "\n", encoding="utf-8")
    args.human_fasta_output.write_text("\n".join(human_parts) + "\n", encoding="utf-8")

    status_counts = metadata["protein_fetch_status"].value_counts().sort_index() if not metadata.empty else {}
    lines = [
        "# Phase 3 NCBI Protein Hit Sequence Fetch Report",
        "",
        f"Protein-review hit rows: {len(genes)}",
        f"Candidate protein sequence rows: {len(metadata)}",
        f"Human reference FASTA entries: {len(human_parts)}",
        "",
        "## Fetch Status",
        "",
    ]
    for status, count in getattr(status_counts, "items", lambda: [])():
        lines.append(f"- {status}: {count}")
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- metadata: `{args.metadata_output.as_posix()}`",
            f"- candidate FASTA: `{args.candidate_fasta_output.as_posix()}`",
            f"- human FASTA: `{args.human_fasta_output.as_posix()}`",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.metadata_output}")


if __name__ == "__main__":
    main()
