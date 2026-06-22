"""Fetch reference bird protein sequences for seed maintenance genes."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import time
import urllib.parse
import urllib.request


DEFAULT_MATRIX = pathlib.Path("data/processed/ortholog_matrix_primary_combined_candidates.tsv")
DEFAULT_SPECIES = "Taeniopygia guttata"
DEFAULT_OUTPUT = pathlib.Path("data/interim/reference_bird_proteins_taeniopygia.tsv")
DEFAULT_FASTA = pathlib.Path("data/interim/reference_bird_proteins_taeniopygia.faa")
DEFAULT_REPORT = pathlib.Path("results/reports/reference_bird_proteins_taeniopygia_report.md")

ELINK_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    fields = list(rows[0].keys()) if rows else [
        "human_gene_symbol",
        "reference_gene_id",
        "reference_protein_uid",
        "reference_protein_accession",
        "sequence_length",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fetch_json(url: str, params: dict[str, str]) -> dict:
    full_url = f"{url}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(full_url, timeout=90) as response:
        text = response.read().decode("utf-8", errors="replace")
        return json.JSONDecoder(strict=False).decode(text)


def fetch_text(url: str, params: dict[str, str]) -> str:
    full_url = f"{url}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(full_url, timeout=90) as response:
        return response.read().decode("utf-8")


def parse_fasta(text: str) -> list[dict[str, str]]:
    records = []
    header = ""
    seq_parts: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header:
                records.append({"header": header, "sequence": "".join(seq_parts)})
            header = line[1:].strip()
            seq_parts = []
        elif line.strip():
            seq_parts.append(line.strip())
    if header:
        records.append({"header": header, "sequence": "".join(seq_parts)})
    return records


def accession_from_header(header: str) -> str:
    return header.split()[0]


def refseq_protein_uids(gene_id: str) -> list[str]:
    payload = fetch_json(
        ELINK_URL,
        {"dbfrom": "gene", "db": "protein", "id": gene_id, "retmode": "json", "tool": "bird_lifespan_feasibility"},
    )
    for linkset in payload.get("linksets", []):
        for linksetdb in linkset.get("linksetdbs", []):
            if linksetdb.get("linkname") == "gene_protein_refseq":
                return linksetdb.get("links", [])
    for linkset in payload.get("linksets", []):
        for linksetdb in linkset.get("linksetdbs", []):
            if linksetdb.get("linkname") == "gene_protein":
                return linksetdb.get("links", [])
    return []


def fetch_proteins(uids: list[str]) -> list[dict[str, str]]:
    if not uids:
        return []
    text = fetch_text(
        EFETCH_URL,
        {"db": "protein", "id": ",".join(uids), "rettype": "fasta", "retmode": "text", "tool": "bird_lifespan_feasibility"},
    )
    return parse_fasta(text)


def reference_gene_rows(matrix_rows: list[dict[str, str]], species: str) -> list[dict[str, str]]:
    rows = [
        row
        for row in matrix_rows
        if row["scientific_name"] == species
        and row["combined_candidate_status"] == "ncbi_gene_candidate"
        and row["ortholog_gene_id"]
    ]
    rows.sort(key=lambda row: row["human_gene_symbol"])
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, default=DEFAULT_MATRIX)
    parser.add_argument("--species", default=DEFAULT_SPECIES)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--fasta", type=pathlib.Path, default=DEFAULT_FASTA)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    matrix_rows = read_tsv(args.matrix)
    refs = reference_gene_rows(matrix_rows, args.species)
    output_rows = []
    fasta_lines = []
    for row in refs:
        uids = refseq_protein_uids(row["ortholog_gene_id"])
        time.sleep(0.35)
        proteins = fetch_proteins(uids)
        time.sleep(0.35)
        if not proteins:
            output_rows.append(
                {
                    "human_gene_symbol": row["human_gene_symbol"],
                    "maintenance_module": row["maintenance_module"],
                    "reference_species": args.species,
                    "reference_gene_id": row["ortholog_gene_id"],
                    "reference_protein_uid": "",
                    "reference_protein_accession": "",
                    "sequence_length": "0",
                    "fetch_status": "no_refseq_protein_link",
                    "header": "",
                    "sequence": "",
                }
            )
            continue
        best = max(proteins, key=lambda record: len(record["sequence"]))
        uid = uids[0] if uids else ""
        output_rows.append(
            {
                "human_gene_symbol": row["human_gene_symbol"],
                "maintenance_module": row["maintenance_module"],
                "reference_species": args.species,
                "reference_gene_id": row["ortholog_gene_id"],
                "reference_protein_uid": uid,
                "reference_protein_accession": accession_from_header(best["header"]),
                "sequence_length": str(len(best["sequence"])),
                "fetch_status": "ok",
                "header": best["header"],
                "sequence": best["sequence"],
            }
        )
        fasta_lines.append(f">{row['human_gene_symbol']}|{accession_from_header(best['header'])}|{args.species}\n{best['sequence']}")

    write_tsv(args.output, output_rows)
    args.fasta.parent.mkdir(parents=True, exist_ok=True)
    args.fasta.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    ok = [row for row in output_rows if row["fetch_status"] == "ok"]
    missing = [row for row in output_rows if row["fetch_status"] != "ok"]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        "\n".join(
            [
                "# Reference Bird Proteins Report",
                "",
                f"Reference species: {args.species}",
                f"Candidate genes from matrix: {len(refs)}",
                f"Proteins fetched: {len(ok)}",
                f"Missing protein links: {len(missing)}",
                "",
                "## Interpretation",
                "These reference bird proteins are used only for local protein-similarity rescue of poorly named assembly annotations. They are not a final orthology database.",
            ]
        ),
        encoding="utf-8",
    )
    print(f"Wrote {args.output}, {args.fasta}, and {args.report}")


if __name__ == "__main__":
    main()
