"""Fetch reference proteins for Week 4 transposon reciprocal searches."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import time
import urllib.parse
import urllib.request


DEFAULT_MATRIX = pathlib.Path("data/processed/ortholog_matrix_primary_diamond_validated_candidates.tsv")
DEFAULT_BATCH = pathlib.Path("data/processed/week4_transposon_validation_batch.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/interim/week4_transposon_reference_proteins.tsv")
DEFAULT_FASTA = pathlib.Path("data/interim/week4_transposon_reference_proteins.faa")
DEFAULT_REPORT = pathlib.Path("results/reports/week4_transposon_reference_proteins_report.md")
DEFAULT_REFERENCE_SPECIES = ["Taeniopygia guttata", "Anolis carolinensis"]

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ELINK_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    fields = list(rows[0].keys()) if rows else ["human_gene_symbol", "reference_species", "fetch_status"]
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
        return response.read().decode("utf-8", errors="replace")


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
    token = header.split()[0]
    if "|" in token:
        parts = token.split("|")
        if len(parts) >= 4:
            return parts[3]
    return token


def refseq_protein_uids(gene_id: str) -> list[str]:
    if not gene_id:
        return []
    payload = fetch_json(
        ELINK_URL,
        {"dbfrom": "gene", "db": "protein", "id": gene_id, "retmode": "json", "tool": "bird_lifespan_week4"},
    )
    for preferred in ["gene_protein_refseq", "gene_protein"]:
        for linkset in payload.get("linksets", []):
            for linksetdb in linkset.get("linksetdbs", []):
                if linksetdb.get("linkname") == preferred:
                    return linksetdb.get("links", [])
    return []


def protein_search_uids(gene: str, species: str) -> list[str]:
    terms = [
        f'({gene}[Gene Name]) AND "{species}"[Organism] AND srcdb_refseq[prop]',
        f'({gene}[Gene Name]) AND "{species}"[Organism]',
        f'({gene}[All Fields]) AND "{species}"[Organism] AND srcdb_refseq[prop]',
    ]
    for term in terms:
        payload = fetch_json(
            ESEARCH_URL,
            {
                "db": "protein",
                "term": term,
                "retmode": "json",
                "retmax": "20",
                "tool": "bird_lifespan_week4",
            },
        )
        ids = payload.get("esearchresult", {}).get("idlist", [])
        if ids:
            return ids
    return []


def fetch_proteins(uids: list[str]) -> list[dict[str, str]]:
    if not uids:
        return []
    text = fetch_text(
        EFETCH_URL,
        {
            "db": "protein",
            "id": ",".join(uids[:20]),
            "rettype": "fasta",
            "retmode": "text",
            "tool": "bird_lifespan_week4",
        },
    )
    return parse_fasta(text)


def choose_best(records: list[dict[str, str]], gene: str) -> dict[str, str] | None:
    if not records:
        return None
    gene_upper = gene.upper()
    preferred = [
        rec
        for rec in records
        if gene_upper in rec["header"].upper() and "ISOFORM X" not in rec["header"].upper()
    ]
    pool = preferred if preferred else records
    return max(pool, key=lambda rec: len(rec["sequence"]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, default=DEFAULT_MATRIX)
    parser.add_argument("--batch", type=pathlib.Path, default=DEFAULT_BATCH)
    parser.add_argument("--reference-species", nargs="*", default=DEFAULT_REFERENCE_SPECIES)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--fasta", type=pathlib.Path, default=DEFAULT_FASTA)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    matrix = read_tsv(args.matrix)
    batch = read_tsv(args.batch)
    genes = sorted({row["human_gene_symbol"] for row in batch})
    matrix_lookup = {
        (row["scientific_name"], row["human_gene_symbol"]): row
        for row in matrix
        if row["scientific_name"] in args.reference_species
    }

    output_rows = []
    fasta_lines = []
    for species in args.reference_species:
        for gene in genes:
            row = matrix_lookup.get((species, gene), {})
            gene_id = row.get("ortholog_gene_id", "")
            uids = refseq_protein_uids(gene_id)
            method = "gene_elink"
            time.sleep(0.34)
            if not uids:
                uids = protein_search_uids(gene, species)
                method = "protein_esearch"
                time.sleep(0.34)
            proteins = fetch_proteins(uids)
            time.sleep(0.34)
            best = choose_best(proteins, gene)
            if best is None:
                output_rows.append(
                    {
                        "human_gene_symbol": gene,
                        "reference_species": species,
                        "reference_gene_id": gene_id,
                        "reference_protein_uid": "",
                        "reference_protein_accession": "",
                        "sequence_length": "0",
                        "fetch_method": method,
                        "fetch_status": "missing",
                        "header": "",
                        "sequence": "",
                    }
                )
                continue
            accession = accession_from_header(best["header"])
            output_rows.append(
                {
                    "human_gene_symbol": gene,
                    "reference_species": species,
                    "reference_gene_id": gene_id,
                    "reference_protein_uid": uids[0] if uids else "",
                    "reference_protein_accession": accession,
                    "sequence_length": str(len(best["sequence"])),
                    "fetch_method": method,
                    "fetch_status": "ok",
                    "header": best["header"],
                    "sequence": best["sequence"],
                }
            )
            fasta_lines.append(f">{gene}|{accession}|{species}\n{best['sequence']}")

    write_tsv(args.output, output_rows)
    args.fasta.parent.mkdir(parents=True, exist_ok=True)
    args.fasta.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")

    ok = [row for row in output_rows if row["fetch_status"] == "ok"]
    missing = [row for row in output_rows if row["fetch_status"] != "ok"]
    lines = [
        "# Week 4 Transposon Reference Proteins Report",
        "",
        f"Reference species: {', '.join(args.reference_species)}",
        f"Genes requested: {len(genes)}",
        f"Reference protein rows: {len(output_rows)}",
        f"Proteins fetched: {len(ok)}",
        f"Missing: {len(missing)}",
        "",
        "## Missing",
        "",
    ]
    if missing:
        for row in missing:
            lines.append(f"- {row['reference_species']} / {row['human_gene_symbol']}")
    else:
        lines.append("- None")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "These proteins are query references for Week 4 reciprocal sequence validation. They are not final curated orthology assertions.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}, {args.fasta}, and {args.report}")


if __name__ == "__main__":
    main()
