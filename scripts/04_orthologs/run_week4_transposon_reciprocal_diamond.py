"""Run reciprocal-style DIAMOND validation for Week 4 transposon batch."""

from __future__ import annotations

import argparse
import csv
import gzip
import pathlib
import subprocess
from collections import defaultdict


DEFAULT_BATCH = pathlib.Path("data/processed/week4_transposon_validation_batch.tsv")
DEFAULT_DOWNLOAD_LOG = pathlib.Path("data/interim/week4_transposon_protein_fasta_download_log.tsv")
DEFAULT_REFERENCE_FASTA = pathlib.Path("data/interim/week4_transposon_reference_proteins.faa")
DEFAULT_WORK_DIR = pathlib.Path("data/interim/week4_diamond")
DEFAULT_FORWARD = pathlib.Path("results/tables/week4_transposon_forward_diamond_hits.tsv")
DEFAULT_RECIPROCAL = pathlib.Path("results/tables/week4_transposon_reciprocal_diamond_hits.tsv")
DEFAULT_VALIDATION = pathlib.Path("data/processed/week4_transposon_reciprocal_validation.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week4_transposon_reciprocal_validation_report.md")
DEFAULT_DIAMOND = pathlib.Path("env/tools/diamond-2.2.1/diamond.exe")

OUTFMT = [
    "6",
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qcovhsp",
    "scovhsp",
    "evalue",
    "bitscore",
    "stitle",
]


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def open_maybe_gzip(path: pathlib.Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def parse_fasta(path: pathlib.Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    with open_maybe_gzip(path) as handle:
        header = ""
        seq_parts: list[str] = []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    seqid = header.split()[0]
                    records[seqid] = (header, "".join(seq_parts))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header:
            seqid = header.split()[0]
            records[seqid] = (header, "".join(seq_parts))
    return records


def query_gene(qseqid: str) -> str:
    return qseqid.split("|")[0]


def run_command(command: list[str]) -> None:
    subprocess.run(command, check=True, capture_output=True, text=True)


def make_db(diamond: pathlib.Path, fasta: pathlib.Path, db: pathlib.Path) -> None:
    db.parent.mkdir(parents=True, exist_ok=True)
    dmnd = pathlib.Path(str(db) + ".dmnd")
    if dmnd.exists() and dmnd.stat().st_size > 0:
        return
    run_command([str(diamond), "makedb", "--in", str(fasta), "--db", str(db)])


def blastp(diamond: pathlib.Path, query: pathlib.Path, db: pathlib.Path, output: pathlib.Path, max_target: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            str(diamond),
            "blastp",
            "--query",
            str(query),
            "--db",
            str(db),
            "--out",
            str(output),
            "--outfmt",
            *OUTFMT,
            "--max-target-seqs",
            str(max_target),
            "--evalue",
            "1e-5",
            "--threads",
            "4",
            "--quiet",
        ]
    )


def parse_diamond(path: pathlib.Path) -> list[dict[str, str]]:
    fields = OUTFMT[1:]
    rows = []
    if not path.exists():
        return rows
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for parts in reader:
            if not parts:
                continue
            if len(parts) < len(fields):
                parts = parts + [""] * (len(fields) - len(parts))
            rows.append(dict(zip(fields, parts, strict=False)))
    return rows


def best_by_gene(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    best: dict[str, dict[str, str]] = {}
    for row in rows:
        gene = query_gene(row["qseqid"])
        bitscore = float(row["bitscore"])
        if gene not in best or bitscore > float(best[gene]["bitscore"]):
            best[gene] = row
    return best


def write_candidate_fasta(
    path: pathlib.Path,
    batch_rows: list[dict[str, str]],
    best_forward: dict[tuple[str, str], dict[str, str]],
    fasta_by_accession: dict[str, pathlib.Path],
) -> dict[tuple[str, str, str], str]:
    needed_by_accession: dict[str, set[str]] = defaultdict(set)
    for row in batch_rows:
        key = (row["best_assembly_accession"], row["human_gene_symbol"])
        hit = best_forward.get(key)
        if hit:
            needed_by_accession[row["best_assembly_accession"]].add(hit["sseqid"])

    sequence_headers: dict[tuple[str, str, str], str] = {}
    lines = []
    for accession, seqids in needed_by_accession.items():
        records = parse_fasta(fasta_by_accession[accession])
        for seqid in sorted(seqids):
            if seqid not in records:
                continue
            original_header, sequence = records[seqid]
            for row in batch_rows:
                key = (row["best_assembly_accession"], row["human_gene_symbol"])
                hit = best_forward.get(key)
                if row["best_assembly_accession"] == accession and hit and hit["sseqid"] == seqid:
                    candidate_id = f"{accession}|{row['human_gene_symbol']}|{seqid}"
                    sequence_headers[(accession, row["human_gene_symbol"], seqid)] = original_header
                    lines.append(f">{candidate_id} {original_header}\n{sequence}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return sequence_headers


def validation_status(forward: dict[str, str] | None, reciprocal: dict[str, str] | None, gene: str) -> tuple[str, str]:
    if forward is None:
        return "not_validated", "no_forward_hit"
    qcov = float(forward["qcovhsp"])
    scov = float(forward["scovhsp"])
    pident = float(forward["pident"])
    if qcov < 50 or scov < 25:
        return "weak_forward_support", "low_forward_coverage"
    if reciprocal is None:
        return "forward_only", "no_reciprocal_hit"
    reciprocal_gene = query_gene(reciprocal["sseqid"])
    if reciprocal_gene != gene:
        return "not_reciprocal", f"reciprocal_best_is_{reciprocal_gene}"
    if pident >= 30 and qcov >= 60 and scov >= 30:
        return "reciprocal_supported", "same_gene_reciprocal_top_hit"
    return "reciprocal_weak", "same_gene_reciprocal_top_hit_but_weak_forward_metrics"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=pathlib.Path, default=DEFAULT_BATCH)
    parser.add_argument("--download-log", type=pathlib.Path, default=DEFAULT_DOWNLOAD_LOG)
    parser.add_argument("--reference-fasta", type=pathlib.Path, default=DEFAULT_REFERENCE_FASTA)
    parser.add_argument("--work-dir", type=pathlib.Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--forward-output", type=pathlib.Path, default=DEFAULT_FORWARD)
    parser.add_argument("--reciprocal-output", type=pathlib.Path, default=DEFAULT_RECIPROCAL)
    parser.add_argument("--validation-output", type=pathlib.Path, default=DEFAULT_VALIDATION)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    parser.add_argument("--diamond", type=pathlib.Path, default=DEFAULT_DIAMOND)
    args = parser.parse_args()

    batch_rows = read_tsv(args.batch)
    downloads = read_tsv(args.download_log)
    fasta_by_accession = {
        row["best_assembly_accession"]: pathlib.Path(row["protein_local_path"])
        for row in downloads
        if row["download_status"] in {"downloaded", "cached"} and row["protein_local_path"]
    }

    all_forward_rows = []
    best_forward: dict[tuple[str, str], dict[str, str]] = {}
    for accession, fasta in fasta_by_accession.items():
        db = args.work_dir / "target_dbs" / accession
        make_db(args.diamond, fasta, db)
        out = args.work_dir / "forward_hits" / f"{accession}.tsv"
        blastp(args.diamond, args.reference_fasta, db, out, max_target=5)
        rows = parse_diamond(out)
        for row in rows:
            row["best_assembly_accession"] = accession
            row["query_gene"] = query_gene(row["qseqid"])
            all_forward_rows.append(row)
        for gene, hit in best_by_gene(rows).items():
            best_forward[(accession, gene)] = hit

    candidate_fasta = args.work_dir / "candidate_target_hits.faa"
    write_candidate_fasta(candidate_fasta, batch_rows, best_forward, fasta_by_accession)
    reference_db = args.work_dir / "reference_db" / "week4_transposon_reference"
    make_db(args.diamond, args.reference_fasta, reference_db)
    reciprocal_raw = args.work_dir / "reciprocal_hits.tsv"
    blastp(args.diamond, candidate_fasta, reference_db, reciprocal_raw, max_target=3)
    reciprocal_rows = parse_diamond(reciprocal_raw)
    for row in reciprocal_rows:
        parts = row["qseqid"].split("|")
        row["best_assembly_accession"] = parts[0] if len(parts) > 0 else ""
        row["query_gene"] = parts[1] if len(parts) > 1 else ""
        row["target_protein_id"] = parts[2] if len(parts) > 2 else ""
        row["reciprocal_gene"] = query_gene(row["sseqid"])
    reciprocal_best: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in reciprocal_rows:
        key = (row["best_assembly_accession"], row["query_gene"], row["target_protein_id"])
        if key not in reciprocal_best or float(row["bitscore"]) > float(reciprocal_best[key]["bitscore"]):
            reciprocal_best[key] = row

    validation_rows = []
    for row in batch_rows:
        accession = row["best_assembly_accession"]
        gene = row["human_gene_symbol"]
        forward = best_forward.get((accession, gene))
        reciprocal = None
        if forward:
            reciprocal = reciprocal_best.get((accession, gene, forward["sseqid"]))
        status, reason = validation_status(forward, reciprocal, gene)
        validation_rows.append(
            {
                **row,
                "week4_diamond_status": status,
                "week4_diamond_reason": reason,
                "forward_query": forward["qseqid"] if forward else "",
                "forward_target_protein_id": forward["sseqid"] if forward else "",
                "forward_pident": forward["pident"] if forward else "",
                "forward_qcovhsp": forward["qcovhsp"] if forward else "",
                "forward_scovhsp": forward["scovhsp"] if forward else "",
                "forward_evalue": forward["evalue"] if forward else "",
                "forward_bitscore": forward["bitscore"] if forward else "",
                "forward_target_title": forward["stitle"] if forward else "",
                "reciprocal_top_query": reciprocal["qseqid"] if reciprocal else "",
                "reciprocal_top_reference": reciprocal["sseqid"] if reciprocal else "",
                "reciprocal_gene": reciprocal["reciprocal_gene"] if reciprocal else "",
                "reciprocal_pident": reciprocal["pident"] if reciprocal else "",
                "reciprocal_qcovhsp": reciprocal["qcovhsp"] if reciprocal else "",
                "reciprocal_scovhsp": reciprocal["scovhsp"] if reciprocal else "",
                "reciprocal_evalue": reciprocal["evalue"] if reciprocal else "",
                "reciprocal_bitscore": reciprocal["bitscore"] if reciprocal else "",
            }
        )

    forward_fields = [
        "best_assembly_accession",
        "query_gene",
        *OUTFMT[1:],
    ]
    reciprocal_fields = [
        "best_assembly_accession",
        "query_gene",
        "target_protein_id",
        "reciprocal_gene",
        *OUTFMT[1:],
    ]
    write_tsv(args.forward_output, all_forward_rows, forward_fields)
    write_tsv(args.reciprocal_output, reciprocal_rows, reciprocal_fields)
    write_tsv(args.validation_output, validation_rows, list(validation_rows[0].keys()))

    status_counts = defaultdict(int)
    gene_counts = defaultdict(int)
    for row in validation_rows:
        status_counts[row["week4_diamond_status"]] += 1
        gene_counts[(row["human_gene_symbol"], row["week4_diamond_status"])] += 1
    lines = [
        "# Week 4 Transposon Reciprocal DIAMOND Validation Report",
        "",
        f"Batch rows validated: {len(validation_rows)}",
        f"Target proteomes searched: {len(fasta_by_accession)}",
        f"Forward hit rows: {len(all_forward_rows)}",
        f"Reciprocal hit rows: {len(reciprocal_rows)}",
        "",
        "## Validation Status",
        "",
    ]
    for status, count in sorted(status_counts.items()):
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Gene x Status", ""])
    for (gene, status), count in sorted(gene_counts.items()):
        lines.append(f"- {gene} / {status}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This is a reciprocal-style sequence validation against local NCBI protein FASTA files. It supports orthology triage but does not replace curated OMA, OrthoDB, or Ensembl Compara calls.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.validation_output}, {args.forward_output}, {args.reciprocal_output}, and {args.report}")


if __name__ == "__main__":
    main()
