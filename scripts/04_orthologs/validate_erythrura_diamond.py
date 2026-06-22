"""Validate Erythrura protein-similarity rescue candidates with DIAMOND."""

from __future__ import annotations

import argparse
import csv
import pathlib
import subprocess
from collections import defaultdict, Counter


DEFAULT_DIAMOND = pathlib.Path("env/tools/diamond-2.2.1/diamond.exe")
DEFAULT_REFERENCE_FASTA = pathlib.Path("data/interim/reference_bird_proteins_taeniopygia.faa")
DEFAULT_TARGET_FASTA = pathlib.Path(
    "data/raw/annotation_rescue/GCA_003676055.1/GCA_003676055.1_GouldianFinch_protein.faa.gz"
)
DEFAULT_RESCUE = pathlib.Path("data/processed/erythrura_protein_similarity_rescue.tsv")
DEFAULT_DB = pathlib.Path("data/interim/diamond/erythrura_gouldiae_protein_db")
DEFAULT_HITS = pathlib.Path("results/tables/erythrura_diamond_hits.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/erythrura_diamond_validation.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/erythrura_diamond_validation_report.md")


HIT_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
]


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_command(command: list[str]) -> None:
    subprocess.run(command, check=True)


def make_db(diamond: pathlib.Path, target_fasta: pathlib.Path, db: pathlib.Path) -> None:
    db.parent.mkdir(parents=True, exist_ok=True)
    if db.with_suffix(".dmnd").exists():
        return
    run_command([str(diamond), "makedb", "--in", str(target_fasta), "-d", str(db)])


def run_diamond(
    diamond: pathlib.Path,
    db: pathlib.Path,
    reference_fasta: pathlib.Path,
    hits_output: pathlib.Path,
    threads: int,
) -> None:
    hits_output.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            str(diamond),
            "blastp",
            "-d",
            str(db),
            "-q",
            str(reference_fasta),
            "-o",
            str(hits_output),
            "-f",
            "6",
            *HIT_FIELDS,
            "--max-target-seqs",
            "10",
            "--evalue",
            "1e-5",
            "--threads",
            str(threads),
        ]
    )


def parse_hits(path: pathlib.Path) -> list[dict[str, str]]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            values = line.rstrip("\n").split("\t")
            if len(values) != len(HIT_FIELDS):
                continue
            rows.append(dict(zip(HIT_FIELDS, values)))
    return rows


def query_gene(qseqid: str) -> str:
    return qseqid.split("|", 1)[0]


def num(row: dict[str, str], field: str) -> float:
    return float(row[field])


def hit_metrics(hit: dict[str, str]) -> dict[str, float]:
    qlen = num(hit, "qlen")
    slen = num(hit, "slen")
    length = num(hit, "length")
    return {
        "pident": num(hit, "pident"),
        "evalue": num(hit, "evalue"),
        "bitscore": num(hit, "bitscore"),
        "alignment_length": length,
        "query_coverage": length / qlen if qlen else 0.0,
        "subject_coverage": length / slen if slen else 0.0,
    }


def validation_status(rank: int | None, metrics: dict[str, float]) -> tuple[str, str]:
    if rank is None:
        return "not_validated", "expected_target_not_in_top10"
    if (
        rank == 1
        and metrics["evalue"] <= 1e-20
        and metrics["pident"] >= 40.0
        and metrics["query_coverage"] >= 0.60
    ):
        return "diamond_validated_high", "expected_target_is_top_hit"
    if (
        rank <= 5
        and metrics["evalue"] <= 1e-10
        and metrics["pident"] >= 30.0
        and metrics["query_coverage"] >= 0.50
    ):
        return "diamond_validated_medium", "expected_target_in_top5"
    return "not_validated", "diamond_threshold_not_met"


def validate(rescue_rows: list[dict[str, str]], hit_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    hits_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for hit in hit_rows:
        hits_by_gene[query_gene(hit["qseqid"])].append(hit)
    for hits in hits_by_gene.values():
        hits.sort(key=lambda row: float(row["bitscore"]), reverse=True)

    output = []
    for row in rescue_rows:
        if row["protein_rescue_status"] != "protein_similarity_candidate":
            continue
        gene = row["human_gene_symbol"]
        expected = row["target_protein_accession"]
        hits = hits_by_gene.get(gene, [])
        top_hit = hits[0] if hits else {}
        expected_hit = None
        expected_rank = None
        for idx, hit in enumerate(hits, start=1):
            if hit["sseqid"] == expected:
                expected_hit = hit
                expected_rank = idx
                break
        metrics = hit_metrics(expected_hit) if expected_hit else {
            "pident": 0.0,
            "evalue": 1.0,
            "bitscore": 0.0,
            "alignment_length": 0.0,
            "query_coverage": 0.0,
            "subject_coverage": 0.0,
        }
        status, reason = validation_status(expected_rank, metrics)
        output.append(
            {
                "scientific_name": row["scientific_name"],
                "human_gene_symbol": gene,
                "maintenance_module": row["maintenance_module"],
                "reference_species": row["reference_species"],
                "reference_protein_accession": row["reference_protein_accession"],
                "expected_target_protein_accession": expected,
                "diamond_top_target": top_hit.get("sseqid", ""),
                "expected_target_rank": str(expected_rank or ""),
                "diamond_pident": f"{metrics['pident']:.3f}",
                "diamond_alignment_length": f"{metrics['alignment_length']:.0f}",
                "diamond_query_coverage": f"{metrics['query_coverage']:.6f}",
                "diamond_subject_coverage": f"{metrics['subject_coverage']:.6f}",
                "diamond_evalue": f"{metrics['evalue']:.3g}",
                "diamond_bitscore": f"{metrics['bitscore']:.3f}",
                "diamond_validation_status": status,
                "diamond_validation_reason": reason,
                "prior_protein_rescue_confidence": row["protein_rescue_confidence"],
            }
        )
    output.sort(key=lambda item: (item["diamond_validation_status"], item["human_gene_symbol"]))
    return output


def write_report(path: pathlib.Path, rows: list[dict[str, str]], hits: list[dict[str, str]]) -> None:
    status_counts = Counter(row["diamond_validation_status"] for row in rows)
    validated = [
        row for row in rows if row["diamond_validation_status"].startswith("diamond_validated")
    ]
    high = [row for row in rows if row["diamond_validation_status"] == "diamond_validated_high"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Erythrura DIAMOND Validation Report",
                "",
                f"DIAMOND hit rows: {len(hits)}",
                f"Protein rescue candidates tested: {len(rows)}",
                f"DIAMOND validated candidates: {len(validated)}",
                f"High-confidence validated candidates: {len(high)}",
                "",
                "## Validation Status",
                *[f"- {status}: {count}" for status, count in sorted(status_counts.items())],
                "",
                "## Interpretation",
                "DIAMOND validation checks whether the k-mer-predicted Erythrura target protein is recovered among the best local protein alignments from Taeniopygia reference proteins. Validated rows are stronger feasibility candidates, but manuscript-grade orthology should still be cross-checked with an orthology database or reciprocal search.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--diamond", type=pathlib.Path, default=DEFAULT_DIAMOND)
    parser.add_argument("--reference-fasta", type=pathlib.Path, default=DEFAULT_REFERENCE_FASTA)
    parser.add_argument("--target-fasta", type=pathlib.Path, default=DEFAULT_TARGET_FASTA)
    parser.add_argument("--rescue", type=pathlib.Path, default=DEFAULT_RESCUE)
    parser.add_argument("--db", type=pathlib.Path, default=DEFAULT_DB)
    parser.add_argument("--hits", type=pathlib.Path, default=DEFAULT_HITS)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    make_db(args.diamond, args.target_fasta, args.db)
    run_diamond(args.diamond, args.db, args.reference_fasta, args.hits, args.threads)
    hits = parse_hits(args.hits)
    rows = validate(read_tsv(args.rescue), hits)
    write_tsv(args.output, rows, list(rows[0].keys()) if rows else [])
    write_report(args.report, rows, hits)
    print(f"Wrote {args.hits}, {args.output}, and {args.report}")


if __name__ == "__main__":
    main()
