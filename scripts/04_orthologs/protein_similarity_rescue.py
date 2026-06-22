"""Protein k-mer similarity rescue for poorly named annotations."""

from __future__ import annotations

import argparse
import csv
import gzip
import pathlib
import re


DEFAULT_REFERENCE = pathlib.Path("data/interim/reference_bird_proteins_taeniopygia.tsv")
DEFAULT_TARGET_FASTA = pathlib.Path("data/raw/annotation_rescue/GCA_003676055.1/GCA_003676055.1_GouldianFinch_protein.faa.gz")
DEFAULT_SPECIES = "Erythrura gouldiae"
DEFAULT_OUTPUT = pathlib.Path("data/processed/erythrura_protein_similarity_rescue.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/erythrura_protein_similarity_rescue_report.md")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    fields = list(rows[0].keys()) if rows else []
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_fasta(path: pathlib.Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    records = []
    header = ""
    seq_parts: list[str] = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
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


def accession(header: str) -> str:
    return header.split()[0]


def kmer_set(seq: str, k: int) -> set[str]:
    seq = re.sub("[^A-Z]", "", seq.upper())
    if len(seq) < k:
        return set()
    return {seq[i : i + k] for i in range(len(seq) - k + 1)}


def score_pair(query_kmers: set[str], target_kmers: set[str], query_len: int, target_len: int) -> tuple[float, float, float]:
    if not query_kmers or not target_kmers:
        return 0.0, 0.0, 0.0
    shared = len(query_kmers & target_kmers)
    query_coverage = shared / len(query_kmers)
    jaccard = shared / len(query_kmers | target_kmers)
    length_ratio = min(query_len, target_len) / max(query_len, target_len) if max(query_len, target_len) else 0.0
    return query_coverage, jaccard, length_ratio


def confidence(query_coverage: float, jaccard: float, length_ratio: float) -> str:
    if query_coverage >= 0.55 and jaccard >= 0.35 and length_ratio >= 0.60:
        return "medium"
    if query_coverage >= 0.35 and jaccard >= 0.20 and length_ratio >= 0.45:
        return "low"
    return "below_threshold"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=pathlib.Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--target-fasta", type=pathlib.Path, default=DEFAULT_TARGET_FASTA)
    parser.add_argument("--species", default=DEFAULT_SPECIES)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    parser.add_argument("--k", type=int, default=5)
    args = parser.parse_args()

    references = [row for row in read_tsv(args.reference) if row.get("fetch_status") == "ok" and row.get("sequence")]
    targets = parse_fasta(args.target_fasta)
    target_kmers = []
    for target in targets:
        seq = target["sequence"]
        target_kmers.append((target, kmer_set(seq, args.k), len(seq)))

    rows = []
    for ref in references:
        query_seq = ref["sequence"]
        query_kmers = kmer_set(query_seq, args.k)
        best = None
        for target, kmers, target_len in target_kmers:
            query_cov, jaccard, length_ratio = score_pair(query_kmers, kmers, len(query_seq), target_len)
            score = (query_cov, jaccard, length_ratio)
            if best is None or score > best["score_tuple"]:
                best = {
                    "target": target,
                    "query_coverage": query_cov,
                    "jaccard": jaccard,
                    "length_ratio": length_ratio,
                    "score_tuple": score,
                }
        if best is None:
            continue
        conf = confidence(best["query_coverage"], best["jaccard"], best["length_ratio"])
        rows.append(
            {
                "scientific_name": args.species,
                "human_gene_symbol": ref["human_gene_symbol"],
                "maintenance_module": ref["maintenance_module"],
                "reference_species": ref["reference_species"],
                "reference_gene_id": ref["reference_gene_id"],
                "reference_protein_accession": ref["reference_protein_accession"],
                "target_protein_accession": accession(best["target"]["header"]),
                "target_header": best["target"]["header"],
                "reference_length": str(len(query_seq)),
                "target_length": str(len(best["target"]["sequence"])),
                "k": str(args.k),
                "query_kmer_coverage": f"{best['query_coverage']:.6f}",
                "kmer_jaccard": f"{best['jaccard']:.6f}",
                "length_ratio": f"{best['length_ratio']:.6f}",
                "protein_rescue_status": "protein_similarity_candidate" if conf != "below_threshold" else "below_threshold",
                "protein_rescue_confidence": conf,
            }
        )

    write_tsv(args.output, rows)
    accepted = [row for row in rows if row["protein_rescue_status"] == "protein_similarity_candidate"]
    medium = [row for row in accepted if row["protein_rescue_confidence"] == "medium"]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        "\n".join(
            [
                "# Erythrura Protein Similarity Rescue Report",
                "",
                f"Reference proteins searched: {len(references)}",
                f"Target proteins: {len(targets)}",
                f"Accepted protein-similarity candidates: {len(accepted)}",
                f"Medium-confidence candidates: {len(medium)}",
                "",
                "## Interpretation",
                "These are k-mer protein similarity candidates against a poorly named GenBank protein set. They are useful to rescue zero symbol coverage, but they should be validated with BLAST/DIAMOND or an orthology database before final mechanism scoring.",
            ]
        ),
        encoding="utf-8",
    )
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
