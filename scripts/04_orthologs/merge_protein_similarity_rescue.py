"""Merge protein-similarity rescue candidates into the primary candidate matrix."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import Counter


DEFAULT_MATRIX = pathlib.Path("data/processed/ortholog_matrix_primary_combined_candidates.tsv")
DEFAULT_RESCUE = pathlib.Path("data/processed/erythrura_protein_similarity_rescue.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/ortholog_matrix_primary_final_candidates.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/ortholog_matrix_primary_final_candidates_report.md")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def key(row: dict[str, str]) -> tuple[str, str]:
    return row["scientific_name"], row["human_gene_symbol"]


def merge_rows(matrix_rows: list[dict[str, str]], rescue_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rescue_hits = {
        key(row): row
        for row in rescue_rows
        if row["protein_rescue_status"] == "protein_similarity_candidate"
    }
    merged = []
    for row in matrix_rows:
        updated = dict(row)
        hit = rescue_hits.get(key(row))
        updated["final_candidate_status"] = row["combined_candidate_status"]
        updated["final_candidate_source"] = row["combined_candidate_source"]
        updated["final_candidate_confidence"] = row["combined_candidate_confidence"]
        if row["combined_candidate_status"] == "not_found" and hit:
            updated["ortholog_query_status"] = "rescued_from_protein_similarity"
            updated["ortholog_gene_id"] = hit["target_protein_accession"]
            updated["ortholog_gene_symbol"] = hit["human_gene_symbol"]
            updated["ortholog_status"] = "candidate_by_reference_bird_protein_similarity"
            updated["ortholog_source_database"] = "local protein FASTA kmer similarity"
            updated["ortholog_source_url"] = ""
            updated["ortholog_confidence"] = hit["protein_rescue_confidence"]
            updated["copy_number_estimate"] = "not_estimated"
            updated["notes"] = (
                f"reference={hit['reference_species']}:{hit['reference_protein_accession']}; "
                f"target={hit['target_protein_accession']}; "
                f"kmer_coverage={hit['query_kmer_coverage']}; "
                f"jaccard={hit['kmer_jaccard']}; "
                f"length_ratio={hit['length_ratio']}"
            )
            updated["final_candidate_status"] = "protein_similarity_candidate"
            updated["final_candidate_source"] = "Taeniopygia reference protein vs Erythrura protein FASTA"
            updated["final_candidate_confidence"] = hit["protein_rescue_confidence"]
        merged.append(updated)
    return merged


def write_report(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    status_counts = Counter(row["final_candidate_status"] for row in rows)
    total_found = sum(count for status, count in status_counts.items() if status != "not_found")
    coverage = total_found / len(rows) if rows else 0.0
    ery = [row for row in rows if row["scientific_name"] == "Erythrura gouldiae"]
    ery_found = [row for row in ery if row["final_candidate_status"] != "not_found"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Primary Final Candidate Matrix Report",
                "",
                f"Input rows: {len(rows)}",
                f"Final candidate hits: {total_found}",
                f"Final candidate coverage: {coverage:.3f}",
                f"Erythrura gouldiae candidate hits: {len(ery_found)}/{len(ery)}",
                "",
                "## Candidate Status",
                *[f"- {status}: {count}" for status, count in sorted(status_counts.items())],
                "",
                "## Interpretation",
                "This final feasibility matrix combines NCBI Gene candidates, GFF symbol rescue, and low-confidence protein-similarity rescue for the poorly named Erythrura gouldiae assembly. Protein rescue rows should be validated with BLAST/DIAMOND or a formal orthology database before biological interpretation.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=pathlib.Path, default=DEFAULT_MATRIX)
    parser.add_argument("--rescue", type=pathlib.Path, default=DEFAULT_RESCUE)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    matrix_rows = read_tsv(args.matrix)
    rescue_rows = read_tsv(args.rescue)
    merged = merge_rows(matrix_rows, rescue_rows)
    fields = list(matrix_rows[0].keys()) + [
        "final_candidate_status",
        "final_candidate_source",
        "final_candidate_confidence",
    ]
    write_tsv(args.output, merged, fields)
    write_report(args.report, merged)
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
