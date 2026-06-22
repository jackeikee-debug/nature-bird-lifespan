"""Run lightweight human-reference sequence confirmation for priority 1 proteins."""

from __future__ import annotations

import argparse
import pathlib
import re

import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices


def parse_fasta(path: pathlib.Path) -> list[dict[str, object]]:
    records = []
    header = None
    seq_parts: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append(record_from_header(header, "".join(seq_parts)))
            header = line[1:].strip()
            seq_parts = []
        elif line.strip():
            seq_parts.append(line.strip())
    if header is not None:
        records.append(record_from_header(header, "".join(seq_parts)))
    return records


def record_from_header(header: str, sequence: str) -> dict[str, object]:
    parts = header.split("|")
    gene = parts[0] if parts else ""
    organism = parts[1] if len(parts) > 1 else ""
    rank = ""
    for part in parts:
        if part.startswith("rank:"):
            rank = part.split(":", 1)[1]
    return {
        "human_gene_symbol": gene,
        "organism": organism,
        "reference_rank": rank,
        "reference_header": header,
        "reference_sequence": sequence,
        "reference_length": len(sequence),
    }


def aligned_identity(alignment, query: str, target: str) -> tuple[int, int, int, int]:
    target_blocks, query_blocks = alignment.aligned
    matches = 0
    aligned_positions = 0
    query_covered = set()
    target_covered = set()
    for (t_start, t_end), (q_start, q_end) in zip(target_blocks, query_blocks):
        block_len = min(t_end - t_start, q_end - q_start)
        for i in range(block_len):
            if target[t_start + i] == query[q_start + i]:
                matches += 1
        aligned_positions += block_len
        query_covered.update(range(q_start, q_end))
        target_covered.update(range(t_start, t_end))
    return matches, aligned_positions, len(query_covered), len(target_covered)


def confirmation_call(row: pd.Series) -> str:
    if row["protein_fetch_status"] != "protein_sequence_found":
        return "no_candidate_sequence"
    if row["best_identity"] >= 0.55 and row["candidate_coverage"] >= 0.7 and row["human_reference_coverage"] >= 0.7:
        if row["sequence_confirmation_route"] == "domain_architecture_plus_reciprocal_sequence":
            return "sequence_supported_domain_check_pending"
        return "sequence_supported"
    if row["best_identity"] >= 0.35 and row["best_aligned_positions"] >= 100:
        return "partial_sequence_support_manual_review"
    return "sequence_not_supported_manual_review"


def gene_decision(row: pd.Series) -> str:
    if row["sequence_supported_or_pending_species"] >= max(4, row["species_with_sequences"] - 1):
        if row["domain_pending_species"] > 0:
            return "domain_architecture_check_required_before_strict_upgrade"
        return "strict_upgrade_candidate_sequence_supported"
    if row["sequence_supported_or_pending_species"] >= 3:
        return "mixed_sequence_support_expand_or_manual_review"
    return "hold_not_strict_not_absence"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-metadata", type=pathlib.Path, required=True)
    parser.add_argument("--human-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--row-output", type=pathlib.Path, required=True)
    parser.add_argument("--gene-summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    metadata = pd.read_csv(args.candidate_metadata, sep="\t")
    refs = pd.DataFrame(parse_fasta(args.human_fasta))
    refs = refs.sort_values(["human_gene_symbol", "reference_length"], ascending=[True, False])
    refs_by_gene = {
        gene: group.head(3).to_dict("records")
        for gene, group in refs.groupby("human_gene_symbol")
    }

    aligner = PairwiseAligner(mode="local")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    rows = []
    seq_rows = metadata[metadata["protein_fetch_status"] == "protein_sequence_found"].copy()
    for _, row in seq_rows.iterrows():
        gene = row["human_gene_symbol"]
        query = str(row["sequence"])
        best = None
        for ref in refs_by_gene.get(gene, []):
            target = str(ref["reference_sequence"])
            if not query or not target:
                continue
            alignment = aligner.align(target, query)[0]
            matches, aligned_positions, query_cov, target_cov = aligned_identity(alignment, query, target)
            result = row.to_dict()
            result.update(
                {
                    "reference_rank": ref["reference_rank"],
                    "reference_header": ref["reference_header"],
                    "reference_length": ref["reference_length"],
                    "alignment_score": float(alignment.score),
                    "best_matches": matches,
                    "best_aligned_positions": aligned_positions,
                    "best_identity": round(matches / aligned_positions, 4) if aligned_positions else 0,
                    "candidate_coverage": round(query_cov / len(query), 4) if query else 0,
                    "human_reference_coverage": round(target_cov / len(target), 4) if target else 0,
                }
            )
            if best is None or result["best_identity"] * result["candidate_coverage"] > best["best_identity"] * best["candidate_coverage"]:
                best = result
        if best is not None:
            rows.append(best)

    no_seq = metadata[metadata["protein_fetch_status"] != "protein_sequence_found"].copy()
    for _, row in no_seq.iterrows():
        result = row.to_dict()
        result.update(
            {
                "reference_rank": "",
                "reference_header": "",
                "reference_length": 0,
                "alignment_score": 0,
                "best_matches": 0,
                "best_aligned_positions": 0,
                "best_identity": 0,
                "candidate_coverage": 0,
                "human_reference_coverage": 0,
            }
        )
        rows.append(result)

    out = pd.DataFrame(rows)
    out["sequence_confirmation_call"] = out.apply(confirmation_call, axis=1)
    out = out.sort_values(
        ["human_gene_symbol", "scientific_name", "protein_rank_for_gene_species"]
    )
    args.row_output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.row_output, sep="\t", index=False)

    best_per_species = (
        out.sort_values(
            ["human_gene_symbol", "scientific_name", "best_identity", "candidate_coverage"],
            ascending=[True, True, False, False],
        )
        .groupby(["human_gene_symbol", "scientific_name"], as_index=False)
        .first()
    )
    gene_summary = (
        best_per_species.groupby("human_gene_symbol", as_index=False)
        .agg(
            species_rows=("scientific_name", "nunique"),
            species_with_sequences=("protein_fetch_status", lambda s: int((s == "protein_sequence_found").sum())),
            sequence_supported_species=("sequence_confirmation_call", lambda s: int((s == "sequence_supported").sum())),
            domain_pending_species=("sequence_confirmation_call", lambda s: int((s == "sequence_supported_domain_check_pending").sum())),
            partial_support_species=("sequence_confirmation_call", lambda s: int((s == "partial_sequence_support_manual_review").sum())),
            no_sequence_species=("sequence_confirmation_call", lambda s: int((s == "no_candidate_sequence").sum())),
            median_identity=("best_identity", "median"),
            min_identity=("best_identity", "min"),
            median_candidate_coverage=("candidate_coverage", "median"),
            min_candidate_coverage=("candidate_coverage", "min"),
        )
    )
    gene_summary["sequence_supported_or_pending_species"] = (
        gene_summary["sequence_supported_species"] + gene_summary["domain_pending_species"]
    )
    gene_summary["sequence_upgrade_decision"] = gene_summary.apply(gene_decision, axis=1)
    args.gene_summary_output.parent.mkdir(parents=True, exist_ok=True)
    gene_summary.to_csv(args.gene_summary_output, sep="\t", index=False)

    row_counts = out["sequence_confirmation_call"].value_counts().sort_index()
    decision_counts = gene_summary["sequence_upgrade_decision"].value_counts().sort_index()
    lines = [
        "# Phase 2 Priority 1 Sequence Similarity Confirmation Report",
        "",
        "## Summary",
        "",
        f"Candidate protein rows assessed: {len(out)}",
        f"Genes assessed: {out['human_gene_symbol'].nunique()}",
        "",
        "## Row-Level Calls",
        "",
    ]
    for call, count in row_counts.items():
        lines.append(f"- {call}: {count}")
    lines.extend(["", "## Gene-Level Decisions", ""])
    for decision, count in decision_counts.items():
        lines.append(f"- {decision}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This is a lightweight sequence-similarity screen against human reference proteins. Non-TDRD genes with broad sequence support can become strict-upgrade candidates; TDRD-family genes remain domain/paralog-check pending even when sequence similarity is strong.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.row_output}")


if __name__ == "__main__":
    main()
