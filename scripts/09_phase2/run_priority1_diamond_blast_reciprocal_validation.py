"""Run DIAMOND and BLASTP reciprocal validation for priority 1 protein candidates."""

from __future__ import annotations

import argparse
import csv
import pathlib
import subprocess
from collections import Counter

import pandas as pd


HIT_FIELDS = [
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
]

BLAST_HIT_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qcovhsp",
    "evalue",
    "bitscore",
]


def run_command(command: list[str]) -> None:
    subprocess.run(command, check=True, text=True, capture_output=True)


def parse_fasta_ids(path: pathlib.Path) -> dict[str, str]:
    headers = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            seqid = line[1:].split()[0]
            headers[seqid] = line[1:]
    return headers


def gene_from_seqid(seqid: str) -> str:
    return str(seqid).split("|", 1)[0]


def species_from_candidate(seqid: str) -> str:
    parts = str(seqid).split("|")
    return parts[1].replace("_", " ") if len(parts) > 1 else ""


def protein_accession_from_candidate(seqid: str) -> str:
    parts = str(seqid).split("|")
    if len(parts) < 5:
        return ""
    return parts[4].split()[0]


def make_diamond_db(diamond: str, fasta: pathlib.Path, db_prefix: pathlib.Path) -> None:
    dmnd = pathlib.Path(str(db_prefix) + ".dmnd")
    if dmnd.exists() and dmnd.stat().st_size > 0:
        return
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    run_command([diamond, "makedb", "--in", str(fasta), "-d", str(db_prefix)])


def run_diamond(diamond: str, query: pathlib.Path, db_prefix: pathlib.Path, out: pathlib.Path, threads: int) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            diamond,
            "blastp",
            "-q",
            str(query),
            "-d",
            str(db_prefix),
            "-o",
            str(out),
            "-f",
            "6",
            *HIT_FIELDS,
            "--max-target-seqs",
            "10",
            "--evalue",
            "1e-5",
            "--threads",
            str(threads),
            "--quiet",
        ]
    )


def make_blast_db(makeblastdb: str, fasta: pathlib.Path, db_prefix: pathlib.Path) -> None:
    pin = pathlib.Path(str(db_prefix) + ".pin")
    if pin.exists() and pin.stat().st_size > 0:
        return
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    run_command([makeblastdb, "-in", str(fasta), "-dbtype", "prot", "-out", str(db_prefix)])


def run_blastp(blastp: str, query: pathlib.Path, db_prefix: pathlib.Path, out: pathlib.Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            blastp,
            "-query",
            str(query),
            "-db",
            str(db_prefix),
            "-out",
            str(out),
            "-outfmt",
            "6 " + " ".join(BLAST_HIT_FIELDS),
            "-max_target_seqs",
            "10",
            "-evalue",
            "1e-5",
        ]
    )


def parse_hits(path: pathlib.Path, tool: str, direction: str) -> pd.DataFrame:
    fields = BLAST_HIT_FIELDS if tool == "blastp" else HIT_FIELDS
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=HIT_FIELDS + ["tool", "direction"])
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for parts in reader:
            if len(parts) < len(fields):
                continue
            row = dict(zip(fields, parts, strict=False))
            row["tool"] = tool
            row["direction"] = direction
            rows.append(row)
    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=HIT_FIELDS + ["tool", "direction"])
    for col in ["pident", "length", "qlen", "slen", "qcovhsp", "scovhsp", "evalue", "bitscore"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if "scovhsp" not in df.columns:
        df["scovhsp"] = (df["length"] / df["slen"] * 100).where(df["slen"] > 0, 0)
    return df[HIT_FIELDS + ["tool", "direction"]]


def best_same_gene_status(forward_hits: pd.DataFrame, reciprocal_hits: pd.DataFrame) -> list[dict[str, object]]:
    rows = []
    if forward_hits.empty:
        return rows
    for qseqid, group in forward_hits.sort_values("bitscore", ascending=False).groupby("qseqid"):
        gene = gene_from_seqid(qseqid)
        same_gene = group[group["sseqid"].map(gene_from_seqid) == gene].copy()
        top = group.iloc[0]
        best_same = same_gene.iloc[0] if not same_gene.empty else None
        reciprocal_same_gene_top = False
        reciprocal_best_gene = ""
        if best_same is not None:
            ref_id = best_same["sseqid"]
            reciprocal_group = reciprocal_hits[reciprocal_hits["qseqid"] == ref_id].sort_values("bitscore", ascending=False)
            if not reciprocal_group.empty:
                reciprocal_top = reciprocal_group.iloc[0]
                reciprocal_best_gene = gene_from_seqid(reciprocal_top["sseqid"])
                reciprocal_same_gene_top = reciprocal_best_gene == gene
        if best_same is None:
            call = "no_same_gene_reference_hit"
            reason = "candidate_best_hits_do_not_include_same_gene_reference"
            metrics = top
            same_gene_rank = ""
        else:
            same_gene_rank = int(group.index.get_loc(best_same.name) + 1)
            metrics = best_same
            if (
                reciprocal_same_gene_top
                and float(best_same["pident"]) >= 35
                and float(best_same["qcovhsp"]) >= 50
                and float(best_same["scovhsp"]) >= 50
            ):
                call = "reciprocal_same_gene_supported"
                reason = "same_gene_reference_hit_and_reciprocal_top_returns_same_gene"
            elif (
                float(best_same["pident"]) >= 30
                and float(best_same["qcovhsp"]) >= 40
                and float(best_same["scovhsp"]) >= 40
            ):
                call = "same_gene_forward_supported_reciprocal_weak"
                reason = "same_gene_reference_hit_but_reciprocal_not_top_or_metrics_weaker"
            else:
                call = "weak_same_gene_support"
                reason = "same_gene_reference_hit_below_threshold"
        rows.append(
            {
                "qseqid": qseqid,
                "human_gene_symbol": gene,
                "scientific_name": species_from_candidate(qseqid),
                "candidate_protein_accession": protein_accession_from_candidate(qseqid),
                "same_gene_hit_rank": same_gene_rank,
                "top_reference_hit": top["sseqid"],
                "top_reference_gene": gene_from_seqid(top["sseqid"]),
                "best_same_gene_reference_hit": "" if best_same is None else best_same["sseqid"],
                "reciprocal_best_gene": reciprocal_best_gene,
                "pident": float(metrics["pident"]),
                "alignment_length": int(metrics["length"]),
                "qcovhsp": float(metrics["qcovhsp"]),
                "scovhsp": float(metrics["scovhsp"]),
                "evalue": float(metrics["evalue"]),
                "bitscore": float(metrics["bitscore"]),
                "validation_call": call,
                "validation_reason": reason,
            }
        )
    return rows


def gene_decision(row: pd.Series) -> str:
    supported = int(row["reciprocal_supported_species"])
    forward = int(row["forward_supported_species"])
    species = int(row["species_rows"])
    if supported >= max(4, species - 1):
        if str(row["human_gene_symbol"]).startswith("TDRD"):
            return "domain_paralog_check_required_after_reciprocal_support"
        return "strict_upgrade_candidate_diamond_blast_supported"
    if supported + forward >= max(4, species - 1):
        return "sequence_supported_reciprocal_incomplete_manual_review"
    if supported + forward >= 3:
        return "mixed_sequence_support_manual_review"
    return "hold_not_strict_not_absence"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--human-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--work-dir", type=pathlib.Path, required=True)
    parser.add_argument("--row-output", type=pathlib.Path, required=True)
    parser.add_argument("--gene-summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--all-hits-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--diamond", default="diamond")
    parser.add_argument("--blastp", default="blastp")
    parser.add_argument("--makeblastdb", default="makeblastdb")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    args.work_dir.mkdir(parents=True, exist_ok=True)
    human_diamond_db = args.work_dir / "diamond_human_refs"
    candidate_diamond_db = args.work_dir / "diamond_candidates"
    make_diamond_db(args.diamond, args.human_fasta, human_diamond_db)
    make_diamond_db(args.diamond, args.candidate_fasta, candidate_diamond_db)
    diamond_forward = args.work_dir / "diamond_candidate_to_human.tsv"
    diamond_reciprocal = args.work_dir / "diamond_human_to_candidate.tsv"
    run_diamond(args.diamond, args.candidate_fasta, human_diamond_db, diamond_forward, args.threads)
    run_diamond(args.diamond, args.human_fasta, candidate_diamond_db, diamond_reciprocal, args.threads)

    human_blast_db = args.work_dir / "blast_human_refs"
    candidate_blast_db = args.work_dir / "blast_candidates"
    make_blast_db(args.makeblastdb, args.human_fasta, human_blast_db)
    make_blast_db(args.makeblastdb, args.candidate_fasta, candidate_blast_db)
    blast_forward = args.work_dir / "blast_candidate_to_human.tsv"
    blast_reciprocal = args.work_dir / "blast_human_to_candidate.tsv"
    run_blastp(args.blastp, args.candidate_fasta, human_blast_db, blast_forward)
    run_blastp(args.blastp, args.human_fasta, candidate_blast_db, blast_reciprocal)

    hit_tables = [
        parse_hits(diamond_forward, "diamond", "candidate_to_human"),
        parse_hits(diamond_reciprocal, "diamond", "human_to_candidate"),
        parse_hits(blast_forward, "blastp", "candidate_to_human"),
        parse_hits(blast_reciprocal, "blastp", "human_to_candidate"),
    ]
    all_hits = pd.concat(hit_tables, ignore_index=True)
    args.all_hits_output.parent.mkdir(parents=True, exist_ok=True)
    all_hits.to_csv(args.all_hits_output, sep="\t", index=False)

    row_tables = []
    for tool in ["diamond", "blastp"]:
        forward = all_hits[(all_hits["tool"] == tool) & (all_hits["direction"] == "candidate_to_human")]
        reciprocal = all_hits[(all_hits["tool"] == tool) & (all_hits["direction"] == "human_to_candidate")]
        rows = pd.DataFrame(best_same_gene_status(forward, reciprocal))
        rows["tool"] = tool
        row_tables.append(rows)
    row_validation = pd.concat(row_tables, ignore_index=True)
    args.row_output.parent.mkdir(parents=True, exist_ok=True)
    row_validation.to_csv(args.row_output, sep="\t", index=False)

    best_per_tool_species = (
        row_validation.sort_values(
            ["tool", "human_gene_symbol", "scientific_name", "validation_call", "bitscore"],
            ascending=[True, True, True, True, False],
        )
        .groupby(["tool", "human_gene_symbol", "scientific_name"], as_index=False)
        .first()
    )
    best_per_species = (
        row_validation.assign(
            support_rank=row_validation["validation_call"].map(
                {
                    "reciprocal_same_gene_supported": 3,
                    "same_gene_forward_supported_reciprocal_weak": 2,
                    "weak_same_gene_support": 1,
                    "no_same_gene_reference_hit": 0,
                }
            )
        )
        .sort_values(["human_gene_symbol", "scientific_name", "support_rank", "bitscore"], ascending=[True, True, False, False])
        .groupby(["human_gene_symbol", "scientific_name"], as_index=False)
        .first()
    )
    gene_summary = (
        best_per_species.groupby("human_gene_symbol", as_index=False)
        .agg(
            species_rows=("scientific_name", "nunique"),
            reciprocal_supported_species=("validation_call", lambda s: int((s == "reciprocal_same_gene_supported").sum())),
            forward_supported_species=("validation_call", lambda s: int((s == "same_gene_forward_supported_reciprocal_weak").sum())),
            weak_or_no_support_species=("validation_call", lambda s: int((~s.isin(["reciprocal_same_gene_supported", "same_gene_forward_supported_reciprocal_weak"])).sum())),
            median_pident=("pident", "median"),
            median_qcovhsp=("qcovhsp", "median"),
            median_scovhsp=("scovhsp", "median"),
            best_tools=("tool", lambda s: ";".join(sorted(set(map(str, s))))),
        )
    )
    gene_summary["diamond_blast_upgrade_decision"] = gene_summary.apply(gene_decision, axis=1)
    args.gene_summary_output.parent.mkdir(parents=True, exist_ok=True)
    gene_summary.to_csv(args.gene_summary_output, sep="\t", index=False)

    call_counts = Counter(row_validation["validation_call"])
    decision_counts = Counter(gene_summary["diamond_blast_upgrade_decision"])
    lines = [
        "# Phase 2 Priority 1 DIAMOND/BLAST Reciprocal Validation Report",
        "",
        "## Summary",
        "",
        f"Candidate protein records: {len(parse_fasta_ids(args.candidate_fasta))}",
        f"Human reference records: {len(parse_fasta_ids(args.human_fasta))}",
        f"All hit rows: {len(all_hits)}",
        f"Row validation records: {len(row_validation)}",
        f"Genes assessed: {gene_summary['human_gene_symbol'].nunique()}",
        "",
        "## Row-Level Calls",
        "",
    ]
    for call, count in sorted(call_counts.items()):
        lines.append(f"- {call}: {count}")
    lines.extend(["", "## Gene-Level Decisions", ""])
    for decision, count in sorted(decision_counts.items()):
        lines.append(f"- {decision}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This validation uses both DIAMOND and BLASTP candidate-to-human and human-to-candidate searches. Non-TDRD genes with reciprocal same-gene support across at least four species can be upgraded as sequence-supported strict candidates. TDRD-family genes remain domain/paralog-check pending even when reciprocal similarity is strong.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.row_output}")


if __name__ == "__main__":
    main()
