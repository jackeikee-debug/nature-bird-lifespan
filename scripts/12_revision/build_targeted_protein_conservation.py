#!/usr/bin/env python3
"""Build a targeted cross-species protein-conservation evidence layer."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import zipfile
from pathlib import Path

import pandas as pd


FOCAL_GENES = [
    "DNMT1",
    "DNMT3A",
    "DNMT3B",
    "HELLS",
    "MBD2",
    "MBD3",
    "MORC3",
    "SAMHD1",
    "SETDB2",
    "UHRF1",
]

# NCBI Gene identifiers for the human reference anchors.
HUMAN_GENE_IDS = {
    "DNMT1": "1786",
    "DNMT3A": "1788",
    "DNMT3B": "1789",
    "HELLS": "3070",
    "MBD2": "8932",
    "MBD3": "53615",
    "MORC3": "23515",
    "SAMHD1": "25939",
    "SETDB2": "83852",
    "UHRF1": "29128",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--matrix",
        default="data/processed/ortholog_matrix_primary_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv",
    )
    parser.add_argument("--datasets", default="tools/ncbi_datasets/datasets.exe")
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument("--cache-dir", default="data/raw/protein_conservation")
    parser.add_argument("--work-dir", default="data/interim/protein_conservation")
    parser.add_argument(
        "--rows-output", default="results/tables/targeted_protein_conservation_rows.tsv"
    )
    parser.add_argument(
        "--gene-output", default="results/tables/targeted_protein_conservation_by_gene.tsv"
    )
    parser.add_argument(
        "--report-output", default="results/reports/targeted_protein_conservation_report.md"
    )
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--force-download", action="store_true")
    return parser.parse_args()


def clean_gene_id(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if not text:
        return None
    match = re.match(r"^(\d+)(?:\.0+)?$", text)
    return match.group(1) if match else None


def read_fasta_text(text: str):
    header = None
    sequence = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(sequence).replace("*", "").upper()
            header = line[1:]
            sequence = []
        else:
            sequence.append(line)
    if header is not None:
        yield header, "".join(sequence).replace("*", "").upper()


def read_fasta(path: Path):
    yield from read_fasta_text(path.read_text(encoding="utf-8"))


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def download_batches(
    gene_ids: list[str], datasets: Path, cache_dir: Path, batch_size: int, force: bool
) -> list[Path]:
    if not datasets.exists():
        raise FileNotFoundError(f"NCBI Datasets executable not found: {datasets}")
    cache_dir.mkdir(parents=True, exist_ok=True)
    batch_dirs = []
    for batch_index, start in enumerate(range(0, len(gene_ids), batch_size), 1):
        ids = gene_ids[start : start + batch_size]
        batch_dir = cache_dir / f"batch_{batch_index:02d}"
        marker = batch_dir / "ncbi_dataset" / "data" / "protein.faa"
        if force and batch_dir.exists():
            shutil.rmtree(batch_dir)
        if marker.exists():
            batch_dirs.append(batch_dir)
            continue
        batch_dir.mkdir(parents=True, exist_ok=True)
        id_file = batch_dir / "gene_ids.txt"
        archive = batch_dir / "ncbi_dataset.zip"
        id_file.write_text("\n".join(ids) + "\n", encoding="ascii")
        command = [
            str(datasets),
            "download",
            "gene",
            "gene-id",
            "--inputfile",
            str(id_file),
            "--include",
            "protein,product-report",
            "--filename",
            str(archive),
            "--no-progressbar",
        ]
        subprocess.run(command, check=True)
        with zipfile.ZipFile(archive) as zipped:
            zipped.extractall(batch_dir)
        if not marker.exists():
            raise RuntimeError(f"Downloaded batch lacks protein.faa: {batch_dir}")
        batch_dirs.append(batch_dir)
    return batch_dirs


def collect_longest_ncbi_proteins(batch_dirs: list[Path]) -> dict[str, dict[str, object]]:
    best: dict[str, dict[str, object]] = {}
    for batch_dir in batch_dirs:
        fasta = batch_dir / "ncbi_dataset" / "data" / "protein.faa"
        for header, sequence in read_fasta(fasta):
            match = re.search(r"\[GeneID=(\d+)\]", header)
            if not match or not sequence:
                continue
            gene_id = match.group(1)
            accession = header.split()[0]
            candidate = {
                "sequence": sequence,
                "accession": accession,
                "header": header,
                "length": len(sequence),
            }
            candidate_rank = (
                int(gene_id in HUMAN_GENE_IDS.values() and accession.startswith("NP_")),
                candidate["length"],
            )
            current_rank = (-1, -1)
            if gene_id in best:
                current_accession = str(best[gene_id]["accession"])
                current_rank = (
                    int(gene_id in HUMAN_GENE_IDS.values() and current_accession.startswith("NP_")),
                    int(best[gene_id]["length"]),
                )
            if candidate_rank > current_rank:
                best[gene_id] = candidate
    return best


def truthy(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.lower().isin({"true", "1", "yes"})


def load_strict_rescue_sequences() -> dict[tuple[str, str], dict[str, object]]:
    specs = [
        (
            Path("data/processed/phase3_gff_rescue_protein_sequences.tsv"),
            Path("data/processed/phase3_gff_rescue_sequence_decisions.tsv"),
            "can_count_as_strict_sequence_after_gff_sequence",
            "protein_id",
            "strict_gff_sequence",
        ),
        (
            Path("data/processed/phase3_cds_translation_rescue_sequences.tsv"),
            Path("data/processed/phase3_cds_translation_rescue_sequence_decisions.tsv"),
            "can_count_as_strict_sequence_after_cds_translation",
            "fetched_accession",
            "strict_cds_translation",
        ),
        (
            Path("data/processed/phase3_full_length_uniprot_rescue_sequences.tsv"),
            Path("data/processed/phase3_full_length_uniprot_rescue_sequence_decisions.tsv"),
            "can_count_as_strict_sequence_after_uniprot_rescue",
            "fetched_accession",
            "strict_uniprot_sequence",
        ),
    ]
    rescue: dict[tuple[str, str], dict[str, object]] = {}
    for sequence_path, decision_path, strict_column, accession_column, source in specs:
        if not sequence_path.exists() or not decision_path.exists():
            continue
        sequences = pd.read_csv(sequence_path, sep="\t", dtype=str)
        decisions = pd.read_csv(decision_path, sep="\t", dtype=str)
        keys = ["scientific_name", "human_gene_symbol"]
        strict = decisions.loc[truthy(decisions[strict_column]), keys].drop_duplicates()
        candidates = sequences.merge(strict, on=keys, how="inner")
        for _, row in candidates.iterrows():
            sequence = str(row.get("sequence", "")).replace("*", "").strip().upper()
            if not sequence or sequence == "NAN":
                continue
            key = (row["scientific_name"], row["human_gene_symbol"])
            candidate = {
                "sequence": sequence,
                "accession": str(row.get(accession_column, "")),
                "header": str(row.get("fasta_header", "")),
                "length": len(sequence),
                "source": source,
            }
            if key not in rescue or candidate["length"] > rescue[key]["length"]:
                rescue[key] = candidate
    return rescue


def align_gene(mafft: Path, input_path: Path, output_path: Path) -> None:
    command = f'"{mafft}" --auto --quiet "{input_path}"'
    with output_path.open("w", encoding="ascii", newline="\n") as handle:
        completed = subprocess.run(
            command,
            shell=True,
            stdout=handle,
            stderr=subprocess.PIPE,
            text=True,
        )
    if completed.returncode != 0:
        raise RuntimeError(f"MAFFT failed for {input_path}: {completed.stderr}")


def alignment_metrics(reference: str, target: str) -> dict[str, float]:
    if len(reference) != len(target):
        raise ValueError("Aligned reference and target have different lengths")
    reference_sites = sum(aa != "-" for aa in reference)
    paired = [(r, t) for r, t in zip(reference, target) if r != "-" and t != "-"]
    aligned_sites = len(paired)
    identical = sum(r == t for r, t in paired)
    coverage = aligned_sites / reference_sites if reference_sites else math.nan
    identity = identical / aligned_sites if aligned_sites else math.nan
    return {
        "human_reference_coverage": coverage,
        "aligned_identity": identity,
        "identity_coverage_product": identity * coverage,
        "aligned_sites": aligned_sites,
        "identical_sites": identical,
    }


def main() -> None:
    args = parse_args()
    matrix = pd.read_csv(args.matrix, sep="\t", dtype=str)
    focal = matrix.loc[matrix["human_gene_symbol"].isin(FOCAL_GENES)].copy()
    focal["ortholog_gene_id_clean"] = focal["ortholog_gene_id"].map(clean_gene_id)

    all_gene_ids = sorted(
        set(focal["ortholog_gene_id_clean"].dropna()) | set(HUMAN_GENE_IDS.values()),
        key=int,
    )
    batch_dirs = download_batches(
        all_gene_ids,
        Path(args.datasets),
        Path(args.cache_dir),
        args.batch_size,
        args.force_download,
    )
    ncbi_sequences = collect_longest_ncbi_proteins(batch_dirs)
    rescue_sequences = load_strict_rescue_sequences()

    missing_human = [gene for gene, gene_id in HUMAN_GENE_IDS.items() if gene_id not in ncbi_sequences]
    if missing_human:
        raise RuntimeError(f"Missing human reference proteins: {', '.join(missing_human)}")

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    result_rows = []

    for gene in FOCAL_GENES:
        gene_rows = focal.loc[focal["human_gene_symbol"].eq(gene)].copy()
        human_reference = ncbi_sequences[HUMAN_GENE_IDS[gene]]
        records = [("REF_Homo_sapiens", human_reference["sequence"])]
        sequence_metadata = {}
        for row_index, (_, row) in enumerate(gene_rows.iterrows(), 1):
            species = row["scientific_name"]
            gene_id = row["ortholog_gene_id_clean"]
            selected = ncbi_sequences.get(gene_id) if gene_id else None
            source = "ncbi_gene_longest_isoform" if selected else ""
            if selected is None:
                selected = rescue_sequences.get((species, gene))
                source = str(selected.get("source", "")) if selected else ""
            record_id = f"S{row_index:03d}"
            if selected:
                records.append((record_id, str(selected["sequence"])))
                sequence_metadata[record_id] = (row, selected, source)
            else:
                result_rows.append(
                    {
                        "scientific_name": species,
                        "opentree_tip_label": row.get("opentree_tip_label", ""),
                        "clade": row["clade"],
                        "flight_status": row["flight_status"],
                        "human_gene_symbol": gene,
                        "alignment_record_id": "",
                        "ortholog_gene_id": gene_id or "",
                        "ortholog_status": row.get("ortholog_status", ""),
                        "ortholog_confidence": row.get("week4_candidate_confidence", ""),
                        "sequence_available": False,
                        "sequence_source": "",
                        "selected_accession": "",
                        "protein_length": math.nan,
                        "human_reference_accession": human_reference["accession"],
                        "human_reference_length": human_reference["length"],
                        "protein_length_ratio": math.nan,
                        "length_completeness": math.nan,
                        "human_reference_coverage": math.nan,
                        "aligned_identity": math.nan,
                        "identity_coverage_product": math.nan,
                        "aligned_sites": math.nan,
                        "identical_sites": math.nan,
                        "sequence_missing_reason": "no_ncbi_gene_protein_or_strict_rescue",
                    }
                )

        input_path = work_dir / f"{gene}.raw.faa"
        aligned_path = work_dir / f"{gene}.aligned.faa"
        write_fasta(input_path, records)
        align_gene(Path(args.mafft), input_path, aligned_path)
        aligned = dict(read_fasta(aligned_path))
        reference_aligned = aligned["REF_Homo_sapiens"]
        for record_id, (row, selected, source) in sequence_metadata.items():
            metrics = alignment_metrics(reference_aligned, aligned[record_id])
            length_ratio = selected["length"] / human_reference["length"]
            length_completeness = min(length_ratio, 1 / length_ratio) if length_ratio > 0 else math.nan
            result_rows.append(
                {
                    "scientific_name": row["scientific_name"],
                    "opentree_tip_label": row.get("opentree_tip_label", ""),
                    "clade": row["clade"],
                    "flight_status": row["flight_status"],
                    "human_gene_symbol": gene,
                    "alignment_record_id": record_id,
                    "ortholog_gene_id": row["ortholog_gene_id_clean"] or "",
                    "ortholog_status": row.get("ortholog_status", ""),
                    "ortholog_confidence": row.get("week4_candidate_confidence", ""),
                    "sequence_available": True,
                    "sequence_source": source,
                    "selected_accession": selected["accession"],
                    "protein_length": selected["length"],
                    "human_reference_accession": human_reference["accession"],
                    "human_reference_length": human_reference["length"],
                    "protein_length_ratio": length_ratio,
                    "length_completeness": length_completeness,
                    **metrics,
                    "sequence_missing_reason": "",
                }
            )

    rows = pd.DataFrame(result_rows)
    rows = rows.sort_values(["human_gene_symbol", "scientific_name"]).reset_index(drop=True)
    rows_output = Path(args.rows_output)
    rows_output.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(rows_output, sep="\t", index=False, na_rep="")

    summaries = []
    for gene, group in rows.groupby("human_gene_symbol", sort=False):
        observed = group.loc[group["sequence_available"]]
        summaries.append(
            {
                "human_gene_symbol": gene,
                "species_total": len(group),
                "sequence_available_n": len(observed),
                "sequence_available_fraction": len(observed) / len(group),
                "ncbi_gene_sequence_n": (observed["sequence_source"] == "ncbi_gene_longest_isoform").sum(),
                "strict_rescue_sequence_n": (observed["sequence_source"] != "ncbi_gene_longest_isoform").sum(),
                "coverage_ge_0_5_n": (observed["human_reference_coverage"] >= 0.5).sum(),
                "coverage_ge_0_8_n": (observed["human_reference_coverage"] >= 0.8).sum(),
                "median_reference_coverage": observed["human_reference_coverage"].median(),
                "median_aligned_identity": observed["aligned_identity"].median(),
                "median_identity_coverage_product": observed["identity_coverage_product"].median(),
                "median_length_completeness": observed["length_completeness"].median(),
            }
        )
    by_gene = pd.DataFrame(summaries)
    gene_output = Path(args.gene_output)
    gene_output.parent.mkdir(parents=True, exist_ok=True)
    by_gene.to_csv(gene_output, sep="\t", index=False, na_rep="")

    source_counts = rows.loc[rows["sequence_available"], "sequence_source"].value_counts()
    report = [
        "# Targeted Protein-Conservation Evidence Layer",
        "",
        "Ten predeclared transposon/repeat/chromatin genes were evaluated across the 68-species primary panel. The longest NCBI Gene-linked protein isoform was used where available; only previously strict reciprocal sequence rescues were allowed as fallbacks. Missing sequence is retained as missing and is not interpreted as gene absence.",
        "",
        f"- Gene-species combinations: {len(rows)}.",
        f"- Protein sequences available: {int(rows['sequence_available'].sum())}/{len(rows)} ({rows['sequence_available'].mean():.1%}).",
        f"- Species with at least five focal proteins: {(rows.groupby('scientific_name')['sequence_available'].sum() >= 5).sum()}/68.",
        f"- Sequences with at least 50% human-reference coverage: {int((rows['human_reference_coverage'] >= 0.5).sum())}.",
        "",
        "## Sequence sources",
        "",
    ]
    report.extend(f"- {source}: {count}." for source, count in source_counts.items())
    report.extend(
        [
            "",
            "## Interpretation guardrails",
            "",
            "Protein identity and coverage are comparative sequence properties, not direct measurements of repair activity, transposon burden, or chromatin state. Gene-level tests are corrected for multiple comparisons, and the module-level aggregate requires at least five observed genes per species. Sequence availability is tested separately as an annotation-bias diagnostic.",
        ]
    )
    report_output = Path(args.report_output)
    report_output.parent.mkdir(parents=True, exist_ok=True)
    report_output.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Wrote {rows_output}, {gene_output}, and {report_output}")


if __name__ == "__main__":
    main()
