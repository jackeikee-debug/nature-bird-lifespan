"""Extract local protein FASTA candidates nominated by partial/family audit."""

from __future__ import annotations

import argparse
import gzip
import pathlib

import pandas as pd


TARGET_CLASSES = {"dnmt1_longer_local_isoform_available"}


def parse_fasta(path: pathlib.Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header = ""
    seq_parts: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if header:
                    records[header] = "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            elif line.strip():
                seq_parts.append(line.strip())
    if header:
        records[header] = "".join(seq_parts)
    return records


def accession_from_header(header: str) -> str:
    return str(header).split()[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=pathlib.Path, required=True)
    parser.add_argument("--metadata-output", type=pathlib.Path, required=True)
    parser.add_argument("--candidate-fasta-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    audit = pd.read_csv(args.audit, sep="\t", dtype=str).fillna("")
    targets = audit[audit["phase3_partial_family_resolution_class"].isin(TARGET_CLASSES)].copy()
    rows = []
    fasta_parts = []
    cache: dict[pathlib.Path, dict[str, str]] = {}
    for _, row in targets.iterrows():
        path = pathlib.Path(row["protein_asset_path"])
        header = row["protein_best_header"]
        status = "protein_asset_missing"
        sequence = ""
        if path.exists():
            if path not in cache:
                cache[path] = parse_fasta(path)
            sequence = cache[path].get(header, "")
            status = "local_alt_protein_found" if sequence else "local_alt_header_not_found"
        accession = accession_from_header(header)
        if sequence:
            safe_species = row["scientific_name"].replace(" ", "_")
            fasta_header = (
                f">{row['human_gene_symbol']}|{safe_species}|local_alt_protein:{accession}|"
                f"rank:1|{accession} {header}"
            )
            fasta_parts.append(fasta_header + "\n" + sequence)
        rows.append(
            {
                "scientific_name": row["scientific_name"],
                "human_gene_symbol": row["human_gene_symbol"],
                "best_assembly_accession": row["best_assembly_accession"],
                "source_resolution_class": row["phase3_partial_family_resolution_class"],
                "original_protein_id": row["protein_id"],
                "original_protein_length": row["protein_length"],
                "local_alt_protein_id": accession,
                "local_alt_header": header,
                "local_alt_protein_length": len(sequence),
                "sequence_fetch_status": status,
                "fetched_accession": accession,
            }
        )

    out = pd.DataFrame(rows)
    args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
    args.candidate_fasta_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.metadata_output, sep="\t", index=False)
    args.candidate_fasta_output.write_text("\n".join(fasta_parts) + ("\n" if fasta_parts else ""), encoding="utf-8")

    lines = [
        "# Phase 3 Local Alternative Protein Candidate Report",
        "",
        f"Audit rows selected: {len(targets)}",
        f"Candidate FASTA records: {len(fasta_parts)}",
        "",
        "## Candidates",
    ]
    for _, row in out.iterrows():
        lines.append(
            f"- {row['scientific_name']} / {row['human_gene_symbol']}: "
            f"{row['original_protein_id']} ({row['original_protein_length']} aa) -> "
            f"{row['local_alt_protein_id']} ({row['local_alt_protein_length']} aa), {row['sequence_fetch_status']}"
        )
    if out.empty:
        lines.append("- none: 0")
    lines.extend(
        [
            "",
            "## Interpretation",
            "These are local assembly protein records nominated by the partial/family audit as longer alternatives to the GFF-selected partial protein. They require reciprocal validation before any scoring change.",
            "",
            "## Outputs",
            f"- metadata: `{args.metadata_output}`",
            f"- candidate FASTA: `{args.candidate_fasta_output}`",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.metadata_output}")


if __name__ == "__main__":
    main()
