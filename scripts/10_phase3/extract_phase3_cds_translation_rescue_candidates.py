"""Extract assembly CDS translations for GFF hits that lacked protein_id."""

from __future__ import annotations

import argparse
import gzip
import pathlib
import re

import pandas as pd
from Bio.Seq import Seq


def parse_fasta(path: pathlib.Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    header = ""
    seq_parts: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
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


def bracket_attr(header: str, key: str) -> str:
    match = re.search(rf"\[{re.escape(key)}=([^\]]+)\]", header)
    return match.group(1) if match else ""


def mrna_to_locus_tag(value: str) -> str:
    match = re.search(r"([A-Z0-9]+_\d+)", str(value))
    return match.group(1) if match else ""


def find_cds_path(gff_path_text: str) -> pathlib.Path | None:
    gff_path = pathlib.Path(gff_path_text)
    directory = gff_path.parent
    candidates = sorted(directory.glob("*_cds_from_genomic.fna.gz"))
    return candidates[0] if candidates else None


def translate_cds(nt: str) -> str:
    seq = Seq(re.sub(r"[^ACGTNacgtn]", "", nt))
    usable = len(seq) - (len(seq) % 3)
    if usable <= 0:
        return ""
    aa = str(seq[:usable].translate(to_stop=False))
    if aa.endswith("*"):
        aa = aa[:-1]
    return aa.replace("*", "X")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff-sequence-metadata", type=pathlib.Path, required=True)
    parser.add_argument("--metadata-output", type=pathlib.Path, required=True)
    parser.add_argument("--candidate-fasta-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    metadata = pd.read_csv(args.gff_sequence_metadata, sep="\t", dtype=str).fillna("")
    targets = metadata[metadata["sequence_fetch_status"] == "no_protein_id_in_gff_cds"].copy()
    cache: dict[pathlib.Path, list[dict[str, str]]] = {}
    rows = []
    fasta_parts = []
    for _, row in targets.iterrows():
        cds_path = find_cds_path(row["gff_local_path"])
        locus_tag = mrna_to_locus_tag(row.get("matched_id", "")) or mrna_to_locus_tag(row.get("gff_mrna_id", ""))
        status = "missing_cds_asset"
        hit = None
        if cds_path and cds_path.exists():
            if cds_path not in cache:
                cache[cds_path] = parse_fasta(cds_path)
            for record in cache[cds_path]:
                if bracket_attr(record["header"], "locus_tag") == locus_tag:
                    hit = record
                    break
            status = "cds_locus_tag_not_found" if hit is None else "cds_translation_found"

        protein_id = bracket_attr(hit["header"], "protein_id") if hit else ""
        protein_product = bracket_attr(hit["header"], "protein") if hit else ""
        partial = bracket_attr(hit["header"], "partial") if hit else ""
        aa = translate_cds(hit["sequence"]) if hit else ""
        if hit and aa:
            safe_species = row["scientific_name"].replace(" ", "_")
            accession = protein_id or locus_tag
            header = (
                f">{row['human_gene_symbol']}|{safe_species}|assembly_cds_translation:{accession}|"
                f"rank:1|{accession} {hit['header']}"
            )
            fasta_parts.append(header + "\n" + aa)
        rows.append(
            {
                **row.to_dict(),
                "cds_local_path": str(cds_path or ""),
                "cds_locus_tag": locus_tag,
                "cds_header": hit["header"] if hit else "",
                "cds_protein_id": protein_id,
                "cds_product": protein_product,
                "cds_partial": partial,
                "sequence": aa,
                "protein_length": len(aa),
                "sequence_fetch_status": status if aa else ("cds_translation_empty" if hit else status),
                "fetched_accession": protein_id or locus_tag,
                "protein_id": protein_id or row.get("protein_id", ""),
            }
        )

    out = pd.DataFrame(rows)
    args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
    args.candidate_fasta_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.metadata_output, sep="\t", index=False)
    args.candidate_fasta_output.write_text("\n".join(fasta_parts) + ("\n" if fasta_parts else ""), encoding="utf-8")

    counts = out["sequence_fetch_status"].value_counts().sort_index() if not out.empty else {}
    lines = [
        "# Phase 3 CDS Translation Rescue Candidate Report",
        "",
        f"GFF metadata rows with no protein_id: {len(targets)}",
        f"CDS translation candidate rows: {len(out)}",
        f"Candidate FASTA records: {len(fasta_parts)}",
        "",
        "## Fetch Status",
    ]
    for status, count in getattr(counts, "items", lambda: [])():
        lines.append(f"- {status}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "This pass recovers candidate protein translations from assembly CDS FASTA by matching GFF mRNA-derived locus tags. It is intended for rows where the GFF hit was real but the matched feature did not expose a CDS protein_id to the earlier extractor.",
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
