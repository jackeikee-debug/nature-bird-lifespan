#!/usr/bin/env python
"""Prepare protein FASTA inputs for targeted gene-family phylogenies.

The goal is not to rebuild genome-wide orthology. These small trees are a
sequence-validation supplement for paralog-prone families in the
transposon/repeat/chromatin module.
"""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(".")
OUT_DIR = ROOT / "data" / "interim" / "phase3" / "gene_family_trees"
MANIFEST = OUT_DIR / "gene_family_tree_sequence_manifest.tsv"

INPUT_FASTAS = [
    ROOT / "data" / "interim" / "phase3" / "phase3_priority1_full_gff_rescue_candidates.faa",
    ROOT / "data" / "interim" / "phase3" / "phase3_priority1_full_gff_rescue_human_refs_combined.faa",
    ROOT / "data" / "interim" / "phase3" / "phase3_full_length_uniprot_rescue_candidates.faa",
    ROOT / "data" / "interim" / "phase3" / "phase3_full_length_uniprot_rescue_human_refs.faa",
    ROOT / "data" / "interim" / "phase3" / "phase3_cds_translation_rescue_candidates.faa",
    ROOT / "data" / "interim" / "phase3" / "phase3_local_alt_protein_candidates.faa",
    ROOT / "data" / "interim" / "phase3" / "manual_refs" / "human_DNMT3B_uniprot_Q9UBC3.faa",
]

FAMILIES = {
    "DNMT_family": {"DNMT1", "DNMT3A", "DNMT3B"},
    "MBD_family": {"MBD2", "MBD3"},
    "HELLS_support_tree": {"HELLS"},
    "UHRF1_support_tree": {"UHRF1"},
    "SETDB2_support_tree": {"SETDB2"},
    "SAMHD1_support_tree": {"SAMHD1"},
}


def read_fasta(path: Path):
    if not path.exists():
        return
    header = None
    seq_parts: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_parts)
            header = line[1:]
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def parse_header(header: str, source_file: Path):
    gene = ""
    species = ""
    accession = ""
    rank = ""

    if header.startswith(("sp|", "tr|")):
      # UniProt format: sp|Q9UBC3|DNM3B_HUMAN ... OS=Homo sapiens ... GN=DNMT3B
        parts = header.split("|")
        accession = parts[1] if len(parts) > 1 else ""
        match = re.search(r"\bGN=([A-Za-z0-9_.-]+)", header)
        gene = match.group(1).upper() if match else parts[2].split("_")[0].upper()
        os_match = re.search(r"\bOS=([^=]+?)\sOX=", header)
        species = os_match.group(1).strip().replace(" ", "_") if os_match else "Homo_sapiens"
        rank = "manual_uniprot"
    else:
        parts = header.split("|")
        gene = parts[0].strip().upper() if parts else ""
        species = parts[1].strip() if len(parts) > 1 else "unknown_species"
        for part in parts[2:]:
            if part.startswith("rank:"):
                rank = part.replace("rank:", "", 1)
            if ":" in part and not accession:
                accession = part.split(":", 1)[1].split()[0]
        if not accession:
            rest = parts[-1] if parts else header
            match = re.search(r"\b([A-Z]{2,4}_[0-9.]+|[NXWYZA][A-Z]{1,3}[0-9]{5,}(?:\.[0-9]+)?|[A-Z0-9]{6,10})\b", rest)
            accession = match.group(1) if match else "no_accession"

    if not rank:
        rank = "1"

    source_name = source_file.stem
    if "gff" in source_name:
        source_class = "gff"
    elif "uniprot" in source_name or header.startswith(("sp|", "tr|")):
        source_class = "uniprot"
    elif "cds" in source_name:
        source_class = "cds_translation"
    elif "local_alt" in source_name:
        source_class = "local_alt"
    else:
        source_class = "reference"

    return {
        "gene": gene,
        "species": species.replace(" ", "_"),
        "accession": accession,
        "rank": rank,
        "source_class": source_class,
        "source_file": str(source_file).replace("\\", "/"),
    }


def clean_id(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text[:120]


def wrap(seq: str, width: int = 70) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    records_by_family: dict[str, list[dict[str, str]]] = {family: [] for family in FAMILIES}
    seen = set()

    for path in INPUT_FASTAS:
        for header, seq in read_fasta(path) or []:
            meta = parse_header(header, path)
            gene = meta["gene"].upper()
            seq = re.sub(r"[^A-Za-z*]", "", seq).upper().replace("*", "")
            if not seq:
                continue
            # Human reference files often contain several isoforms. Keep rank 1
            # unless this is a manually added reference accession.
            if meta["species"] == "Homo_sapiens" and meta["rank"] not in {"1", "manual_uniprot"}:
                continue
            for family, genes in FAMILIES.items():
                if gene not in genes:
                    continue
                dedup_key = (family, gene, meta["species"], seq)
                if dedup_key in seen:
                    continue
                seen.add(dedup_key)
                base_id = clean_id(f"{gene}__{meta['species']}__{meta['source_class']}__{meta['accession']}")
                record = {
                    "family": family,
                    "output_id": base_id,
                    "gene": gene,
                    "species": meta["species"],
                    "source_class": meta["source_class"],
                    "source_file": meta["source_file"],
                    "accession": meta["accession"],
                    "rank": meta["rank"],
                    "length_aa": str(len(seq)),
                    "original_header": header,
                    "sequence": seq,
                }
                records_by_family[family].append(record)

    manifest_rows = []
    for family, records in records_by_family.items():
        id_counts: dict[str, int] = {}
        records = sorted(records, key=lambda r: (r["gene"], r["species"], r["source_class"], -int(r["length_aa"])))
        out_path = OUT_DIR / f"{family}_raw.faa"
        with out_path.open("w", encoding="ascii") as handle:
            for record in records:
                base = record["output_id"]
                id_counts[base] = id_counts.get(base, 0) + 1
                if id_counts[base] > 1:
                    record["output_id"] = f"{base}_{id_counts[base]}"
                handle.write(f">{record['output_id']}\n{wrap(record['sequence'])}\n")
                manifest_rows.append({k: v for k, v in record.items() if k != "sequence"})

    headers = [
        "family",
        "output_id",
        "gene",
        "species",
        "source_class",
        "source_file",
        "accession",
        "rank",
        "length_aa",
        "original_header",
    ]
    with MANIFEST.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(headers) + "\n")
        for row in manifest_rows:
            handle.write("\t".join(str(row.get(h, "")) for h in headers) + "\n")

    for family, records in records_by_family.items():
        genes = sorted({r["gene"] for r in records})
        species = sorted({r["species"] for r in records})
        print(f"{family}: {len(records)} sequences, {len(species)} species, genes={','.join(genes)}")
    print(f"Wrote {MANIFEST}")


if __name__ == "__main__":
    main()
