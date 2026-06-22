"""Build a manifest for Phase 3 assembly protein/CDS/RNA sequence assets."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


FILE_SUFFIXES = {
    "protein": "_protein.faa.gz",
    "cds": "_cds_from_genomic.fna.gz",
    "rna": "_rna.fna.gz",
}


def asset_url(base_url: str, file_type: str) -> str:
    base_url = str(base_url).strip()
    if not base_url:
        return ""
    prefix = base_url.rsplit("/", 1)[0]
    basename = base_url.rsplit("/", 1)[1].removesuffix("_genomic.gff.gz")
    return f"{prefix}/{basename}{FILE_SUFFIXES[file_type]}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff-input", type=pathlib.Path, required=True)
    parser.add_argument("--annotation-manifest", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, default=pathlib.Path("data/raw/annotation_rescue"))
    args = parser.parse_args()

    gff_input = pd.read_csv(args.gff_input, sep="\t", dtype=str).fillna("")
    annotation = pd.read_csv(args.annotation_manifest, sep="\t", dtype=str).fillna("")
    target_species = set(gff_input["scientific_name"])
    target_assemblies = set(gff_input["best_assembly_accession"])
    rows = annotation[
        annotation["scientific_name"].isin(target_species)
        & annotation["best_assembly_accession"].isin(target_assemblies)
    ].copy()

    manifest_rows = []
    for _, row in rows.sort_values(["scientific_name"]).iterrows():
        for file_type in ["protein", "cds", "rna"]:
            url = asset_url(row["gff_url"], file_type)
            local_name = url.rsplit("/", 1)[1] if url else ""
            local_path = args.output_dir / row["best_assembly_accession"] / local_name if local_name else pathlib.Path("")
            manifest_rows.append(
                {
                    "scientific_name": row["scientific_name"],
                    "clade": row.get("clade", ""),
                    "flight_status": row.get("flight_status", ""),
                    "species_taxid": row.get("species_taxid", ""),
                    "best_assembly_accession": row["best_assembly_accession"],
                    "file_type": file_type,
                    "asset_url": url,
                    "asset_local_path": str(local_path),
                    "asset_status_pre_download": "cached" if local_path and local_path.exists() else "missing",
                }
            )
    manifest = pd.DataFrame(manifest_rows)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.output, sep="\t", index=False)

    counts = manifest.groupby(["file_type", "asset_status_pre_download"], as_index=False).size()
    lines = [
        "# Phase 3 Assembly Sequence Asset Manifest Report",
        "",
        f"Target species: {len(target_species)}",
        f"Target assemblies: {len(target_assemblies)}",
        f"Manifest rows: {len(manifest)}",
        "",
        "## Asset Status Before Download",
    ]
    for _, row in counts.iterrows():
        lines.append(f"- {row['file_type']} / {row['asset_status_pre_download']}: {row['size']}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "This manifest complements the GFF-only rescue pass by adding assembly-level protein, CDS, and RNA FASTA assets. These files are needed for full-length extraction and for rows where GFF annotations lack a CDS protein_id.",
            "",
            "## Outputs",
            f"- manifest: `{args.output}`",
        ]
    )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
