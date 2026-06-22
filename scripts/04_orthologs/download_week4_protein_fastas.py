"""Download protein FASTA files for Week 4 validation batches."""

from __future__ import annotations

import argparse
import csv
import pathlib
import urllib.request


DEFAULT_MANIFEST = pathlib.Path("data/processed/week4_transposon_protein_fasta_manifest.tsv")
DEFAULT_OUTPUT_DIR = pathlib.Path("data/raw/week4_protein_fastas")
DEFAULT_LOG = pathlib.Path("data/interim/week4_transposon_protein_fasta_download_log.tsv")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    fields = list(rows[0].keys()) if rows else ["download_status", "error"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def local_fasta_path(output_dir: pathlib.Path, row: dict[str, str]) -> pathlib.Path:
    accession = row["best_assembly_accession"]
    basename = row["protein_url"].rstrip("/").split("/")[-1]
    return output_dir / accession / basename


def download_stream(url: str, path: pathlib.Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with urllib.request.urlopen(url, timeout=180) as response, tmp_path.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
    tmp_path.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=pathlib.Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--log", type=pathlib.Path, default=DEFAULT_LOG)
    parser.add_argument("--max-files", type=int, default=None)
    args = parser.parse_args()

    rows = read_tsv(args.manifest)
    rows = [row for row in rows if row.get("download_priority") == "download"]
    if args.max_files is not None:
        rows = rows[: args.max_files]

    log_rows = []
    for row in rows:
        local_path = local_fasta_path(args.output_dir, row)
        status = "cached"
        error = ""
        if not local_path.exists() or local_path.stat().st_size == 0:
            try:
                download_stream(row["protein_url"], local_path)
                status = "downloaded"
            except Exception as exc:  # noqa: BLE001
                status = "error"
                error = str(exc)
        log_rows.append(
            {
                "scientific_name": row["scientific_name"],
                "clade": row["clade"],
                "species_taxid": row["species_taxid"],
                "best_assembly_accession": row["best_assembly_accession"],
                "protein_url": row["protein_url"],
                "protein_local_path": str(local_path),
                "download_status": status,
                "file_size_bytes": str(local_path.stat().st_size if local_path.exists() else 0),
                "error": error,
            }
        )
        write_tsv(args.log, log_rows)
    write_tsv(args.log, log_rows)
    print(f"Wrote {args.log}")


if __name__ == "__main__":
    main()
