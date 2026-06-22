"""Download GFF files listed in the annotation rescue manifest."""

from __future__ import annotations

import argparse
import csv
import pathlib
import urllib.request


DEFAULT_MANIFEST = pathlib.Path("data/processed/annotation_rescue_manifest.tsv")
DEFAULT_OUTPUT_DIR = pathlib.Path("data/raw/annotation_rescue")
DEFAULT_LOG = pathlib.Path("data/interim/annotation_rescue_download_log.tsv")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    fields = list(rows[0].keys()) if rows else ["scientific_name", "gff_local_path", "download_status", "error"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def local_gff_path(output_dir: pathlib.Path, row: dict[str, str]) -> pathlib.Path:
    accession = row["best_assembly_accession"]
    basename = row["gff_url"].rstrip("/").split("/")[-1]
    return output_dir / accession / basename


def download(url: str, path: pathlib.Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=120) as response:
        path.write_bytes(response.read())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=pathlib.Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--log", type=pathlib.Path, default=DEFAULT_LOG)
    parser.add_argument("--max-files", type=int, default=None)
    args = parser.parse_args()

    rows = read_tsv(args.manifest)
    if args.max_files is not None:
        rows = rows[: args.max_files]
    log_rows = []
    for row in rows:
        local_path = local_gff_path(args.output_dir, row)
        status = "cached"
        error = ""
        if not local_path.exists():
            try:
                download(row["gff_url"], local_path)
                status = "downloaded"
            except Exception as exc:  # noqa: BLE001 - keep per-file errors.
                status = "error"
                error = str(exc)
        log_rows.append(
            {
                "scientific_name": row["scientific_name"],
                "best_assembly_accession": row["best_assembly_accession"],
                "gff_url": row["gff_url"],
                "gff_local_path": str(local_path),
                "download_status": status,
                "error": error,
            }
        )
    write_tsv(args.log, log_rows)
    print(f"Wrote {args.log}")


if __name__ == "__main__":
    main()
