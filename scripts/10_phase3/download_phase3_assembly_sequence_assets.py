"""Download Phase 3 assembly sequence assets listed in a manifest."""

from __future__ import annotations

import argparse
import pathlib
import time
import urllib.error
import urllib.request

import pandas as pd


def download_stream(url: str, path: pathlib.Path, timeout: int) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with urllib.request.urlopen(url, timeout=timeout) as response, tmp.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
    tmp.replace(path)
    return path.stat().st_size


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--log", type=pathlib.Path, required=True)
    parser.add_argument("--file-types", nargs="+", default=["protein", "cds", "rna"])
    parser.add_argument("--max-files", type=int, default=None)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--sleep", type=float, default=0.1)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str).fillna("")
    rows = manifest[manifest["file_type"].isin(args.file_types)].copy()
    if args.max_files is not None:
        rows = rows.head(args.max_files)

    log_rows = []
    for _, row in rows.iterrows():
        path = pathlib.Path(row["asset_local_path"])
        status = "cached"
        error = ""
        size = path.stat().st_size if path.exists() else 0
        if not path.exists():
            try:
                size = download_stream(row["asset_url"], path, args.timeout)
                status = "downloaded"
                if args.sleep:
                    time.sleep(args.sleep)
            except (urllib.error.URLError, TimeoutError, OSError, ValueError) as exc:
                status = "error"
                error = str(exc)
                size = 0
        log_rows.append(
            {
                **row.to_dict(),
                "download_status": status,
                "downloaded_bytes": str(size),
                "error": error,
            }
        )

    log = pd.DataFrame(log_rows)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    log.to_csv(args.log, sep="\t", index=False)

    counts = log["download_status"].value_counts().sort_index().to_dict() if not log.empty else {}
    print(f"Wrote {args.log}")
    print("; ".join(f"{key}={value}" for key, value in counts.items()))


if __name__ == "__main__":
    main()
