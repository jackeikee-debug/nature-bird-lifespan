"""Extract the main tab-delimited AnAge file from the downloaded zip."""

from __future__ import annotations

import argparse
import pathlib
import zipfile


DEFAULT_INPUT = pathlib.Path("data/raw/anage/anage_data.zip")
DEFAULT_OUTPUT = pathlib.Path("data/interim/anage_raw.tsv")


def extract(input_zip: pathlib.Path, output: pathlib.Path) -> None:
    if not input_zip.exists():
        raise FileNotFoundError(
            f"Missing {input_zip}. Run scripts/00_download/download_anage.py first."
        )
    with zipfile.ZipFile(input_zip) as archive:
        candidates = [
            name
            for name in archive.namelist()
            if name.lower().endswith((".txt", ".tsv")) and not name.endswith("/")
        ]
        if not candidates:
            raise RuntimeError(f"No tab-delimited data file found in {input_zip}")
        data_name = max(candidates, key=lambda name: archive.getinfo(name).file_size)
        text = archive.read(data_name).decode("utf-8-sig", errors="replace")

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    extract(args.input, args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()

