"""Download the latest stable AnAge dataset zip.

The script intentionally stores the original zip unchanged under data/raw so
that downstream parsing remains reproducible.
"""

from __future__ import annotations

import argparse
import pathlib
import urllib.request


DEFAULT_URL = "https://genomics.senescence.info/species/dataset.zip"
DEFAULT_OUTPUT = pathlib.Path("data/raw/anage/anage_data.zip")


def download(url: str, output: pathlib.Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "bird-lifespan-genome-maintenance/0.1"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        output.write_bytes(response.read())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", default=DEFAULT_URL)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    download(args.url, args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()

