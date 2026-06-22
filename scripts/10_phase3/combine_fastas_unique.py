"""Combine FASTA files while keeping unique headers."""

from __future__ import annotations

import argparse
import pathlib


def read_records(path: pathlib.Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    seq_parts: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(seq_parts)))
            header = line[1:].strip()
            seq_parts = []
        elif line.strip():
            seq_parts.append(line.strip())
    if header:
        records.append((header, "".join(seq_parts)))
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    seen = set()
    records = []
    input_counts = []
    for path in args.inputs:
        path_records = read_records(path) if path.exists() else []
        input_counts.append((path, len(path_records)))
        for header, seq in path_records:
            if header in seen:
                continue
            seen.add(header)
            records.append((header, seq))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "\n".join(f">{header}\n{seq}" for header, seq in records) + ("\n" if records else ""),
        encoding="utf-8",
    )
    by_gene: dict[str, int] = {}
    for header, _ in records:
        gene = header.split("|", 1)[0]
        by_gene[gene] = by_gene.get(gene, 0) + 1
    lines = [
        "# Combined FASTA Report",
        "",
        f"Input files: {len(args.inputs)}",
        f"Unique FASTA records: {len(records)}",
        "",
        "## Input Counts",
    ]
    for path, count in input_counts:
        lines.append(f"- {path}: {count}")
    lines.extend(["", "## Records by Gene"])
    for gene, count in sorted(by_gene.items()):
        lines.append(f"- {gene}: {count}")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
