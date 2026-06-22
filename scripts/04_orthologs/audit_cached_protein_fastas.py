"""Audit cached protein FASTA gzip files for readability."""

from __future__ import annotations

import argparse
import gzip
import pathlib

import pandas as pd


DEFAULT_MANIFEST = pathlib.Path("data/processed/week4_transposon_ncbi_crosscheck_protein_fasta_manifest.tsv")
DEFAULT_FASTA_DIR = pathlib.Path("data/raw/week4_protein_fastas")
DEFAULT_OUTPUT = pathlib.Path("data/interim/week4_transposon_ncbi_crosscheck_fasta_integrity.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week4_transposon_ncbi_crosscheck_fasta_integrity_report.md")


def local_fasta_path(fasta_dir: pathlib.Path, accession: str, protein_url: str) -> pathlib.Path:
    basename = protein_url.rstrip("/").split("/")[-1]
    return fasta_dir / accession / basename


def audit_file(path: pathlib.Path) -> tuple[str, int, int, str]:
    if not path.exists():
        return "missing", 0, 0, ""
    headers = 0
    seq_lines = 0
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    headers += 1
                elif line.strip():
                    seq_lines += 1
        if headers == 0:
            return "invalid_no_headers", headers, seq_lines, ""
        return "ok", headers, seq_lines, ""
    except Exception as exc:  # noqa: BLE001
        return "error", headers, seq_lines, f"{type(exc).__name__}: {exc}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--fasta-dir", type=pathlib.Path, default=DEFAULT_FASTA_DIR)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t")
    rows = []
    for _, row in manifest.iterrows():
        path = local_fasta_path(args.fasta_dir, row["best_assembly_accession"], row["protein_url"])
        status, headers, seq_lines, error = audit_file(path)
        rows.append(
            {
                "scientific_name": row["scientific_name"],
                "best_assembly_accession": row["best_assembly_accession"],
                "protein_url": row["protein_url"],
                "protein_local_path": str(path),
                "integrity_status": status,
                "fasta_headers": headers,
                "sequence_lines": seq_lines,
                "file_size_bytes": path.stat().st_size if path.exists() else 0,
                "error": error,
            }
        )
    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    counts = out["integrity_status"].value_counts().sort_index()
    lines = [
        "# Protein FASTA Integrity Report",
        "",
        f"Files audited: {len(out)}",
        "",
        "## Status Counts",
        "",
    ]
    for status, count in counts.items():
        lines.append(f"- {status}: {count}")
    bad = out[out["integrity_status"] != "ok"]
    lines.extend(["", "## Non-OK Files", ""])
    if bad.empty:
        lines.append("- None")
    else:
        for _, rec in bad.iterrows():
            lines.append(f"- {rec['best_assembly_accession']} / {rec['scientific_name']}: {rec['integrity_status']} {rec['error']}")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
