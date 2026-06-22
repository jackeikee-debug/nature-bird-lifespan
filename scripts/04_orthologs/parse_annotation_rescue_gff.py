"""Search downloaded GFF annotations for maintenance gene symbols."""

from __future__ import annotations

import argparse
import csv
import gzip
import pathlib
import re
import urllib.parse
from collections import Counter


DEFAULT_MANIFEST = pathlib.Path("data/processed/annotation_rescue_manifest.tsv")
DEFAULT_DOWNLOAD_LOG = pathlib.Path("data/interim/annotation_rescue_download_log.tsv")
DEFAULT_GENES = pathlib.Path("data/processed/maintenance_gene_sets.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/annotation_rescue_gene_hits.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/annotation_rescue_gene_hits_report.md")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    fields = list(rows[0].keys()) if rows else [
        "scientific_name",
        "human_gene_symbol",
        "rescue_status",
        "matched_attribute",
        "feature_id",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_attributes(text: str) -> dict[str, str]:
    attrs = {}
    for part in text.split(";"):
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        attrs[key] = urllib.parse.unquote(value)
    return attrs


def gene_patterns(symbol: str) -> list[re.Pattern[str]]:
    escaped = re.escape(symbol)
    return [
        re.compile(rf"(^|[,;\s]){escaped}($|[,;\s])", re.IGNORECASE),
        re.compile(rf"\b{escaped}\b", re.IGNORECASE),
    ]


def feature_matches_symbol(attrs: dict[str, str], symbol: str) -> tuple[bool, str]:
    fields = ["gene", "Name", "gene_synonym", "standard_name", "product", "description"]
    haystack = []
    for field in fields:
        if attrs.get(field):
            haystack.append(f"{field}={attrs[field]}")
    text = ";".join(haystack)
    for pattern in gene_patterns(symbol):
        if pattern.search(text):
            return True, text
    return False, ""


def search_gff(gff_path: pathlib.Path, symbols: list[str]) -> dict[str, dict[str, str]]:
    hits: dict[str, dict[str, str]] = {}
    with gzip.open(gff_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            feature_type = parts[2]
            if feature_type not in {"gene", "mRNA", "transcript", "CDS"}:
                continue
            attrs = parse_attributes(parts[8])
            for symbol in symbols:
                if symbol in hits:
                    continue
                matched, attribute_text = feature_matches_symbol(attrs, symbol)
                if matched:
                    hits[symbol] = {
                        "matched_feature_type": feature_type,
                        "matched_attribute": attribute_text,
                        "feature_id": attrs.get("ID", ""),
                        "feature_name": attrs.get("Name", attrs.get("gene", "")),
                    }
    return hits


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--download-log", type=pathlib.Path, default=DEFAULT_DOWNLOAD_LOG)
    parser.add_argument("--genes", type=pathlib.Path, default=DEFAULT_GENES)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    manifest = {row["scientific_name"]: row for row in read_tsv(args.manifest)}
    downloads = read_tsv(args.download_log)
    genes = read_tsv(args.genes)
    gene_symbols = [row["human_gene_symbol"] for row in genes]
    gene_modules = {row["human_gene_symbol"]: row["maintenance_module"] for row in genes}
    rows = []
    for row in downloads:
        species = row["scientific_name"]
        gff_path = pathlib.Path(row["gff_local_path"])
        hits = {}
        if row["download_status"] in {"downloaded", "cached"} and gff_path.exists():
            hits = search_gff(gff_path, gene_symbols)
        for symbol in gene_symbols:
            hit = hits.get(symbol, {})
            rows.append(
                {
                    "scientific_name": species,
                    "clade": manifest.get(species, {}).get("clade", ""),
                    "best_assembly_accession": row["best_assembly_accession"],
                    "human_gene_symbol": symbol,
                    "maintenance_module": gene_modules[symbol],
                    "rescue_status": "gff_symbol_hit" if symbol in hits else "not_found_in_gff_attributes",
                    "matched_feature_type": hit.get("matched_feature_type", ""),
                    "matched_attribute": hit.get("matched_attribute", ""),
                    "feature_id": hit.get("feature_id", ""),
                    "feature_name": hit.get("feature_name", ""),
                    "gff_local_path": row["gff_local_path"],
                }
            )
    write_tsv(args.output, rows)
    found = [row for row in rows if row["rescue_status"] == "gff_symbol_hit"]
    by_species = Counter(row["scientific_name"] for row in found)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        "\n".join(
            [
                "# Annotation Rescue Gene Hits Report",
                "",
                f"Species parsed: {len(downloads)}",
                f"Species-gene rows: {len(rows)}",
                f"GFF symbol hits: {len(found)}",
                "",
                "## Hits by Species",
                *([f"- {species}: {count}" for species, count in sorted(by_species.items())] or ["- none: 0"]),
                "",
                "## Interpretation",
                "These are direct GFF attribute hits for maintenance gene symbols. They rescue annotation presence for species where NCBI Gene symbol/taxid search returned zero candidates, but final orthology still requires cross-database validation or protein-level reciprocal search.",
            ]
        ),
        encoding="utf-8",
    )
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
