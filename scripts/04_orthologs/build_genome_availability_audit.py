"""Query NCBI Assembly metadata for PGLS species.

This script builds a lightweight genome availability audit without downloading
genomes. It uses NCBI E-utilities and caches JSON responses under data/raw.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re
import time
import urllib.parse
import urllib.request
import urllib.error
from typing import Any


DEFAULT_INPUT = pathlib.Path("data/processed/pgls_trait_table.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/genome_availability_audit.tsv")
DEFAULT_SUMMARY = pathlib.Path("results/reports/genome_availability_report.md")
DEFAULT_CACHE = pathlib.Path("data/raw/ncbi_assembly_cache")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("_")


def fetch_json(url: str, timeout: int) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "bird-lifespan-genome-audit/0.1"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def cached_fetch(path: pathlib.Path, url: str, timeout: int, delay: float) -> dict[str, Any]:
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    attempts = 0
    while True:
        attempts += 1
        try:
            data = fetch_json(url, timeout)
            break
        except urllib.error.HTTPError as exc:
            if exc.code != 429 or attempts >= 6:
                raise
            wait_seconds = max(delay, 2.0) * attempts
            print(f"NCBI rate limit hit; sleeping {wait_seconds:.1f}s before retry {attempts}/5")
            time.sleep(wait_seconds)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    time.sleep(delay)
    return data


def esearch_assembly(name: str, cache_dir: pathlib.Path, timeout: int, delay: float, retmax: int) -> list[str]:
    term = f'"{name}"[Organism] AND latest[filter]'
    params = urllib.parse.urlencode(
        {
            "db": "assembly",
            "term": term,
            "retmode": "json",
            "retmax": str(retmax),
        }
    )
    url = f"{EUTILS}/esearch.fcgi?{params}"
    data = cached_fetch(cache_dir / f"{slug(name)}.esearch.json", url, timeout, delay)
    return data.get("esearchresult", {}).get("idlist", [])


def esummary_assembly(
    name: str,
    ids: list[str],
    cache_dir: pathlib.Path,
    timeout: int,
    delay: float,
) -> list[dict[str, Any]]:
    if not ids:
        return []
    params = urllib.parse.urlencode(
        {
            "db": "assembly",
            "id": ",".join(ids),
            "retmode": "json",
        }
    )
    url = f"{EUTILS}/esummary.fcgi?{params}"
    data = cached_fetch(cache_dir / f"{slug(name)}.esummary.json", url, timeout, delay)
    result = data.get("result", {})
    rows = []
    for uid in result.get("uids", []):
        if uid in result:
            rows.append(result[uid])
    return rows


def as_int(value: Any) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def has_annotation(row: dict[str, Any]) -> bool:
    if row.get("annotrpturl"):
        return True
    busco = row.get("busco") or {}
    if busco.get("refseqannotationrelease"):
        return True
    properties = set(row.get("propertylist") or [])
    return any("annotation" in item.lower() for item in properties)


def assembly_score(row: dict[str, Any]) -> tuple:
    refseq_category = (row.get("refseq_category") or "").lower()
    accession = row.get("assemblyaccession") or ""
    status = row.get("assemblystatus") or ""
    category_score = {
        "reference genome": 0,
        "representative genome": 1,
        "na": 4,
        "": 5,
    }.get(refseq_category, 3)
    source_score = 0 if accession.startswith("GCF_") else 1
    status_score = {
        "Complete Genome": 0,
        "Chromosome": 1,
        "Scaffold": 2,
        "Contig": 3,
    }.get(status, 4)
    annotation_score = 0 if has_annotation(row) else 1
    scaffold_n50 = as_int(row.get("scaffoldn50"))
    contig_n50 = as_int(row.get("contign50"))
    release = row.get("asmupdatedate") or row.get("seqreleasedate") or ""
    return (
        category_score,
        source_score,
        status_score,
        annotation_score,
        -scaffold_n50,
        -contig_n50,
        release,
    )


def best_assembly(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    if not rows:
        return None
    return sorted(rows, key=assembly_score)[0]


def genome_tier(row: dict[str, Any] | None) -> str:
    if row is None:
        return "none"
    status = row.get("assemblystatus") or ""
    accession = row.get("assemblyaccession") or ""
    annotated = has_annotation(row)
    refseq = accession.startswith("GCF_")
    if refseq and annotated and status in {"Complete Genome", "Chromosome"}:
        return "tier1_refseq_annotated_chromosome"
    if annotated and status in {"Complete Genome", "Chromosome", "Scaffold"}:
        return "tier2_annotated"
    if status in {"Complete Genome", "Chromosome", "Scaffold"}:
        return "tier3_assembly_only"
    return "tier4_low_contiguity_or_unclear"


def flatten_species(row: dict[str, str], assemblies: list[dict[str, Any]]) -> dict[str, str]:
    best = best_assembly(assemblies)
    busco = (best or {}).get("busco") or {}
    synonym = (best or {}).get("synonym") or {}
    return {
        "scientific_name": row["scientific_name"],
        "opentree_tip_label": row["opentree_tip_label"],
        "clade": row["clade"],
        "flight_status": row["flight_status"],
        "ncbi_query_name": row.get("opentree_matched_name") or row["scientific_name"],
        "assembly_count_latest": str(len(assemblies)),
        "genome_available": "yes" if best else "no",
        "genome_analysis_tier": genome_tier(best),
        "best_assembly_accession": "" if not best else best.get("assemblyaccession", ""),
        "refseq_accession": "" if not best else synonym.get("refseq", ""),
        "genbank_accession": "" if not best else synonym.get("genbank", ""),
        "refseq_category": "" if not best else best.get("refseq_category", ""),
        "assembly_name": "" if not best else best.get("assemblyname", ""),
        "assembly_status": "" if not best else best.get("assemblystatus", ""),
        "assembly_type": "" if not best else best.get("assemblytype", ""),
        "taxid": "" if not best else str(best.get("taxid", "")),
        "species_taxid": "" if not best else str(best.get("speciestaxid", "")),
        "ncbi_species_name": "" if not best else best.get("speciesname", ""),
        "has_annotation_report": "yes" if has_annotation(best or {}) else "no",
        "annotation_report_url": "" if not best else best.get("annotrpturl", ""),
        "scaffold_n50": "" if not best else str(best.get("scaffoldn50", "")),
        "contig_n50": "" if not best else str(best.get("contign50", "")),
        "coverage": "" if not best else str(best.get("coverage", "")),
        "busco_lineage": str(busco.get("buscolineage", "")),
        "busco_complete": str(busco.get("complete", "")),
        "busco_singlecopy": str(busco.get("singlecopy", "")),
        "busco_duplicated": str(busco.get("duplicated", "")),
        "busco_fragmented": str(busco.get("fragmented", "")),
        "busco_missing": str(busco.get("missing", "")),
        "ftp_path_refseq": "" if not best else best.get("ftppath_refseq", ""),
        "ftp_path_genbank": "" if not best else best.get("ftppath_genbank", ""),
        "assembly_release_date": "" if not best else best.get("asmreleasedate_refseq", "") or best.get("asmreleasedate_genbank", ""),
        "assembly_update_date": "" if not best else best.get("asmupdatedate", ""),
        "manual_review_reason": review_reason(row, best),
    }


def review_reason(species_row: dict[str, str], assembly: dict[str, Any] | None) -> str:
    reasons = []
    if assembly is None:
        return "no_ncbi_assembly_found"
    if genome_tier(assembly).startswith("tier4"):
        reasons.append("low_or_unclear_assembly_tier")
    if not has_annotation(assembly):
        reasons.append("no_annotation_report")
    if assembly.get("assemblystatus") not in {"Complete Genome", "Chromosome", "Scaffold"}:
        reasons.append("low_assembly_level")
    if species_row.get("opentree_matched_name") and species_row["opentree_matched_name"] != species_row["scientific_name"]:
        reasons.append("opentree_synonym_name_used")
    if assembly.get("speciesname") and assembly.get("speciesname") != species_row.get("opentree_matched_name", species_row["scientific_name"]):
        reasons.append("ncbi_species_name_differs")
    return ";".join(reasons) if reasons else "none"


def write_report(rows: list[dict[str, str]], output: pathlib.Path) -> None:
    from collections import Counter

    tiers = Counter(row["genome_analysis_tier"] for row in rows)
    clade_tier = Counter((row["clade"], row["genome_analysis_tier"]) for row in rows)
    missing = [row for row in rows if row["genome_available"] == "no"]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        "\n".join(
            [
                "# Genome Availability Report",
                "",
                f"Species queried: {len(rows)}",
                f"Species with at least one latest NCBI assembly: {sum(row['genome_available'] == 'yes' for row in rows)}",
                f"Species without NCBI assembly: {len(missing)}",
                "",
                "## Genome Analysis Tier Counts",
                *[f"- {tier}: {count}" for tier, count in sorted(tiers.items())],
                "",
                "## Clade by Tier",
                *[f"- {clade} / {tier}: {count}" for (clade, tier), count in sorted(clade_tier.items())],
                "",
                "## Missing Assemblies",
                *[f"- {row['scientific_name']}" for row in missing[:50]],
                "",
                "## Notes",
                "- Tier 1: RefSeq, annotated, chromosome or complete genome.",
                "- Tier 2: annotated assembly, including scaffold-level assemblies.",
                "- Tier 3: assembly available but no clear annotation report.",
                "- Tier 4: low-contiguity or unclear assembly metadata.",
                "- This audit does not download genome FASTA or GFF files.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary", type=pathlib.Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--cache-dir", type=pathlib.Path, default=DEFAULT_CACHE)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--delay", type=float, default=0.75)
    parser.add_argument("--retmax", type=int, default=20)
    args = parser.parse_args()

    output_rows = []
    for index, species in enumerate(read_tsv(args.input), start=1):
        query_name = species.get("opentree_matched_name") or species["scientific_name"]
        ids = esearch_assembly(query_name, args.cache_dir, args.timeout, args.delay, args.retmax)
        assemblies = esummary_assembly(query_name, ids, args.cache_dir, args.timeout, args.delay)
        output_rows.append(flatten_species(species, assemblies))
        print(f"[{index}] {species['scientific_name']}: {len(assemblies)} assemblies")

    fields = list(output_rows[0].keys())
    write_tsv(args.output, output_rows, fields)
    write_report(output_rows, args.summary)
    print(f"Wrote {args.output} and {args.summary}")


if __name__ == "__main__":
    main()
