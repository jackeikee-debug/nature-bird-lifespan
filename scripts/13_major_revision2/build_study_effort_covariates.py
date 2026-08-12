"""Build reproducible phenotype study effort covariates for the 68-species panel.

The primary local proxies come directly from the AnAge record used to define
maximum lifespan: number of cited reference identifiers, sample-size category,
and data-quality category. An optional exact-binomial PubMed count is cached so
the analysis can distinguish database-record depth from broader publication
attention.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import pathlib
import re
import time
import urllib.parse
import urllib.request
import zipfile


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_anage(path: pathlib.Path) -> list[dict[str, str]]:
    with zipfile.ZipFile(path) as archive:
        raw = archive.read("anage_data.txt").decode("utf-8-sig")
    return list(csv.DictReader(io.StringIO(raw), delimiter="\t"))


def normalize_name(value: str) -> str:
    return re.sub(r"\s+", " ", value.strip()).casefold()


def count_references(value: str) -> int:
    return len({token.strip() for token in value.split(",") if token.strip()})


def safe_log1p(value: int | float | None) -> float | None:
    if value is None:
        return None
    return math.log1p(value)


def fetch_pubmed_count(name: str, cache_dir: pathlib.Path, delay: float) -> int:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / (re.sub(r"[^A-Za-z0-9_.-]+", "_", name) + ".json")
    if cache_path.exists():
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
        return int(payload["esearchresult"]["count"])
    term = f'"{name}"[Title/Abstract]'
    query = urllib.parse.urlencode(
        {"db": "pubmed", "term": term, "retmode": "json", "rettype": "count"}
    )
    request = urllib.request.Request(
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{query}",
        headers={"User-Agent": "vertebrate-lifespan-study-effort/1.0"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = json.loads(response.read().decode("utf-8"))
    cache_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    time.sleep(delay)
    return int(payload["esearchresult"]["count"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", type=pathlib.Path, default=pathlib.Path("data/processed/species_master.tsv"))
    parser.add_argument(
        "--panel",
        type=pathlib.Path,
        default=pathlib.Path("data/processed/maintenance_lifespan_phase2_W3_full_background_expanded.tsv"),
        help="Table defining the analyzed species panel; duplicate species rows are collapsed.",
    )
    parser.add_argument("--anage", type=pathlib.Path, default=pathlib.Path("data/raw/anage/anage_data.zip"))
    parser.add_argument("--output", type=pathlib.Path, default=pathlib.Path("data/processed/jme_revision2_study_effort_covariates.tsv"))
    parser.add_argument("--cache-dir", type=pathlib.Path, default=pathlib.Path("data/raw/pubmed_species_effort_cache"))
    parser.add_argument("--skip-pubmed", action="store_true")
    parser.add_argument("--delay", type=float, default=0.34)
    args = parser.parse_args()

    species = read_tsv(args.species)
    panel_names = {row["scientific_name"] for row in read_tsv(args.panel)}
    species = [row for row in species if row["scientific_name"] in panel_names]
    if len(species) != 68:
        raise RuntimeError(f"Expected 68 panel species, found {len(species)}")
    anage = read_anage(args.anage)
    by_name = {
        normalize_name(f"{row['Genus']} {row['Species']}"): row
        for row in anage
        if row.get("Genus") and row.get("Species")
    }
    sample_rank = {"unknown": 0, "small": 1, "medium": 2, "large": 3}
    quality_rank = {"low": 0, "questionable": 1, "acceptable": 2, "high": 3}

    output = []
    for index, row in enumerate(species, start=1):
        match_name = row.get("anage_matched_name") or row["scientific_name"]
        record = by_name.get(normalize_name(match_name))
        refs = count_references(record.get("References", "")) if record else 0
        sample = (record.get("Sample size", "") if record else "").strip().casefold()
        quality = (record.get("Data quality", "") if record else "").strip().casefold()
        pubmed_count = None
        if not args.skip_pubmed:
            pubmed_count = fetch_pubmed_count(row["scientific_name"], args.cache_dir, args.delay)
        output.append(
            {
                "scientific_name": row["scientific_name"],
                "anage_matched_name": match_name,
                "anage_reference_count": refs,
                "log1p_anage_reference_count": safe_log1p(refs),
                "anage_sample_size_category": sample,
                "anage_sample_size_rank": sample_rank.get(sample, ""),
                "anage_data_quality": quality,
                "anage_data_quality_rank": quality_rank.get(quality, ""),
                "pubmed_exact_binomial_count": "" if pubmed_count is None else pubmed_count,
                "log1p_pubmed_exact_binomial_count": "" if pubmed_count is None else safe_log1p(pubmed_count),
                "pubmed_query": f'"{row["scientific_name"]}"[Title/Abstract]',
                "anage_record_found": "yes" if record else "no",
            }
        )
        print(f"[{index}/{len(species)}] {row['scientific_name']}: AnAge refs={refs}, PubMed={pubmed_count}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
