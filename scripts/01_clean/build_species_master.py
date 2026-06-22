"""Build the first analysis-ready species master table.

Inputs:
  - config/species_seed.tsv: project-specific labels and inclusion notes
  - data/raw/anage/anage_data.zip: zipped AnAge tab-delimited data

Output:
  - data/processed/species_master.tsv
"""

from __future__ import annotations

import argparse
import csv
import io
import pathlib
import re
import zipfile


DEFAULT_SEED = pathlib.Path("config/species_seed_expanded.tsv")
DEFAULT_ALIASES = pathlib.Path("config/taxonomy_aliases.tsv")
DEFAULT_SUPPLEMENTS = pathlib.Path("config/manual_life_history_supplements.tsv")
DEFAULT_ANAGE_ZIP = pathlib.Path("data/raw/anage/anage_data.zip")
DEFAULT_INTERIM = pathlib.Path("data/interim/anage_raw.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/species_master.tsv")

OUTPUT_FIELDS = [
    "scientific_name",
    "anage_matched_name",
    "common_name",
    "clade",
    "flight_status",
    "body_mass_g",
    "max_lifespan_years",
    "sexual_maturity_years",
    "habitat",
    "diet",
    "genome_available",
    "annotation_quality",
    "inclusion_role",
    "anage_record_present",
    "anage_body_mass_raw",
    "anage_lifespan_raw",
    "anage_maturity_raw",
    "life_history_data_source",
    "manual_supplement_source",
    "notes",
]


def normalize_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    text = value.strip()
    if not text or text.lower() in {"na", "n/a", "unknown", "not specified"}:
        return None
    match = re.search(r"-?\d+(?:\.\d+)?", text.replace(",", ""))
    if not match:
        return None
    return float(match.group(0))


def find_first(row: dict[str, str], candidates: list[str]) -> str:
    for key in candidates:
        if key in row and row[key].strip():
            return row[key].strip()
    return ""


def find_first_item(row: dict[str, str], candidates: list[str]) -> tuple[str, str]:
    for key in candidates:
        if key in row and row[key].strip():
            return key, row[key].strip()
    return "", ""


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_aliases(path: pathlib.Path) -> dict[str, str]:
    if not path.exists():
        return {}
    aliases = {}
    for row in read_tsv(path):
        aliases[row["project_name"].strip().lower()] = row["anage_name"].strip()
    return aliases


def read_supplements(path: pathlib.Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    supplements = {}
    for row in read_tsv(path):
        row = {key: (value or "") for key, value in row.items()}
        name = row.get("scientific_name", "").strip().lower()
        if name:
            supplements[name] = row
    return supplements


def read_anage_zip(path: pathlib.Path, interim_output: pathlib.Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Run scripts/00_download/download_anage.py first."
        )

    with zipfile.ZipFile(path) as archive:
        candidates = [
            name
            for name in archive.namelist()
            if name.lower().endswith((".txt", ".tsv")) and not name.endswith("/")
        ]
        if not candidates:
            raise RuntimeError(f"No tab-delimited data file found in {path}")
        data_name = max(candidates, key=lambda name: archive.getinfo(name).file_size)
        raw_bytes = archive.read(data_name)

    text = raw_bytes.decode("utf-8-sig", errors="replace")
    interim_output.parent.mkdir(parents=True, exist_ok=True)
    interim_output.write_text(text, encoding="utf-8")

    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    rows = []
    for row in reader:
        rows.append({normalize_header(key): value for key, value in row.items() if key})
    return rows


def index_anage(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    indexed = {}
    for row in rows:
        genus = find_first(row, ["genus"])
        species = find_first(row, ["species"])
        scientific_name = find_first(
            row,
            ["scientific_name", "binomial_name", "latin_name", "species_name"],
        )
        if not scientific_name and genus and species:
            scientific_name = f"{genus} {species}"
        if scientific_name:
            indexed[scientific_name.lower()] = row
    return indexed


def convert_mass_to_g(raw: str) -> float | None:
    value = parse_float(raw)
    if value is None:
        return None
    lowered = raw.lower()
    if "kg" in lowered:
        return value * 1000
    if "mg" in lowered:
        return value / 1000
    return value


def build_species_master(
    seed_path: pathlib.Path,
    aliases_path: pathlib.Path,
    supplements_path: pathlib.Path,
    anage_zip: pathlib.Path,
    interim_output: pathlib.Path,
) -> list[dict[str, str]]:
    seed_rows = read_tsv(seed_path)
    aliases = read_aliases(aliases_path)
    supplements = read_supplements(supplements_path)
    anage_rows = read_anage_zip(anage_zip, interim_output)
    anage_by_name = index_anage(anage_rows)

    output_rows = []
    for seed in seed_rows:
        scientific_name = seed["scientific_name"].strip()
        anage_matched_name = aliases.get(scientific_name.lower(), scientific_name)
        anage = anage_by_name.get(anage_matched_name.lower(), {})

        common_name = find_first(anage, ["common_name"])
        body_mass_raw = find_first(anage, ["body_mass_g", "body_mass", "adult_weight_g"])
        lifespan_raw = find_first(
            anage,
            ["maximum_longevity_yrs", "maximum_lifespan_yrs", "max_lifespan_yrs"],
        )
        maturity_key, maturity_raw = find_first_item(
            anage,
            [
                "female_maturity_days",
                "male_maturity_days",
                "sexual_maturity_days",
                "female_maturity_yrs",
                "male_maturity_yrs",
            ],
        )

        mass_g = convert_mass_to_g(body_mass_raw)
        lifespan_years = parse_float(lifespan_raw)
        maturity_value = parse_float(maturity_raw)
        maturity_years = ""
        if maturity_value is not None:
            maturity_years = maturity_value / 365.25 if "day" in maturity_key else maturity_value

        supplement = supplements.get(scientific_name.lower())
        supplement_source = ""
        source_parts = []
        if anage:
            source_parts.append("AnAge")
        if supplement and supplement.get("use_in_model", "").strip().lower() == "yes":
            supplement_source = supplement.get("source_url", "") or supplement.get("source_name", "")
            if mass_g is None:
                mass_g = parse_float(supplement.get("body_mass_g", ""))
                if mass_g is not None:
                    source_parts.append("manual_body_mass")
            if lifespan_years is None:
                lifespan_years = parse_float(supplement.get("max_lifespan_years", ""))
                if lifespan_years is not None:
                    source_parts.append("manual_lifespan")
            if maturity_years == "":
                supplemented_maturity = parse_float(
                    supplement.get("sexual_maturity_years", "")
                )
                if supplemented_maturity is not None:
                    maturity_years = supplemented_maturity
                    source_parts.append("manual_maturity")

        output_rows.append(
            {
                "scientific_name": scientific_name,
                "anage_matched_name": anage_matched_name if anage else "",
                "common_name": common_name,
                "clade": seed["clade"],
                "flight_status": seed["flight_status"],
                "body_mass_g": "" if mass_g is None else f"{mass_g:.6g}",
                "max_lifespan_years": "" if lifespan_years is None else f"{lifespan_years:.6g}",
                "sexual_maturity_years": ""
                if maturity_years == ""
                else f"{maturity_years:.6g}",
                "habitat": seed["habitat"],
                "diet": seed["diet"],
                "genome_available": seed["genome_available"],
                "annotation_quality": seed["annotation_quality"],
                "inclusion_role": seed["inclusion_role"],
                "anage_record_present": "yes" if anage else "no",
                "anage_body_mass_raw": body_mass_raw,
                "anage_lifespan_raw": lifespan_raw,
                "anage_maturity_raw": maturity_raw,
                "life_history_data_source": "+".join(source_parts) if source_parts else "",
                "manual_supplement_source": supplement_source,
                "notes": seed.get("notes", ""),
            }
        )
    return output_rows


def write_tsv(rows: list[dict[str, str]], output: pathlib.Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=pathlib.Path, default=DEFAULT_SEED)
    parser.add_argument("--aliases", type=pathlib.Path, default=DEFAULT_ALIASES)
    parser.add_argument("--supplements", type=pathlib.Path, default=DEFAULT_SUPPLEMENTS)
    parser.add_argument("--anage-zip", type=pathlib.Path, default=DEFAULT_ANAGE_ZIP)
    parser.add_argument("--interim-output", type=pathlib.Path, default=DEFAULT_INTERIM)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    rows = build_species_master(
        args.seed,
        args.aliases,
        args.supplements,
        args.anage_zip,
        args.interim_output,
    )
    write_tsv(rows, args.output)
    matched = sum(row["anage_record_present"] == "yes" for row in rows)
    print(f"Wrote {args.output} with {len(rows)} species; AnAge matches: {matched}")


if __name__ == "__main__":
    main()
