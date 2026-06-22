"""Suggest AnAge-backed replacements for species with incomplete life-history data."""

from __future__ import annotations

import argparse
import csv
import pathlib


DEFAULT_TAXONOMY_AUDIT = pathlib.Path("data/processed/taxonomy_audit.tsv")
DEFAULT_ANAGE = pathlib.Path("data/interim/anage_raw.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/replacement_candidates.tsv")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_float(value: str | None) -> float | None:
    if value is None or value.strip() == "":
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if parsed > 0 else None


def scientific_name(row: dict[str, str]) -> str:
    return f"{row.get('Genus', '').strip()} {row.get('Species', '').strip()}".strip()


def model_ready(row: dict[str, str]) -> bool:
    mass = parse_float(row.get("Body mass (g)")) or parse_float(row.get("Adult weight (g)"))
    lifespan = parse_float(row.get("Maximum longevity (yrs)"))
    return bool(scientific_name(row) and mass and lifespan)


def match_rank(target: dict[str, str], candidate: dict[str, str]) -> tuple[int, str]:
    target_name = target["scientific_name"]
    target_genus = target_name.split()[0] if target_name.split() else ""
    if candidate.get("Genus", "") == target_genus:
        return 1, "same_genus"
    if target.get("anage_family") and candidate.get("Family") == target.get("anage_family"):
        return 2, "same_family"
    if target.get("anage_order") and candidate.get("Order") == target.get("anage_order"):
        return 3, "same_order"
    return 99, "none"


def infer_target_taxonomy(target: dict[str, str], anage_rows: list[dict[str, str]]) -> dict[str, str]:
    if target.get("anage_family") or target.get("anage_order"):
        return target
    genus = target["scientific_name"].split()[0] if target["scientific_name"].split() else ""
    genus_matches = [row for row in anage_rows if row.get("Genus") == genus]
    if genus_matches:
        target = dict(target)
        target["anage_order"] = genus_matches[0].get("Order", "")
        target["anage_family"] = genus_matches[0].get("Family", "")
    return target


def build_candidates(
    taxonomy_rows: list[dict[str, str]],
    anage_rows: list[dict[str, str]],
    max_per_species: int,
) -> list[dict[str, str]]:
    ready = [row for row in anage_rows if model_ready(row)]
    targets = [
        row
        for row in taxonomy_rows
        if row["missing_max_lifespan"] == "yes"
        or row["missing_body_mass"] == "yes"
        or row["anage_record_present"] == "no"
    ]

    rows = []
    for target in targets:
        target = infer_target_taxonomy(target, anage_rows)
        scored = []
        for candidate in ready:
            rank, rank_label = match_rank(target, candidate)
            if rank == 99:
                continue
            quality = candidate.get("Data quality", "").lower()
            quality_rank = {"high": 0, "acceptable": 1, "low": 2, "questionable": 3}.get(
                quality, 4
            )
            sample = candidate.get("Sample size", "").lower()
            sample_rank = {"huge": 0, "large": 1, "medium": 2, "small": 3}.get(sample, 4)
            scored.append((rank, quality_rank, sample_rank, scientific_name(candidate), rank_label, candidate))

        for rank, _, _, _, rank_label, candidate in sorted(scored)[:max_per_species]:
            candidate_name = scientific_name(candidate)
            use_case = "replacement_candidate"
            if target["scientific_name"] == "Python bivittatus" and candidate_name == "Python molurus":
                use_case = "legacy_taxonomy_candidate_not_default_alias"
            rows.append(
                {
                    "target_species": target["scientific_name"],
                    "target_issue": target["taxonomy_issue"],
                    "candidate_species": candidate_name,
                    "candidate_common_name": candidate.get("Common name", ""),
                    "match_rank": rank_label,
                    "candidate_class": candidate.get("Class", ""),
                    "candidate_order": candidate.get("Order", ""),
                    "candidate_family": candidate.get("Family", ""),
                    "body_mass_g": candidate.get("Body mass (g)", "")
                    or candidate.get("Adult weight (g)", ""),
                    "max_lifespan_years": candidate.get("Maximum longevity (yrs)", ""),
                    "data_quality": candidate.get("Data quality", ""),
                    "sample_size": candidate.get("Sample size", ""),
                    "specimen_origin": candidate.get("Specimen origin", ""),
                    "source": candidate.get("Source", ""),
                    "recommended_use": use_case,
                    "notes": "Use only after deciding to replace target species in the panel; do not impute target species values.",
                }
            )
    return rows


def write_tsv(rows: list[dict[str, str]], output: pathlib.Path) -> None:
    fields = [
        "target_species",
        "target_issue",
        "candidate_species",
        "candidate_common_name",
        "match_rank",
        "candidate_class",
        "candidate_order",
        "candidate_family",
        "body_mass_g",
        "max_lifespan_years",
        "data_quality",
        "sample_size",
        "specimen_origin",
        "source",
        "recommended_use",
        "notes",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--taxonomy-audit", type=pathlib.Path, default=DEFAULT_TAXONOMY_AUDIT)
    parser.add_argument("--anage", type=pathlib.Path, default=DEFAULT_ANAGE)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-per-species", type=int, default=8)
    args = parser.parse_args()

    rows = build_candidates(
        read_tsv(args.taxonomy_audit),
        read_tsv(args.anage),
        args.max_per_species,
    )
    write_tsv(rows, args.output)
    print(f"Wrote {args.output} with {len(rows)} replacement candidates")


if __name__ == "__main__":
    main()

