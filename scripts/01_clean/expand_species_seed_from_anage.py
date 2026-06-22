"""Create an expanded species seed table from AnAge.

The expanded seed keeps all manually curated species from config/species_seed.tsv
and then adds model-ready AnAge records with broad clade and flight labels.
Automatically added rows deliberately keep habitat, diet, and genome annotation
fields as unknown until they are curated from dedicated sources.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import defaultdict, deque


DEFAULT_ANAGE = pathlib.Path("data/interim/anage_raw.tsv")
DEFAULT_MANUAL_SEED = pathlib.Path("config/species_seed.tsv")
DEFAULT_OUTPUT = pathlib.Path("config/species_seed_expanded.tsv")

FIELDNAMES = [
    "scientific_name",
    "clade",
    "flight_status",
    "habitat",
    "diet",
    "genome_available",
    "annotation_quality",
    "inclusion_role",
    "notes",
]

TARGETS = {
    "Aves": 100,
    "Mammalia_Chiroptera": 40,
    "Mammalia_nonChiroptera": 60,
    "Reptilia": 40,
}

FLIGHTLESS_BIRD_ORDERS = {
    "Apterygiformes",
    "Casuariiformes",
    "Rheiformes",
    "Sphenisciformes",
    "Struthioniformes",
}


def parse_float(value: str) -> float | None:
    text = value.strip().replace(",", "")
    if not text:
        return None
    try:
        parsed = float(text)
    except ValueError:
        return None
    return parsed if parsed > 0 else None


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def scientific_name(row: dict[str, str]) -> str:
    return f"{row['Genus'].strip()} {row['Species'].strip()}".strip()


def is_model_ready(row: dict[str, str]) -> bool:
    lifespan = parse_float(row.get("Maximum longevity (yrs)", ""))
    mass = parse_float(row.get("Body mass (g)", "")) or parse_float(
        row.get("Adult weight (g)", "")
    )
    name = scientific_name(row)
    return bool(name and lifespan and mass)


def project_group(row: dict[str, str]) -> str | None:
    if row["Class"] == "Aves":
        return "Aves"
    if row["Class"] == "Mammalia" and row["Order"] == "Chiroptera":
        return "Mammalia_Chiroptera"
    if row["Class"] == "Mammalia":
        return "Mammalia_nonChiroptera"
    if row["Class"] == "Reptilia":
        return "Reptilia"
    return None


def flight_status(row: dict[str, str], group: str) -> str:
    if group == "Aves":
        return "flightless" if row["Order"] in FLIGHTLESS_BIRD_ORDERS else "flighted"
    if group == "Mammalia_Chiroptera":
        return "powered_flight"
    return "non_flying"


def auto_row(row: dict[str, str], group: str) -> dict[str, str]:
    order = row["Order"].strip()
    family = row["Family"].strip()
    return {
        "scientific_name": scientific_name(row),
        "clade": group,
        "flight_status": flight_status(row, group),
        "habitat": "unknown",
        "diet": "unknown",
        "genome_available": "unknown",
        "annotation_quality": "unknown",
        "inclusion_role": "auto_anage_candidate",
        "notes": f"auto-selected from AnAge; order={order}; family={family}",
    }


def quality_rank(row: dict[str, str]) -> tuple[int, float, str]:
    quality = row.get("Data quality", "").strip().lower()
    quality_score = {"acceptable": 0, "questionable": 1}.get(quality, 2)
    sample_raw = row.get("Sample size", "").strip()
    sample = parse_float(sample_raw) or 0
    return quality_score, -sample, scientific_name(row)


def round_robin_select(rows: list[dict[str, str]], target: int) -> list[dict[str, str]]:
    buckets: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        buckets[(row["Order"], row["Family"])].append(row)

    queues = []
    for key in sorted(buckets):
        bucket = sorted(buckets[key], key=quality_rank)
        queues.append(deque(bucket))

    selected = []
    while queues and len(selected) < target:
        next_queues = []
        for queue in queues:
            if queue and len(selected) < target:
                selected.append(queue.popleft())
            if queue:
                next_queues.append(queue)
        queues = next_queues
    return selected


def build_expanded_seed(
    anage_path: pathlib.Path,
    manual_seed_path: pathlib.Path,
    targets: dict[str, int],
) -> list[dict[str, str]]:
    manual_rows = read_tsv(manual_seed_path)
    seen = {row["scientific_name"].lower() for row in manual_rows}
    anage_rows = [row for row in read_tsv(anage_path) if is_model_ready(row)]

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in anage_rows:
        group = project_group(row)
        if group and scientific_name(row).lower() not in seen:
            grouped[group].append(row)

    expanded = list(manual_rows)
    for group, target in targets.items():
        current = sum(row["clade"] == group for row in expanded)
        needed = max(0, target - current)
        for row in round_robin_select(grouped[group], needed):
            expanded.append(auto_row(row, group))
            seen.add(scientific_name(row).lower())
    return expanded


def write_tsv(rows: list[dict[str, str]], output: pathlib.Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--anage", type=pathlib.Path, default=DEFAULT_ANAGE)
    parser.add_argument("--manual-seed", type=pathlib.Path, default=DEFAULT_MANUAL_SEED)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--birds", type=int, default=TARGETS["Aves"])
    parser.add_argument("--bats", type=int, default=TARGETS["Mammalia_Chiroptera"])
    parser.add_argument("--mammals", type=int, default=TARGETS["Mammalia_nonChiroptera"])
    parser.add_argument("--reptiles", type=int, default=TARGETS["Reptilia"])
    args = parser.parse_args()

    targets = {
        "Aves": args.birds,
        "Mammalia_Chiroptera": args.bats,
        "Mammalia_nonChiroptera": args.mammals,
        "Reptilia": args.reptiles,
    }
    rows = build_expanded_seed(args.anage, args.manual_seed, targets)
    write_tsv(rows, args.output)
    counts = {group: sum(row["clade"] == group for row in rows) for group in targets}
    print(f"Wrote {args.output} with {len(rows)} species; counts={counts}")


if __name__ == "__main__":
    main()

