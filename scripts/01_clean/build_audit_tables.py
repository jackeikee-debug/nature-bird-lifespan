"""Build taxonomy and outlier audit tables for the Week 1 species panel."""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
import statistics
from collections import defaultdict


DEFAULT_SPECIES_MASTER = pathlib.Path("data/processed/species_master.tsv")
DEFAULT_RESIDUALS = pathlib.Path("data/processed/lifespan_residuals.tsv")
DEFAULT_ANAGE = pathlib.Path("data/interim/anage_raw.tsv")
DEFAULT_ALIASES = pathlib.Path("config/taxonomy_aliases.tsv")
DEFAULT_SUPPLEMENTS = pathlib.Path("config/manual_life_history_supplements.tsv")
DEFAULT_TAXONOMY_OUTPUT = pathlib.Path("data/processed/taxonomy_audit.tsv")
DEFAULT_OUTLIER_OUTPUT = pathlib.Path("data/processed/outlier_audit.tsv")


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
    if not math.isfinite(parsed):
        return None
    return parsed


def scientific_name(row: dict[str, str]) -> str:
    return f"{row.get('Genus', '').strip()} {row.get('Species', '').strip()}".strip()


def index_anage(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {scientific_name(row).lower(): row for row in rows if scientific_name(row)}


def read_aliases(path: pathlib.Path) -> dict[str, str]:
    if not path.exists():
        return {}
    return {
        row["project_name"].strip().lower(): row["anage_name"].strip()
        for row in read_tsv(path)
    }


def read_supplements(path: pathlib.Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    return {
        (row.get("scientific_name") or "").strip().lower(): {
            key: (value or "") for key, value in row.items()
        }
        for row in read_tsv(path)
        if (row.get("scientific_name") or "").strip()
    }


def candidate_names(name: str, anage_rows: list[dict[str, str]]) -> str:
    parts = name.split()
    if not parts:
        return ""
    genus = parts[0].lower()
    genus_matches = [
        scientific_name(row)
        for row in anage_rows
        if row.get("Genus", "").strip().lower() == genus
    ]
    if genus_matches:
        return "; ".join(genus_matches[:8])
    return ""


def taxonomy_issue(
    row: dict[str, str],
    anage: dict[str, str] | None,
    aliases: dict[str, str],
    supplement: dict[str, str] | None,
) -> tuple[str, str]:
    name = row["scientific_name"]
    has_mass = bool(row.get("body_mass_g", "").strip())
    has_lifespan = bool(row.get("max_lifespan_years", "").strip())
    if row.get("anage_record_present") == "no":
        if supplement and supplement.get("use_in_model", "").lower() == "yes":
            return "manual_supplement_available", "rebuild species_master with supplement enabled"
        return "missing_same_species_anage_record", "find primary source or replace with AnAge-covered species"
    if name.lower() in aliases:
        if not has_lifespan or not has_mass:
            return "alias_matched_but_incomplete", "audit missing AnAge fields"
        return "alias_applied", "keep alias and document taxonomy"
    if not has_lifespan and has_mass:
        if supplement and supplement.get("use_in_model", "").lower() == "yes":
            return "lifespan_supplement_available", "rebuild species_master with supplement enabled"
        return "missing_max_lifespan", "find primary longevity source or exclude from residual model"
    if has_lifespan and not has_mass:
        return "missing_body_mass", "find body mass source or exclude from residual model"
    if anage and anage.get("Data quality", "").strip().lower() == "questionable":
        return "questionable_anage_quality", "manual outlier and source audit recommended"
    return "ok", "none"


def write_taxonomy_audit(
    species_rows: list[dict[str, str]],
    anage_rows: list[dict[str, str]],
    aliases: dict[str, str],
    supplements: dict[str, dict[str, str]],
    output: pathlib.Path,
) -> None:
    anage_by_name = index_anage(anage_rows)
    fields = [
        "scientific_name",
        "anage_matched_name",
        "alias_applied",
        "anage_record_present",
        "hagrid",
        "anage_class",
        "anage_order",
        "anage_family",
        "anage_common_name",
        "clade",
        "flight_status",
        "body_mass_g",
        "max_lifespan_years",
        "sexual_maturity_years",
        "data_quality",
        "sample_size",
        "specimen_origin",
        "source",
        "missing_body_mass",
        "missing_max_lifespan",
        "missing_sexual_maturity",
        "manual_supplement_status",
        "manual_supplement_source",
        "taxonomy_issue",
        "suggested_action",
        "near_anage_candidates",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in species_rows:
            name = row["scientific_name"]
            matched_name = row.get("anage_matched_name", "") or aliases.get(
                name.lower(), name
            )
            anage = anage_by_name.get(matched_name.lower())
            supplement = supplements.get(name.lower())
            issue, action = taxonomy_issue(row, anage, aliases, supplement)
            writer.writerow(
                {
                    "scientific_name": name,
                    "anage_matched_name": row.get("anage_matched_name", ""),
                    "alias_applied": "yes" if name.lower() in aliases else "no",
                    "anage_record_present": row.get("anage_record_present", ""),
                    "hagrid": "" if not anage else anage.get("HAGRID", ""),
                    "anage_class": "" if not anage else anage.get("Class", ""),
                    "anage_order": "" if not anage else anage.get("Order", ""),
                    "anage_family": "" if not anage else anage.get("Family", ""),
                    "anage_common_name": "" if not anage else anage.get("Common name", ""),
                    "clade": row.get("clade", ""),
                    "flight_status": row.get("flight_status", ""),
                    "body_mass_g": row.get("body_mass_g", ""),
                    "max_lifespan_years": row.get("max_lifespan_years", ""),
                    "sexual_maturity_years": row.get("sexual_maturity_years", ""),
                    "data_quality": "" if not anage else anage.get("Data quality", ""),
                    "sample_size": "" if not anage else anage.get("Sample size", ""),
                    "specimen_origin": "" if not anage else anage.get("Specimen origin", ""),
                    "source": "" if not anage else anage.get("Source", ""),
                    "missing_body_mass": "yes" if not row.get("body_mass_g", "").strip() else "no",
                    "missing_max_lifespan": "yes" if not row.get("max_lifespan_years", "").strip() else "no",
                    "missing_sexual_maturity": "yes" if not row.get("sexual_maturity_years", "").strip() else "no",
                    "manual_supplement_status": ""
                    if not supplement
                    else supplement.get("use_in_model", ""),
                    "manual_supplement_source": ""
                    if not supplement
                    else supplement.get("source_url", "") or supplement.get("source_name", ""),
                    "taxonomy_issue": issue,
                    "suggested_action": action,
                    "near_anage_candidates": ""
                    if anage
                    else candidate_names(name, anage_rows),
                }
            )


def zscore(value: float | None, values: list[float]) -> float | None:
    if value is None or len(values) < 2:
        return None
    sd = statistics.stdev(values)
    if sd == 0:
        return None
    return (value - statistics.fmean(values)) / sd


def outlier_flags(row: dict[str, str], z_mass: float | None, z_life: float | None, z_resid: float | None) -> str:
    flags = []
    residual = parse_float(row.get("lifespan_residual_log10"))
    ratio = parse_float(row.get("lifespan_residual_ratio"))
    if z_resid is not None and z_resid >= 2:
        flags.append("high_positive_residual_z")
    if z_resid is not None and z_resid <= -2:
        flags.append("high_negative_residual_z")
    if ratio is not None and ratio >= 3:
        flags.append("lifespan_3x_above_body_mass_expectation")
    if ratio is not None and ratio <= 0.33:
        flags.append("lifespan_3x_below_body_mass_expectation")
    if z_mass is not None and abs(z_mass) >= 2:
        flags.append("body_mass_extreme")
    if z_life is not None and abs(z_life) >= 2:
        flags.append("lifespan_extreme")
    if residual is not None and abs(residual) >= 0.5:
        flags.append("large_log10_residual")
    return ";".join(flags) if flags else "none"


def write_outlier_audit(
    residual_rows: list[dict[str, str]],
    taxonomy_rows: list[dict[str, str]],
    output: pathlib.Path,
) -> None:
    taxonomy_by_name = {row["scientific_name"]: row for row in taxonomy_rows}
    all_log_mass = [parse_float(row.get("log10_body_mass_g")) for row in residual_rows]
    all_log_life = [parse_float(row.get("log10_max_lifespan_years")) for row in residual_rows]
    all_resid = [parse_float(row.get("lifespan_residual_log10")) for row in residual_rows]
    all_log_mass = [value for value in all_log_mass if value is not None]
    all_log_life = [value for value in all_log_life if value is not None]
    all_resid = [value for value in all_resid if value is not None]

    by_clade: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in residual_rows:
        clade = row.get("clade", "")
        for key in [
            "log10_body_mass_g",
            "log10_max_lifespan_years",
            "lifespan_residual_log10",
        ]:
            value = parse_float(row.get(key))
            if value is not None:
                by_clade[clade][key].append(value)

    fields = [
        "rank_abs_residual",
        "scientific_name",
        "clade",
        "flight_status",
        "body_mass_g",
        "max_lifespan_years",
        "lifespan_residual_log10",
        "lifespan_residual_ratio",
        "z_log_body_mass_all",
        "z_log_lifespan_all",
        "z_residual_all",
        "z_residual_within_clade",
        "data_quality",
        "sample_size",
        "specimen_origin",
        "source",
        "outlier_flags",
        "audit_priority",
    ]

    enriched = []
    for row in residual_rows:
        tax = taxonomy_by_name.get(row["scientific_name"], {})
        log_mass = parse_float(row.get("log10_body_mass_g"))
        log_life = parse_float(row.get("log10_max_lifespan_years"))
        residual = parse_float(row.get("lifespan_residual_log10"))
        z_mass = zscore(log_mass, all_log_mass)
        z_life = zscore(log_life, all_log_life)
        z_resid = zscore(residual, all_resid)
        clade_resid = by_clade[row.get("clade", "")]["lifespan_residual_log10"]
        z_clade_resid = zscore(residual, clade_resid)
        flags = outlier_flags(row, z_mass, z_life, z_resid)
        priority = "high" if flags != "none" or tax.get("data_quality") == "questionable" else "standard"
        enriched.append(
            {
                "rank_abs_residual": "",
                "scientific_name": row["scientific_name"],
                "clade": row.get("clade", ""),
                "flight_status": row.get("flight_status", ""),
                "body_mass_g": row.get("body_mass_g", ""),
                "max_lifespan_years": row.get("max_lifespan_years", ""),
                "lifespan_residual_log10": row.get("lifespan_residual_log10", ""),
                "lifespan_residual_ratio": row.get("lifespan_residual_ratio", ""),
                "z_log_body_mass_all": "" if z_mass is None else f"{z_mass:.6f}",
                "z_log_lifespan_all": "" if z_life is None else f"{z_life:.6f}",
                "z_residual_all": "" if z_resid is None else f"{z_resid:.6f}",
                "z_residual_within_clade": ""
                if z_clade_resid is None
                else f"{z_clade_resid:.6f}",
                "data_quality": tax.get("data_quality", ""),
                "sample_size": tax.get("sample_size", ""),
                "specimen_origin": tax.get("specimen_origin", ""),
                "source": tax.get("source", ""),
                "outlier_flags": flags,
                "audit_priority": priority,
            }
        )

    enriched.sort(
        key=lambda row: abs(parse_float(row["lifespan_residual_log10"]) or 0),
        reverse=True,
    )
    for index, row in enumerate(enriched, start=1):
        row["rank_abs_residual"] = str(index)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(enriched)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-master", type=pathlib.Path, default=DEFAULT_SPECIES_MASTER)
    parser.add_argument("--residuals", type=pathlib.Path, default=DEFAULT_RESIDUALS)
    parser.add_argument("--anage", type=pathlib.Path, default=DEFAULT_ANAGE)
    parser.add_argument("--aliases", type=pathlib.Path, default=DEFAULT_ALIASES)
    parser.add_argument("--supplements", type=pathlib.Path, default=DEFAULT_SUPPLEMENTS)
    parser.add_argument("--taxonomy-output", type=pathlib.Path, default=DEFAULT_TAXONOMY_OUTPUT)
    parser.add_argument("--outlier-output", type=pathlib.Path, default=DEFAULT_OUTLIER_OUTPUT)
    args = parser.parse_args()

    species_rows = read_tsv(args.species_master)
    residual_rows = read_tsv(args.residuals)
    anage_rows = read_tsv(args.anage)
    aliases = read_aliases(args.aliases)
    supplements = read_supplements(args.supplements)

    write_taxonomy_audit(
        species_rows,
        anage_rows,
        aliases,
        supplements,
        args.taxonomy_output,
    )
    taxonomy_rows = read_tsv(args.taxonomy_output)
    write_outlier_audit(residual_rows, taxonomy_rows, args.outlier_output)
    print(f"Wrote {args.taxonomy_output} and {args.outlier_output}")


if __name__ == "__main__":
    main()
