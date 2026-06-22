"""Build a concrete external-query plan for ambiguous Week 4 orthology rows."""

from __future__ import annotations

import argparse
import json
import pathlib
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd


DEFAULT_QUEUE = pathlib.Path("data/processed/week4_external_orthology_queue.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/week4_external_orthology_query_plan.tsv")
DEFAULT_SPECIES_MAP = pathlib.Path("data/interim/week4_ensembl_species_probe.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week4_external_orthology_query_plan_report.md")

ENSEMBL_REST = "https://rest.ensembl.org"


def ensembl_species_name(scientific_name: str) -> str:
    return scientific_name.lower().replace(" ", "_")


def fetch_json(url: str) -> tuple[int, object | None, str]:
    request = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=45) as response:
            text = response.read().decode("utf-8", errors="replace")
            return response.status, json.loads(text), ""
    except urllib.error.HTTPError as exc:
        return exc.code, None, exc.reason
    except Exception as exc:  # noqa: BLE001
        return 0, None, str(exc)


def probe_ensembl_species(species_names: list[str]) -> pd.DataFrame:
    rows = []
    for scientific_name in species_names:
        species = ensembl_species_name(scientific_name)
        url = f"{ENSEMBL_REST}/info/assembly/{urllib.parse.quote(species)}?"
        status, payload, error = fetch_json(url)
        rows.append(
            {
                "scientific_name": scientific_name,
                "ensembl_species": species,
                "ensembl_probe_url": url,
                "ensembl_http_status": status,
                "ensembl_available": "yes" if status == 200 else "no",
                "ensembl_assembly_name": payload.get("assembly_name", "") if isinstance(payload, dict) else "",
                "ensembl_probe_error": error,
            }
        )
        time.sleep(0.12)
    return pd.DataFrame(rows)


def query_urls(row: pd.Series, ensembl_available: str, ensembl_species: str) -> dict[str, str]:
    gene = row["human_gene_symbol"]
    protein = row.get("forward_target_protein_id", "")
    species = row["scientific_name"]
    urls = {
        "oma_search_url": f"https://omabrowser.org/oma/search/?type=Protein&q={urllib.parse.quote(str(protein or gene))}",
        "orthodb_search_url": f"https://www.orthodb.org/?query={urllib.parse.quote(str(protein or gene))}",
        "ncbi_protein_url": f"https://www.ncbi.nlm.nih.gov/protein/{urllib.parse.quote(str(protein))}" if protein else "",
        "uniprot_search_url": (
            "https://www.uniprot.org/uniprotkb?"
            + urllib.parse.urlencode({"query": f"{gene} {species}"})
        ),
    }
    if ensembl_available == "yes":
        urls["ensembl_xrefs_symbol_url"] = (
            f"{ENSEMBL_REST}/xrefs/symbol/{urllib.parse.quote(ensembl_species)}/"
            f"{urllib.parse.quote(gene)}?content-type=application/json"
        )
        urls["ensembl_homology_symbol_url"] = (
            f"{ENSEMBL_REST}/homology/symbol/{urllib.parse.quote(ensembl_species)}/"
            f"{urllib.parse.quote(gene)}?type=orthologues;sequence=none;content-type=application/json"
        )
    else:
        urls["ensembl_xrefs_symbol_url"] = ""
        urls["ensembl_homology_symbol_url"] = ""
    return urls


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=pathlib.Path, default=DEFAULT_QUEUE)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--species-map", type=pathlib.Path, default=DEFAULT_SPECIES_MAP)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    parser.add_argument("--skip-ensembl-probe", action="store_true")
    args = parser.parse_args()

    queue = pd.read_csv(args.queue, sep="\t")
    species_names = sorted(queue["scientific_name"].unique())
    if args.skip_ensembl_probe and args.species_map.exists():
        species_map = pd.read_csv(args.species_map, sep="\t")
    else:
        species_map = probe_ensembl_species(species_names)
        args.species_map.parent.mkdir(parents=True, exist_ok=True)
        species_map.to_csv(args.species_map, sep="\t", index=False)

    plan = queue.merge(
        species_map[["scientific_name", "ensembl_species", "ensembl_available", "ensembl_assembly_name"]],
        on="scientific_name",
        how="left",
    )
    url_rows = plan.apply(
        lambda row: query_urls(row, row.get("ensembl_available", "no"), row.get("ensembl_species", "")),
        axis=1,
        result_type="expand",
    )
    plan = pd.concat([plan, url_rows], axis=1)
    plan["manual_decision_field"] = ""
    plan["external_database_evidence"] = ""
    plan["recommended_action"] = plan["week4_diamond_status"].map(
        {
            "not_reciprocal": "check_PIWI_family_or_paralogy_before_counting",
            "weak_forward_support": "check_domain_architecture_or_external_orthology_before_counting",
            "reciprocal_weak": "cross_check_before_counting_as_weak_positive",
            "not_validated": "search_external_databases_or_keep_unresolved",
        }
    ).fillna("manual_review")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plan.to_csv(args.output, sep="\t", index=False)

    availability = species_map["ensembl_available"].value_counts().sort_index()
    priority_counts = plan["external_validation_priority"].value_counts().sort_index()
    lines = [
        "# Week 4 External Orthology Query Plan Report",
        "",
        f"Rows planned: {len(plan)}",
        f"Species probed for Ensembl: {len(species_map)}",
        "",
        "## Ensembl Species Probe",
        "",
    ]
    for status, count in availability.items():
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Rows by Priority", ""])
    for priority, count in priority_counts.items():
        lines.append(f"- priority {priority}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This table is the manual/semiautomated worklist for external orthology evidence. Ensembl is only queried for species that the REST assembly endpoint recognizes; OMA, OrthoDB, NCBI Protein, and UniProt search URLs are provided for all rows.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}, {args.species_map}, and {args.report}")


if __name__ == "__main__":
    main()
