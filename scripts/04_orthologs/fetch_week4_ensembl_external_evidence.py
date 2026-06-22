"""Fetch Ensembl REST evidence for rows where target species are available."""

from __future__ import annotations

import argparse
import json
import pathlib
import time
import urllib.parse
import urllib.error
import urllib.request

import pandas as pd


DEFAULT_PLAN = pathlib.Path("data/processed/week4_external_orthology_query_plan.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/week4_ensembl_external_evidence.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week4_ensembl_external_evidence_report.md")


def fetch_json(url: str) -> tuple[int, object | None, str]:
    if not isinstance(url, str) or not url:
        return 0, None, "missing_url"
    request = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            text = response.read().decode("utf-8", errors="replace")
            return response.status, json.loads(text), ""
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace") if exc.fp else ""
        return exc.code, None, body[:300] or exc.reason
    except Exception as exc:  # noqa: BLE001
        return 0, None, str(exc)


def summarize_xrefs(payload: object | None) -> tuple[int, str, str]:
    if not isinstance(payload, list):
        return 0, "", ""
    ids = []
    descriptions = []
    for item in payload:
        if isinstance(item, dict):
            ids.append(str(item.get("id", "")))
            descriptions.append(str(item.get("description", "")))
    return len(payload), ";".join(ids[:5]), ";".join(descriptions[:3])


def summarize_homology(payload: object | None) -> tuple[int, int, str, str]:
    if not isinstance(payload, dict):
        return 0, 0, "", ""
    data = payload.get("data", [])
    if not data:
        return 0, 0, "", ""
    homologies = data[0].get("homologies", []) if isinstance(data[0], dict) else []
    orthologue_count = 0
    species = []
    types = []
    for hom in homologies:
        if not isinstance(hom, dict):
            continue
        hom_type = str(hom.get("type", ""))
        target = hom.get("target", {})
        if "ortholog" in hom_type.lower() or "orthologue" in hom_type.lower():
            orthologue_count += 1
        types.append(hom_type)
        if isinstance(target, dict):
            species.append(str(target.get("species", "")))
    return len(homologies), orthologue_count, ";".join(sorted(set(species))[:10]), ";".join(sorted(set(types))[:10])


def summarize_accession_xrefs(payload: object | None) -> tuple[int, str, str, str]:
    if not isinstance(payload, list):
        return 0, "", "", ""
    ids = []
    species = []
    dbnames = []
    for item in payload:
        if isinstance(item, dict):
            ids.append(str(item.get("id", "")))
            species.append(str(item.get("species", "")))
            dbnames.append(str(item.get("dbname", "")))
    return len(payload), ";".join(ids[:5]), ";".join(sorted(set(species))[:5]), ";".join(sorted(set(dbnames))[:5])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=pathlib.Path, default=DEFAULT_PLAN)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    plan = pd.read_csv(args.plan, sep="\t")
    rows = plan[plan["ensembl_available"] == "yes"].copy()
    out = []
    for _, row in rows.iterrows():
        x_status, x_payload, x_error = fetch_json(row.get("ensembl_xrefs_symbol_url", ""))
        time.sleep(0.15)
        h_status, h_payload, h_error = fetch_json(row.get("ensembl_homology_symbol_url", ""))
        time.sleep(0.15)
        accession = str(row.get("forward_target_protein_id", ""))
        accession_url = (
            "https://rest.ensembl.org/xrefs/id/"
            f"{urllib.parse.quote(accession)}?content-type=application/json"
            if accession
            else ""
        )
        a_status, a_payload, a_error = fetch_json(accession_url)
        time.sleep(0.15)
        x_count, x_ids, x_desc = summarize_xrefs(x_payload)
        h_count, h_orth, h_species, h_types = summarize_homology(h_payload)
        a_count, a_ids, a_species, a_dbnames = summarize_accession_xrefs(a_payload)
        out.append(
            {
                "human_gene_symbol": row["human_gene_symbol"],
                "scientific_name": row["scientific_name"],
                "ensembl_species": row["ensembl_species"],
                "week4_diamond_status": row["week4_diamond_status"],
                "week4_diamond_reason": row["week4_diamond_reason"],
                "forward_target_protein_id": row["forward_target_protein_id"],
                "reciprocal_gene": row["reciprocal_gene"],
                "xrefs_http_status": x_status,
                "xrefs_count": x_count,
                "xrefs_ids": x_ids,
                "xrefs_descriptions": x_desc,
                "xrefs_error": x_error,
                "homology_http_status": h_status,
                "homology_count": h_count,
                "homology_orthologue_count": h_orth,
                "homology_target_species_sample": h_species,
                "homology_types_sample": h_types,
                "homology_error": h_error,
                "accession_xrefs_http_status": a_status,
                "accession_xrefs_count": a_count,
                "accession_xrefs_ids": a_ids,
                "accession_xrefs_species": a_species,
                "accession_xrefs_dbnames": a_dbnames,
                "accession_xrefs_error": a_error,
                "ensembl_support_summary": (
                    "symbol_and_accession_found"
                    if x_count > 0 and a_count > 0
                    else "accession_found"
                    if a_count > 0
                    else "symbol_found"
                    if x_count > 0
                    else "no_symbol_or_accession"
                ),
            }
        )
    evidence = pd.DataFrame(out)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(args.output, sep="\t", index=False)

    symbol_counts = evidence["ensembl_support_summary"].value_counts().sort_index() if not evidence.empty else {}
    lines = [
        "# Week 4 Ensembl External Evidence Report",
        "",
        f"Rows queried: {len(evidence)}",
        "",
        "## Ensembl Evidence",
        "",
    ]
    if len(evidence):
        for status, count in symbol_counts.items():
            lines.append(f"- {status}: {count}")
    else:
        lines.append("- No Ensembl-available rows")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This is a narrow Ensembl REST evidence pass for queue species recognized by Ensembl. It records symbol-level availability and homology metadata, but ambiguous PIWI/TRIM rows still require manual orthology interpretation.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
