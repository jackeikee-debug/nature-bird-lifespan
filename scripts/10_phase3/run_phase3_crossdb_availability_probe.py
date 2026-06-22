"""Probe OMA, OrthoDB, and Ensembl availability for Phase 3 priority-1 rows."""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd


OMA_XREF = "https://omabrowser.org/api/xref/"
ORTHODB_SEARCH = "https://data.orthodb.org/current/search"
ENSEMBL_SPECIES = "https://rest.ensembl.org/info/species"
ENSEMBL_HOMOLOGY = "https://rest.ensembl.org/homology/symbol/human/{symbol}"

PREFERRED_ORTHODB_LEVELS = ["Aves", "Sauropsida", "Tetrapoda", "Vertebrata"]


def safe_key(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(value))


def urlopen_text(url: str, headers: dict[str, str] | None = None, timeout: int = 60) -> str:
    request = urllib.request.Request(url, headers=headers or {})
    last_error: Exception | None = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code in {429, 500, 502, 503, 504}:
                time.sleep(4 * (attempt + 1))
                continue
            return json.dumps({"_http_error": exc.code, "_message": str(exc)})
        except urllib.error.URLError as exc:
            last_error = exc
            time.sleep(3 * (attempt + 1))
    return json.dumps({"_error": str(last_error)})


def cached_json(url: str, cache_path: pathlib.Path, headers: dict[str, str] | None = None, sleep: float = 0.0):
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists():
        text = cache_path.read_text(encoding="utf-8")
    else:
        text = urlopen_text(url, headers=headers)
        cache_path.write_text(text, encoding="utf-8")
        if sleep:
            time.sleep(sleep)
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return {"_raw": text}


def norm_symbol(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(text).lower())


def oma_gene_hits(symbol: str, cache_dir: pathlib.Path, sleep: float) -> list[dict]:
    url = OMA_XREF + "?" + urllib.parse.urlencode({"search": symbol})
    data = cached_json(url, cache_dir / "oma" / f"xref_{safe_key(symbol)}.json", sleep=sleep)
    return data if isinstance(data, list) else []


def classify_oma(row: pd.Series, gene_hits: list[dict]) -> tuple[str, str, str, str, str]:
    taxid = str(row["species_taxid"])
    symbol = row["human_gene_symbol"]
    target = [hit for hit in gene_hits if str(hit.get("genome", {}).get("taxon_id", "")) == taxid]
    if not target:
        return "oma_no_target_taxid_hit", "OMA xref returned no hit for the target species taxid.", "", "", ""
    exact = []
    for hit in target:
        xref = str(hit.get("xref", ""))
        if norm_symbol(xref.split("{", 1)[0]) == norm_symbol(symbol) and str(hit.get("seq_match", "")) == "exact":
            exact.append(hit)
    best = exact[0] if exact else target[0]
    call = "oma_target_exact_symbol_support" if exact else "oma_target_broad_hit_manual_review"
    reason = "OMA xref contains exact symbol hit for target taxid." if exact else "OMA target taxid hit exists but exact symbol rule did not match."
    return (
        call,
        reason,
        str(best.get("omaid", "")),
        str(best.get("xref", "")),
        str(best.get("genome", {}).get("species", "")),
    )


def orthodb_gene_support(symbol: str, cache_dir: pathlib.Path, sleep: float) -> tuple[str, str, str, str, str]:
    url = ORTHODB_SEARCH + "?" + urllib.parse.urlencode({"query": symbol})
    data = cached_json(url, cache_dir / "orthodb" / f"search_{safe_key(symbol)}.json", sleep=sleep)
    groups = data.get("bigdata", []) if isinstance(data, dict) else []
    if not groups:
        return "orthodb_no_search_hit", "OrthoDB search returned no orthogroups.", "", "", ""
    symbol_norm = norm_symbol(symbol)
    candidates = []
    for group in groups:
        name = str(group.get("name", ""))
        level = str(group.get("level_name", ""))
        if symbol_norm in norm_symbol(name) or symbol.lower() in name.lower():
            candidates.append(group)
    if not candidates:
        candidates = groups[:10]
    for level in PREFERRED_ORTHODB_LEVELS:
        for group in candidates:
            if str(group.get("level_name", "")) == level:
                return (
                    "orthodb_preferred_level_support",
                    f"OrthoDB search found a {level}-level orthogroup for this symbol/name.",
                    str(group.get("public_id", group.get("id", ""))),
                    str(group.get("name", "")),
                    level,
                )
    group = candidates[0]
    return (
        "orthodb_broad_group_support",
        "OrthoDB search found a broad orthogroup but not at preferred bird/reptile/vertebrate levels.",
        str(group.get("public_id", group.get("id", ""))),
        str(group.get("name", "")),
        str(group.get("level_name", "")),
    )


def ensembl_species_map(cache_dir: pathlib.Path, sleep: float) -> dict[str, dict]:
    url = ENSEMBL_SPECIES + "?" + urllib.parse.urlencode({"content-type": "application/json"})
    data = cached_json(url, cache_dir / "ensembl" / "species.json", sleep=sleep)
    species = data.get("species", []) if isinstance(data, dict) else []
    return {str(item.get("taxon_id", "")): item for item in species}


def ensembl_homology(symbol: str, taxid: str, cache_dir: pathlib.Path, sleep: float) -> tuple[str, str, int]:
    params = {"type": "orthologues", "target_taxon": taxid, "format": "condensed"}
    url = ENSEMBL_HOMOLOGY.format(symbol=urllib.parse.quote(symbol)) + "?" + urllib.parse.urlencode(params)
    data = cached_json(
        url,
        cache_dir / "ensembl" / f"homology_{safe_key(symbol)}_{safe_key(taxid)}.json",
        headers={"Content-Type": "application/json"},
        sleep=sleep,
    )
    if isinstance(data, dict) and data.get("_http_error"):
        return "ensembl_homology_http_error", str(data.get("_message", "")), 0
    homologies = []
    for item in data.get("data", []) if isinstance(data, dict) else []:
        homologies.extend(item.get("homologies", []))
    if homologies:
        return "ensembl_target_homology_support", "Ensembl REST returned target_taxon homology records.", len(homologies)
    return "ensembl_no_target_homology", "Ensembl REST returned no target_taxon homology records.", 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--cache-dir", type=pathlib.Path, required=True)
    parser.add_argument("--sleep", type=float, default=0.05)
    args = parser.parse_args()

    audit = pd.read_csv(args.audit, sep="\t", dtype=str).fillna("")
    genes = sorted(audit["human_gene_symbol"].unique())
    oma_hits = {gene: oma_gene_hits(gene, args.cache_dir, args.sleep) for gene in genes}
    orthodb_hits = {gene: orthodb_gene_support(gene, args.cache_dir, args.sleep) for gene in genes}
    ensembl_map = ensembl_species_map(args.cache_dir, args.sleep)

    rows = []
    for _, row in audit.iterrows():
        oma_call, oma_reason, oma_id, oma_xref, oma_species = classify_oma(row, oma_hits[row["human_gene_symbol"]])
        odb_call, odb_reason, odb_id, odb_name, odb_level = orthodb_hits[row["human_gene_symbol"]]
        taxid = str(row["species_taxid"])
        ensembl_species = ensembl_map.get(taxid, {})
        if ensembl_species:
            ens_species_call = "ensembl_species_present"
            ens_species_name = str(ensembl_species.get("name", ""))
            ens_call, ens_reason, ens_count = ensembl_homology(row["human_gene_symbol"], taxid, args.cache_dir, args.sleep)
        else:
            ens_species_call = "ensembl_species_absent"
            ens_species_name = ""
            ens_call, ens_reason, ens_count = "ensembl_not_queried_species_absent", "Target species taxid is absent from Ensembl species list.", 0
        rows.append(
            {
                "scientific_name": row["scientific_name"],
                "human_gene_symbol": row["human_gene_symbol"],
                "species_taxid": taxid,
                "next_rescue_route": row["next_rescue_route"],
                "phase3_gff_sequence_decision": row["phase3_gff_sequence_decision"],
                "oma_call": oma_call,
                "oma_reason": oma_reason,
                "oma_id_top": oma_id,
                "oma_xref_top": oma_xref,
                "oma_species_top": oma_species,
                "orthodb_call": odb_call,
                "orthodb_reason": odb_reason,
                "orthodb_group_id": odb_id,
                "orthodb_group_name": odb_name,
                "orthodb_level": odb_level,
                "ensembl_species_call": ens_species_call,
                "ensembl_species_name": ens_species_name,
                "ensembl_homology_call": ens_call,
                "ensembl_homology_reason": ens_reason,
                "ensembl_homology_count": ens_count,
            }
        )
    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    summary = (
        out.groupby(["next_rescue_route", "oma_call", "ensembl_species_call"], as_index=False)
        .agg(rows=("human_gene_symbol", "count"), species=("scientific_name", "nunique"), genes=("human_gene_symbol", "nunique"))
        .sort_values(["next_rescue_route", "oma_call", "ensembl_species_call"])
    )
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary_output, sep="\t", index=False)

    oma_support = int((out["oma_call"] == "oma_target_exact_symbol_support").sum())
    ensembl_present = int((out["ensembl_species_call"] == "ensembl_species_present").sum())
    odb_pref = int((out["orthodb_call"] == "orthodb_preferred_level_support").sum())
    if oma_support > 0:
        oma_sentence = "OMA provides species-gene level support for some target rows and should be used for targeted rescue."
    else:
        oma_sentence = "OMA did not provide target-taxid support for this priority-1 set; it is not currently useful for direct species-gene rescue of these rows."
    if ensembl_present > 0:
        ensembl_sentence = "Ensembl covers a small subset of target species, but homology records must be checked row by row before use."
    else:
        ensembl_sentence = "Ensembl target-species coverage is absent for this set, so Compara is not useful for direct rescue here."
    lines = [
        "# Phase 3 Cross-Database Availability Probe Report",
        "",
        f"Rows probed: {len(out)}",
        f"OMA exact target-taxid symbol support rows: {oma_support}",
        f"OrthoDB preferred-level orthogroup support rows: {odb_pref}",
        f"Rows whose target species is present in Ensembl: {ensembl_present}",
        "",
        "## Interpretation",
        f"{oma_sentence} OrthoDB is currently the strongest cross-database layer for these rows, but mainly as gene-family/orthogroup evidence rather than species-gene strict rescue. {ensembl_sentence}",
        "",
        "## Outputs",
        f"- row table: `{args.output}`",
        f"- summary: `{args.summary_output}`",
    ]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
