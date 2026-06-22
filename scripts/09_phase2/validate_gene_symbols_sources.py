"""Validate expanded-panel gene symbols and attach lightweight source evidence."""

from __future__ import annotations

import argparse
import json
import pathlib
import time
from typing import Any

import pandas as pd
import requests


MYGENE_URL = "https://mygene.info/v3/query"

SOURCE_KEYWORDS = {
    "transposon_repeat_suppression": [
        "piRNA",
        "transposon",
        "retrotransposon",
        "LINE-1",
        "repeat",
        "RNA interference",
        "retroelement",
    ],
    "chromatin_repression_heterochromatin": [
        "chromatin",
        "heterochromatin",
        "histone",
        "methyltransferase",
        "polycomb",
        "silencing",
    ],
    "DNA_repair_replication_stress": [
        "DNA repair",
        "double-strand break",
        "homologous recombination",
        "nucleotide excision",
        "base excision",
        "checkpoint",
    ],
    "proteostasis_autophagy_mitophagy": [
        "chaperone",
        "proteasome",
        "ubiquitin",
        "autophagy",
        "mitophagy",
        "mitochondrial",
    ],
    "inflammation_innate_immune_restraint": [
        "inflammatory",
        "inflammasome",
        "interferon",
        "innate immune",
        "NF-kappa",
        "cytosolic DNA",
    ],
    "cancer_surveillance_senescence": [
        "tumor suppressor",
        "cell cycle",
        "senescence",
        "apoptosis",
        "telomere",
        "genome stability",
    ],
}


MANUAL_ALIASES = {
    "STING1": "TMEM173",
    "PARK7": "PARK7",
    "PRKN": "PRKN",
}


def query_mygene(symbol: str, cache_dir: pathlib.Path, sleep_s: float) -> dict[str, Any]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{symbol}.json"
    if cache_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))

    params = {
        "q": f"symbol:{symbol} OR alias:{symbol}",
        "species": "human",
        "fields": "symbol,name,alias,summary,entrezgene,HGNC,go.BP,go.MF,go.CC,uniprot,ensembl.gene,type_of_gene,taxid",
        "size": 5,
    }
    response = requests.get(MYGENE_URL, params=params, timeout=30)
    response.raise_for_status()
    payload = response.json()
    cache_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    if sleep_s:
        time.sleep(sleep_s)
    return payload


def best_hit(symbol: str, payload: dict[str, Any]) -> dict[str, Any] | None:
    hits = payload.get("hits", [])
    if not hits:
        return None
    exact = [h for h in hits if str(h.get("symbol", "")).upper() == symbol.upper()]
    if exact:
        return exact[0]
    alias_matches = []
    for hit in hits:
        aliases = hit.get("alias", [])
        if isinstance(aliases, str):
            aliases = [aliases]
        if symbol.upper() in {str(a).upper() for a in aliases}:
            alias_matches.append(hit)
    if alias_matches:
        return alias_matches[0]
    return hits[0]


def flatten_go_terms(hit: dict[str, Any]) -> str:
    terms = []
    go_bp = hit.get("go", {}).get("BP", []) if isinstance(hit.get("go"), dict) else []
    if isinstance(go_bp, dict):
        go_bp = [go_bp]
    for term in go_bp[:12]:
        label = term.get("term")
        if label:
            terms.append(label)
    return ";".join(dict.fromkeys(terms))


def evidence_tags(module: str, summary: str, go_terms: str, submodule: str) -> str:
    text = f"{summary} {go_terms} {submodule}".lower()
    tags = []
    for keyword in SOURCE_KEYWORDS.get(module, []):
        if keyword.lower() in text:
            tags.append(keyword)
    if submodule:
        tags.append(f"submodule:{submodule}")
    return ";".join(dict.fromkeys(tags))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--review-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    parser.add_argument("--cache-dir", type=pathlib.Path, default=pathlib.Path("data/interim/phase2/mygene"))
    parser.add_argument("--sleep", type=float, default=0.05)
    args = parser.parse_args()

    panel = pd.read_csv(args.panel, sep="\t")
    rows = []
    for _, row in panel.iterrows():
        symbol = row["human_gene_symbol"]
        query_symbol = MANUAL_ALIASES.get(symbol, symbol)
        payload = query_mygene(query_symbol, args.cache_dir, args.sleep)
        hit = best_hit(query_symbol, payload)
        if hit is None:
            current_symbol = ""
            validation_status = "not_found"
            name = ""
            entrez = ""
            summary = ""
            go_terms = ""
            aliases = ""
        else:
            current_symbol = hit.get("symbol", "")
            validation_status = (
                "validated_exact"
                if current_symbol.upper() == query_symbol.upper()
                else "validated_alias_or_updated"
            )
            name = hit.get("name", "")
            entrez = hit.get("entrezgene", "")
            summary = hit.get("summary", "")
            go_terms = flatten_go_terms(hit)
            hit_alias = hit.get("alias", [])
            if isinstance(hit_alias, str):
                hit_alias = [hit_alias]
            aliases = ";".join(str(a) for a in hit_alias[:20])

        tags = evidence_tags(
            row["maintenance_module_v2"],
            summary,
            go_terms,
            row["submodule_v2"],
        )
        source_status = "has_source_evidence" if tags else "needs_source_review"
        rows.append(
            {
                **row.to_dict(),
                "query_symbol": symbol,
                "mygene_query_symbol": query_symbol,
                "validated_symbol": current_symbol,
                "symbol_validation_status": validation_status,
                "entrezgene": entrez,
                "gene_name": name,
                "mygene_aliases": aliases,
                "go_bp_terms_top": go_terms,
                "source_evidence_tags": tags,
                "source_evidence_status": source_status,
                "validation_decision": (
                    "pass"
                    if validation_status.startswith("validated")
                    and source_status == "has_source_evidence"
                    else "review"
                ),
            }
        )

    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    review = out[out["validation_decision"] != "pass"].copy()
    args.review_output.parent.mkdir(parents=True, exist_ok=True)
    review.to_csv(args.review_output, sep="\t", index=False)

    n = len(out)
    n_symbol = int(out["symbol_validation_status"].str.startswith("validated").sum())
    n_source = int((out["source_evidence_status"] == "has_source_evidence").sum())
    n_pass = int((out["validation_decision"] == "pass").sum())
    gate = "pass" if n_pass / n >= 0.90 else ("caution" if n_pass / n >= 0.75 else "fail")
    lines = [
        "# Phase 2 P2.1 Gene Symbol and Source Validation Report",
        "",
        f"Input genes: {n}",
        f"Validated current/alias symbols: {n_symbol} ({n_symbol / n:.1%})",
        f"Genes with source evidence tags: {n_source} ({n_source / n:.1%})",
        f"Genes passing both checks: {n_pass} ({n_pass / n:.1%})",
        "",
        f"Decision gate: **{gate}**",
        "",
        "Gate definition: pass requires at least 90 percent of genes to have validated symbols and at least one source evidence tag.",
        "",
        f"Review rows: {len(review)}",
    ]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} with gate={gate}")


if __name__ == "__main__":
    main()
