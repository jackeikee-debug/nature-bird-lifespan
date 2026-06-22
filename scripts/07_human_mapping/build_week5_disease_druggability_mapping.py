"""Build Week 5 disease association and tractability mapping."""

from __future__ import annotations

import argparse
import json
import pathlib
import time
from collections import Counter
from typing import Any

import pandas as pd
import requests


DEFAULT_HUMAN_MAPPING = pathlib.Path("data/processed/human_mapping.tsv")
DEFAULT_RAW_DIR = pathlib.Path("data/raw/open_targets")
DEFAULT_INTERIM = pathlib.Path("data/interim/human_mapping/week5_open_targets_gene_summary.tsv")
DEFAULT_DISEASE_OUTPUT = pathlib.Path("data/processed/human_disease_mapping.tsv")
DEFAULT_GENE_OUTPUT = pathlib.Path("data/processed/human_translation_priority.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week5_disease_druggability_report.md")

MYGENE_URL = "https://mygene.info/v3/query"
OPENTARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

TARGET_QUERY = """
query targetDiseaseAssociations($ensemblId: String!, $page: Pagination) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    associatedDiseases(page: $page) {
      count
      rows {
        score
        disease {
          id
          name
          therapeuticAreas { name }
        }
        datatypeScores { id score }
        datasourceScores { id score }
      }
    }
    tractability { modality label value }
  }
}
"""

DISEASE_AREA_KEYWORDS = {
    "cancer": ["cancer", "carcinoma", "neoplasm", "tumor", "tumour", "leukemia", "lymphoma", "melanoma"],
    "nervous_system": ["nervous", "neuro", "parkinson", "alzheimer", "ataxia"],
    "immune_inflammation": ["immune", "inflammatory", "inflammation", "autoimmune"],
    "reproductive_or_breast": ["reproductive", "breast", "fertility", "gonadal"],
    "genetic_congenital": ["genetic", "familial", "congenital", "orphan"],
    "metabolic_endocrine": ["metabolic", "endocrine", "diabetes"],
}

TRANSLATION_MODULE_BONUS = {
    "transposon_suppression": 2,
    "DNA_repair": 2,
    "cancer_surveillance": 2,
    "inflammation_control": 1,
    "mitochondrial_quality_control": 1,
    "proteostasis": 1,
    "autophagy": 1,
}


def flatten(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, (int, float)):
        return str(value)
    if isinstance(value, list):
        return "; ".join(flatten(v) for v in value if v is not None)
    if isinstance(value, dict):
        return "; ".join(flatten(v) for v in value.values() if v is not None)
    return str(value)


def query_mygene_ensembl(symbol: str) -> dict[str, str]:
    response = requests.get(
        MYGENE_URL,
        params={
            "q": f"symbol:{symbol}",
            "species": "human",
            "fields": "symbol,entrezgene,ensembl.gene",
            "size": 5,
        },
        timeout=60,
    )
    response.raise_for_status()
    hits = response.json().get("hits", [])
    exact = [hit for hit in hits if str(hit.get("symbol", "")).upper() == symbol.upper()]
    hit = exact[0] if exact else hits[0] if hits else {}
    ensembl_value = hit.get("ensembl", {})
    if isinstance(ensembl_value, list):
        ensembl_ids = [str(x.get("gene", "")) for x in ensembl_value if isinstance(x, dict)]
    elif isinstance(ensembl_value, dict):
        ensembl_ids = [str(ensembl_value.get("gene", ""))]
    else:
        ensembl_ids = []
    ensembl_ids = [x for x in ensembl_ids if x.startswith("ENSG")]
    return {
        "symbol": str(hit.get("symbol", "")),
        "entrezgene": str(hit.get("entrezgene", "")),
        "ensembl_gene_id": ensembl_ids[0] if ensembl_ids else "",
        "ensembl_gene_ids_all": ";".join(sorted(set(ensembl_ids))),
    }


def query_open_targets(ensembl_id: str, cache_path: pathlib.Path, top_n: int) -> dict[str, Any]:
    if cache_path.exists() and cache_path.stat().st_size > 0:
        return json.loads(cache_path.read_text(encoding="utf-8"))
    response = requests.post(
        OPENTARGETS_GRAPHQL,
        json={
            "query": TARGET_QUERY,
            "variables": {
                "ensemblId": ensembl_id,
                "page": {"index": 0, "size": top_n},
            },
        },
        timeout=90,
    )
    response.raise_for_status()
    data = response.json()
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return data


def score_pairs(items: list[dict[str, Any]]) -> str:
    pairs = []
    for item in items or []:
        key = item.get("id", "")
        score = item.get("score", "")
        if key != "" and score != "":
            pairs.append(f"{key}:{float(score):.3f}")
    return ";".join(pairs)


def tractability_summary(items: list[dict[str, Any]]) -> tuple[str, str]:
    positives = []
    labels = []
    for item in items or []:
        label = str(item.get("label", ""))
        modality = str(item.get("modality", ""))
        value = item.get("value", False)
        if label:
            labels.append(f"{modality}:{label}:{value}")
        if value is True:
            positives.append(f"{modality}:{label}")
    return ";".join(positives), ";".join(labels)


def small_molecule_tractability(tractability_positive: str) -> str:
    parts = [part for part in str(tractability_positive).split(";") if part]
    sm_parts = [part for part in parts if part.startswith("SM:")]
    return ";".join(sm_parts)


def classify_area(disease_name: str, therapeutic_areas: str) -> str:
    text = f"{disease_name} {therapeutic_areas}".lower()
    hits = []
    for area, keywords in DISEASE_AREA_KEYWORDS.items():
        if any(keyword in text for keyword in keywords):
            hits.append(area)
    return ";".join(hits) if hits else "other"


def gene_translation_priority(row: pd.Series) -> str:
    score = 0
    score += TRANSLATION_MODULE_BONUS.get(row["maintenance_module"], 0)
    if row["max_open_targets_score"] >= 0.70:
        score += 3
    elif row["max_open_targets_score"] >= 0.40:
        score += 2
    elif row["max_open_targets_score"] >= 0.20:
        score += 1
    if row["small_molecule_tractability_count"] > 0:
        score += 2
    if row["cancer_top_disease_count"] > 0:
        score += 1
    if row["immune_inflammation_top_disease_count"] > 0:
        score += 1
    if row["nervous_system_top_disease_count"] > 0:
        score += 1
    if "high_transposon" in str(row.get("week5_v1_priority", "")):
        score += 1
    if score >= 7:
        return "high_translation_priority"
    if score >= 5:
        return "medium_translation_priority"
    if score >= 3:
        return "supporting_translation_context"
    return "low_translation_priority"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--human-mapping", type=pathlib.Path, default=DEFAULT_HUMAN_MAPPING)
    parser.add_argument("--raw-dir", type=pathlib.Path, default=DEFAULT_RAW_DIR)
    parser.add_argument("--interim", type=pathlib.Path, default=DEFAULT_INTERIM)
    parser.add_argument("--disease-output", type=pathlib.Path, default=DEFAULT_DISEASE_OUTPUT)
    parser.add_argument("--gene-output", type=pathlib.Path, default=DEFAULT_GENE_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    parser.add_argument("--top-n", type=int, default=10)
    args = parser.parse_args()

    human = pd.read_csv(args.human_mapping, sep="\t")
    args.raw_dir.mkdir(parents=True, exist_ok=True)
    args.interim.parent.mkdir(parents=True, exist_ok=True)

    gene_meta_rows = []
    disease_rows = []
    failures = []
    for _, gene in human.iterrows():
        symbol = str(gene["human_gene_symbol"]).upper()
        try:
            meta = query_mygene_ensembl(symbol)
            time.sleep(0.05)
            ensembl_candidates = [
                x for x in str(meta.get("ensembl_gene_ids_all", "")).split(";")
                if x.startswith("ENSG")
            ]
            ensembl_id = meta.get("ensembl_gene_id", "")
            if ensembl_id and ensembl_id not in ensembl_candidates:
                ensembl_candidates.insert(0, ensembl_id)
            gene_meta = {
                "human_gene_symbol": symbol,
                "maintenance_module": gene["maintenance_module"],
                "week5_v1_priority": gene["interpretation_priority"],
                **meta,
            }
            if not ensembl_candidates:
                gene_meta["open_targets_status"] = "no_ensembl_id"
                gene_meta_rows.append(gene_meta)
                continue
            target = None
            raw = {}
            tried_ids = []
            for candidate_id in ensembl_candidates:
                tried_ids.append(candidate_id)
                raw = query_open_targets(candidate_id, args.raw_dir / f"{symbol}__{candidate_id}.json", args.top_n)
                target = raw.get("data", {}).get("target")
                if target:
                    ensembl_id = candidate_id
                    break
                time.sleep(0.05)
            gene_meta["open_targets_ensembl_ids_tried"] = ";".join(tried_ids)
            if not target:
                gene_meta["open_targets_status"] = "no_target"
                gene_meta_rows.append(gene_meta)
                continue
            assoc = target.get("associatedDiseases") or {}
            tract_pos, tract_all = tractability_summary(target.get("tractability") or [])
            sm_tract = small_molecule_tractability(tract_pos)
            rows = assoc.get("rows") or []
            gene_meta.update(
                {
                    "open_targets_status": "ok",
                    "ensembl_gene_id": ensembl_id,
                    "open_targets_target_id": target.get("id", ""),
                    "open_targets_approved_symbol": target.get("approvedSymbol", ""),
                    "open_targets_disease_count": assoc.get("count", 0),
                    "tractability_positive": tract_pos,
                    "small_molecule_tractability": sm_tract,
                    "tractability_all": tract_all,
                }
            )
            gene_meta_rows.append(gene_meta)
            for rank, disease_row in enumerate(rows, start=1):
                disease = disease_row.get("disease") or {}
                therapeutic_areas = ";".join(
                    area.get("name", "")
                    for area in disease.get("therapeuticAreas") or []
                    if area.get("name", "")
                )
                disease_name = str(disease.get("name", ""))
                disease_rows.append(
                    {
                        "human_gene_symbol": symbol,
                        "maintenance_module": gene["maintenance_module"],
                        "ensembl_gene_id": ensembl_id,
                        "rank": rank,
                        "disease_id": disease.get("id", ""),
                        "disease_name": disease_name,
                        "therapeutic_areas": therapeutic_areas,
                        "disease_area_class": classify_area(disease_name, therapeutic_areas),
                        "open_targets_score": disease_row.get("score", ""),
                        "datatype_scores": score_pairs(disease_row.get("datatypeScores") or []),
                        "datasource_scores": score_pairs(disease_row.get("datasourceScores") or []),
                        "tractability_positive": tract_pos,
                        "small_molecule_tractability": sm_tract,
                    }
                )
            time.sleep(0.15)
        except Exception as exc:
            failures.append({"human_gene_symbol": symbol, "error": repr(exc)})
            gene_meta_rows.append(
                {
                    "human_gene_symbol": symbol,
                    "maintenance_module": gene["maintenance_module"],
                    "week5_v1_priority": gene["interpretation_priority"],
                    "open_targets_status": "error",
                    "error": repr(exc),
                }
            )

    gene_meta_df = pd.DataFrame(gene_meta_rows)
    disease_df = pd.DataFrame(disease_rows)
    args.interim.parent.mkdir(parents=True, exist_ok=True)
    gene_meta_df.to_csv(args.interim, sep="\t", index=False)
    args.disease_output.parent.mkdir(parents=True, exist_ok=True)
    disease_df.to_csv(args.disease_output, sep="\t", index=False)

    if disease_df.empty:
        summary = gene_meta_df.copy()
        summary["max_open_targets_score"] = 0.0
        summary["top_disease_names"] = ""
    else:
        disease_df["open_targets_score"] = pd.to_numeric(disease_df["open_targets_score"], errors="coerce")
        area_counts = {}
        for area in DISEASE_AREA_KEYWORDS:
            area_counts[area] = disease_df.assign(
                has_area=disease_df["disease_area_class"].str.contains(area, regex=False, na=False)
            ).groupby("human_gene_symbol")["has_area"].sum()
        grouped = disease_df.groupby("human_gene_symbol")
        summary = gene_meta_df.merge(
            grouped.agg(
                max_open_targets_score=("open_targets_score", "max"),
                mean_top_open_targets_score=("open_targets_score", "mean"),
                top_disease_names=("disease_name", lambda x: "; ".join(x.head(5))),
                top_disease_area_classes=("disease_area_class", lambda x: "; ".join(x.head(10))),
            ).reset_index(),
            on="human_gene_symbol",
            how="left",
        )
        for area, series in area_counts.items():
            summary = summary.merge(
                series.rename(f"{area}_top_disease_count").reset_index(),
                on="human_gene_symbol",
                how="left",
            )
    for col in [
        "max_open_targets_score",
        "mean_top_open_targets_score",
        "cancer_top_disease_count",
        "nervous_system_top_disease_count",
        "immune_inflammation_top_disease_count",
        "reproductive_or_breast_top_disease_count",
        "genetic_congenital_top_disease_count",
        "metabolic_endocrine_top_disease_count",
    ]:
        if col not in summary.columns:
            summary[col] = 0
        summary[col] = pd.to_numeric(summary[col], errors="coerce").fillna(0)
    summary["tractability_positive"] = summary.get("tractability_positive", "").fillna("")
    summary["small_molecule_tractability"] = summary.get("small_molecule_tractability", "").fillna("")
    summary["tractability_positive_count"] = summary["tractability_positive"].apply(
        lambda x: 0 if not str(x) else len([part for part in str(x).split(";") if part])
    )
    summary["small_molecule_tractability_count"] = summary["small_molecule_tractability"].apply(
        lambda x: 0 if not str(x) else len([part for part in str(x).split(";") if part])
    )
    summary["week5_translation_priority"] = summary.apply(gene_translation_priority, axis=1)
    priority_rank = {
        "high_translation_priority": 1,
        "medium_translation_priority": 2,
        "supporting_translation_context": 3,
        "low_translation_priority": 4,
    }
    summary["priority_rank"] = summary["week5_translation_priority"].map(priority_rank).fillna(99)
    summary = summary.sort_values(
        ["priority_rank", "max_open_targets_score", "human_gene_symbol"],
        ascending=[True, False, True],
    ).drop(columns=["priority_rank"])
    args.gene_output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.gene_output, sep="\t", index=False)

    priority_counts = Counter(summary["week5_translation_priority"])
    status_counts = Counter(gene_meta_df.get("open_targets_status", []))
    module_priority = (
        summary.groupby(["maintenance_module", "week5_translation_priority"])
        .size()
        .reset_index(name="n")
        .sort_values(["maintenance_module", "week5_translation_priority"])
    )
    trans = summary[summary["maintenance_module"] == "transposon_suppression"]
    report_lines = [
        "# Week 5 Disease and Druggability Mapping Report",
        "",
        "Built Open Targets disease-association and tractability summaries for the 41 human maintenance genes.",
        "",
        "## Sources",
        "",
        "- MyGene.info for human Ensembl gene IDs.",
        "- Open Targets Platform GraphQL API for associated diseases and target tractability.",
        "",
        "## Outputs",
        "",
        f"- Disease rows: {len(disease_df)}",
        f"- Gene summary rows: {len(summary)}",
        f"- Disease output: `{args.disease_output}`",
        f"- Gene output: `{args.gene_output}`",
        f"- Open Targets cache: `{args.raw_dir}`",
        "",
        "## Open Targets Status",
        "",
        *[f"- {key}: {value}" for key, value in sorted(status_counts.items())],
        "",
        "## Translation Priority Counts",
        "",
        *[f"- {key}: {value}" for key, value in sorted(priority_counts.items())],
        "",
        "## Module x Priority",
        "",
        module_priority.to_markdown(index=False),
        "",
        "## Transposon Module Snapshot",
        "",
    ]
    for _, row in trans.iterrows():
        report_lines.append(
            f"- {row['human_gene_symbol']}: max_score={row['max_open_targets_score']:.3f}; "
            f"small_molecule_tractability={row['small_molecule_tractability'] or 'none'}; "
            f"top_diseases={row['top_disease_names']}; "
            f"priority={row['week5_translation_priority']}."
        )
    if failures:
        report_lines.extend(["", "## Failures", ""])
        for failure in failures:
            report_lines.append(f"- {failure['human_gene_symbol']}: {failure['error']}")
    report_lines.extend(
        [
            "",
            "## Caveats",
            "",
            "- Open Targets association scores are disease-target evidence scores, not causal ageing evidence.",
            "- Top disease rows are useful for translational framing and triage, but disease claims should be manually curated before manuscript submission.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.disease_output}, {args.gene_output}, and {args.report}")


if __name__ == "__main__":
    main()
