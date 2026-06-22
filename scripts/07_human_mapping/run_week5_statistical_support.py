"""Run Week 5 enrichment, background, and permutation tests."""

from __future__ import annotations

import argparse
import io
import json
import math
import pathlib
import random
import re
import time
import zipfile
from collections import Counter
from typing import Any

import numpy as np
import pandas as pd
import requests
from scipy.stats import fisher_exact, mannwhitneyu


DEFAULT_HUMAN = pathlib.Path("data/processed/human_mapping.tsv")
DEFAULT_TRANSLATION = pathlib.Path("data/processed/human_translation_priority.tsv")
DEFAULT_HGNC = pathlib.Path("data/raw/hgnc/hgnc_complete_set.txt")
DEFAULT_BACKGROUND = pathlib.Path("data/processed/week5_background_human_gene_evidence.tsv")
DEFAULT_MODULE_OUT = pathlib.Path("results/tables/week5_module_enrichment_tests.tsv")
DEFAULT_BACKGROUND_OUT = pathlib.Path("results/tables/week5_background_comparison_tests.tsv")
DEFAULT_PERM_OUT = pathlib.Path("results/tables/week5_permutation_tests.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week5_statistical_support_report.md")

HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
MYGENE_URL = "https://mygene.info/v3/query"
OPENTARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
OPENTARGETS_CACHE = pathlib.Path("data/raw/open_targets_background")
MYGENE_CACHE = pathlib.Path("data/raw/mygene_background")

HAGR_ZIPS = {
    "genage": pathlib.Path("data/raw/hagr/genage_human.zip"),
    "longevitymap": pathlib.Path("data/raw/hagr/longevity_map.zip"),
    "cellage": pathlib.Path("data/raw/hagr/cellage.zip"),
}

TARGET_QUERY = """
query targetDiseaseAssociations($ensemblId: String!, $page: Pagination) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    associatedDiseases(page: $page) {
      count
      rows {
        score
        disease {
          id
          name
          therapeuticAreas { name }
        }
      }
    }
    tractability { modality label value }
  }
}
"""

CATEGORY_PATTERNS = {
    "genome_instability": [
        r"dna repair",
        r"double.strand",
        r"genome instability",
        r"genomic instability",
        r"damage response",
        r"homologous recombination",
        r"non.homologous",
        r"telomere",
    ],
    "repeat_control": [
        r"transpos",
        r"retrotranspos",
        r"retroelement",
        r"repeat element",
        r"piwi",
        r"pirna",
        r"line-1",
        r"line1",
        r"endogenous retrovirus",
        r"heterochromatin",
        r"setdb1",
        r"trim28",
    ],
    "senescence": [r"senescence", r"senescent", r"cell cycle arrest", r"sasp"],
}

DISEASE_AREA_KEYWORDS = {
    "cancer_disease": ["cancer", "carcinoma", "neoplasm", "tumor", "tumour", "leukemia", "lymphoma", "melanoma"],
    "neuro_disease": ["nervous", "neuro", "parkinson", "alzheimer", "ataxia"],
    "immune_inflammation_disease": ["immune", "inflammatory", "inflammation", "autoimmune"],
}

BINARY_FLAGS = [
    "genage",
    "longevitymap",
    "cellage",
    "genome_instability",
    "repeat_control",
    "senescence",
    "cancer_disease",
    "neuro_disease",
    "immune_inflammation_disease",
    "small_molecule",
]

LEAD_GENES = {"TRIM28", "SETDB1", "MOV10"}
TRANSPOSON_GENES = {"TRIM28", "SETDB1", "MOV10", "PIWIL1", "PIWIL2"}


def download(url: str, path: pathlib.Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    path.write_bytes(response.content)


def read_hagr_symbols() -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for name, path in HAGR_ZIPS.items():
        symbols: set[str] = set()
        if not path.exists():
            out[name] = symbols
            continue
        with zipfile.ZipFile(path) as zf:
            for member in zf.namelist():
                if member.endswith("/") or not member.lower().endswith((".csv", ".tsv", ".txt")):
                    continue
                text = zf.read(member).decode("utf-8", errors="replace")
                sep = "\t" if text[:4096].count("\t") >= text[:4096].count(",") else ","
                df = pd.read_csv(io.StringIO(text), sep=sep)
                lower_cols = {str(c).lower(): c for c in df.columns}
                sym_col = None
                for candidate in ["symbol", "gene symbol", "gene(s)"]:
                    if candidate in lower_cols:
                        sym_col = lower_cols[candidate]
                        break
                if sym_col is None:
                    continue
                for value in df[sym_col].dropna().astype(str):
                    for part in re.split(r"[,;/\s]+", value.upper()):
                        part = part.strip()
                        if part and part != "NAN":
                            symbols.add(part)
        out[name] = symbols
    return out


def flatten(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, (int, float)):
        return str(value)
    if isinstance(value, list):
        return " | ".join(flatten(v) for v in value)
    if isinstance(value, dict):
        return " | ".join(flatten(v) for v in value.values())
    return str(value)


def query_mygene(symbol: str) -> dict[str, Any]:
    MYGENE_CACHE.mkdir(parents=True, exist_ok=True)
    cache = MYGENE_CACHE / f"{symbol}.json"
    if cache.exists() and cache.stat().st_size > 0:
        return json.loads(cache.read_text(encoding="utf-8"))
    response = requests.get(
        MYGENE_URL,
        params={
            "q": f"symbol:{symbol}",
            "species": "human",
            "fields": "symbol,name,summary,entrezgene,ensembl.gene,go.BP,go.CC,go.MF,pathway,phenotype",
            "size": 5,
        },
        timeout=60,
    )
    response.raise_for_status()
    hits = response.json().get("hits", [])
    exact = [hit for hit in hits if str(hit.get("symbol", "")).upper() == symbol.upper()]
    hit = exact[0] if exact else hits[0] if hits else {}
    cache.write_text(json.dumps(hit, indent=2), encoding="utf-8")
    time.sleep(0.03)
    return hit


def ensembl_ids_from_hit(hit: dict[str, Any]) -> list[str]:
    value = hit.get("ensembl", {})
    if isinstance(value, list):
        ids = [str(v.get("gene", "")) for v in value if isinstance(v, dict)]
    elif isinstance(value, dict):
        ids = [str(value.get("gene", ""))]
    else:
        ids = []
    return sorted({x for x in ids if x.startswith("ENSG")})


def query_open_targets(symbol: str, ensembl_ids: list[str], top_n: int) -> tuple[dict[str, Any], str]:
    OPENTARGETS_CACHE.mkdir(parents=True, exist_ok=True)
    for ensembl_id in ensembl_ids:
        cache = OPENTARGETS_CACHE / f"{symbol}__{ensembl_id}.json"
        if cache.exists() and cache.stat().st_size > 0:
            data = json.loads(cache.read_text(encoding="utf-8"))
        else:
            response = requests.post(
                OPENTARGETS_GRAPHQL,
                json={
                    "query": TARGET_QUERY,
                    "variables": {"ensemblId": ensembl_id, "page": {"index": 0, "size": top_n}},
                },
                timeout=90,
            )
            response.raise_for_status()
            data = response.json()
            cache.write_text(json.dumps(data, indent=2), encoding="utf-8")
            time.sleep(0.08)
        target = data.get("data", {}).get("target")
        if target:
            return target, ensembl_id
    return {}, ""


def has_keyword(text: str, patterns: list[str]) -> int:
    text_l = text.lower()
    return int(any(re.search(pattern, text_l) for pattern in patterns))


def small_molecule_count(tractability: list[dict[str, Any]]) -> int:
    return sum(
        1
        for item in tractability or []
        if item.get("value") is True and str(item.get("modality", "")) == "SM"
    )


def disease_flags(rows: list[dict[str, Any]]) -> dict[str, int]:
    text = " ".join(
        str((row.get("disease") or {}).get("name", ""))
        + " "
        + " ".join(str(area.get("name", "")) for area in ((row.get("disease") or {}).get("therapeuticAreas") or []))
        for row in rows or []
    ).lower()
    return {
        flag: int(any(keyword in text for keyword in keywords))
        for flag, keywords in DISEASE_AREA_KEYWORDS.items()
    }


def build_background(args: argparse.Namespace, maintenance_symbols: set[str]) -> pd.DataFrame:
    if args.background.exists() and not args.refresh_background:
        bg = pd.read_csv(args.background, sep="\t")
        if len(bg) >= args.background_size:
            return bg.head(args.background_size).copy()
    download(HGNC_URL, args.hgnc)
    hgnc = pd.read_csv(args.hgnc, sep="\t", dtype=str)
    protein = hgnc[hgnc["locus_group"].fillna("").str.lower().eq("protein-coding gene")].copy()
    protein = protein[~protein["symbol"].isin(maintenance_symbols)]
    rng = np.random.default_rng(args.seed)
    selected = protein.sample(n=min(args.background_size, len(protein)), random_state=args.seed)
    hagr = read_hagr_symbols()
    rows = []
    for idx, symbol in enumerate(selected["symbol"].astype(str).str.upper().tolist(), start=1):
        hit = query_mygene(symbol)
        text = flatten(hit)
        target, ensembl_id = query_open_targets(symbol, ensembl_ids_from_hit(hit), args.open_targets_top_n)
        assoc_rows = ((target.get("associatedDiseases") or {}).get("rows") or []) if target else []
        max_score = max([float(row.get("score") or 0) for row in assoc_rows], default=0.0)
        flags = disease_flags(assoc_rows)
        row = {
            "human_gene_symbol": symbol,
            "background_source": "HGNC_protein_coding_random",
            "ensembl_gene_id": ensembl_id,
            "open_targets_status": "ok" if target else "no_target",
            "genage": int(symbol in hagr["genage"]),
            "longevitymap": int(symbol in hagr["longevitymap"]),
            "cellage": int(symbol in hagr["cellage"]),
            "genome_instability": has_keyword(text, CATEGORY_PATTERNS["genome_instability"]),
            "repeat_control": has_keyword(text, CATEGORY_PATTERNS["repeat_control"]),
            "senescence": has_keyword(text, CATEGORY_PATTERNS["senescence"]),
            "small_molecule": int(small_molecule_count(target.get("tractability") or []) > 0) if target else 0,
            "max_open_targets_score": max_score,
            **flags,
        }
        rows.append(row)
        if idx % 50 == 0:
            print(f"Annotated {idx}/{len(selected)} background genes")
    bg = pd.DataFrame(rows)
    args.background.parent.mkdir(parents=True, exist_ok=True)
    bg.to_csv(args.background, sep="\t", index=False)
    return bg


def maintenance_evidence(human: pd.DataFrame, translation: pd.DataFrame) -> pd.DataFrame:
    merged = human.merge(
        translation[
            [
                "human_gene_symbol",
                "max_open_targets_score",
                "small_molecule_tractability_count",
                "cancer_top_disease_count",
                "nervous_system_top_disease_count",
                "immune_inflammation_top_disease_count",
            ]
        ],
        on="human_gene_symbol",
        how="left",
    )
    out = pd.DataFrame(
        {
            "human_gene_symbol": merged["human_gene_symbol"],
            "maintenance_module": merged["maintenance_module"],
            "genage": (merged["genage_human_evidence"] == "yes").astype(int),
            "longevitymap": (merged["longevitymap_evidence"] == "yes").astype(int),
            "cellage": (merged["cellage_evidence"] == "yes").astype(int),
            "genome_instability": (merged["genome_instability_evidence"] == "keyword_supported").astype(int),
            "repeat_control": (merged["transposon_or_repeat_evidence"] == "keyword_supported").astype(int),
            "senescence": (merged["senescence_evidence"] == "keyword_supported").astype(int),
            "cancer_disease": (pd.to_numeric(merged["cancer_top_disease_count"], errors="coerce").fillna(0) > 0).astype(int),
            "neuro_disease": (pd.to_numeric(merged["nervous_system_top_disease_count"], errors="coerce").fillna(0) > 0).astype(int),
            "immune_inflammation_disease": (pd.to_numeric(merged["immune_inflammation_top_disease_count"], errors="coerce").fillna(0) > 0).astype(int),
            "small_molecule": (pd.to_numeric(merged["small_molecule_tractability_count"], errors="coerce").fillna(0) > 0).astype(int),
            "max_open_targets_score": pd.to_numeric(merged["max_open_targets_score"], errors="coerce").fillna(0),
        }
    )
    return out


def bh_adjust(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        value = min(prev, p[idx] * n / (i + 1))
        ranked[idx] = value
        prev = value
    return ranked.tolist()


def fisher_row(group_a: pd.DataFrame, group_b: pd.DataFrame, flag: str, label_a: str, label_b: str) -> dict[str, Any]:
    a_hit = int(group_a[flag].sum())
    a_no = int(len(group_a) - a_hit)
    b_hit = int(group_b[flag].sum())
    b_no = int(len(group_b) - b_hit)
    odds, p = fisher_exact([[a_hit, a_no], [b_hit, b_no]], alternative="greater")
    return {
        "test": "fisher_exact_greater",
        "comparison_a": label_a,
        "comparison_b": label_b,
        "feature": flag,
        "a_hits": a_hit,
        "a_n": len(group_a),
        "b_hits": b_hit,
        "b_n": len(group_b),
        "odds_ratio": odds,
        "p_value": p,
    }


def composite_score(df: pd.DataFrame) -> pd.Series:
    binary_sum = df[BINARY_FLAGS].sum(axis=1)
    return binary_sum + pd.to_numeric(df["max_open_targets_score"], errors="coerce").fillna(0)


def permutation_test(df: pd.DataFrame, observed_genes: set[str], sample_size: int, n_perm: int, seed: int) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    scores = composite_score(df)
    symbols = df["human_gene_symbol"].astype(str).tolist()
    score_by_symbol = dict(zip(symbols, scores))
    observed = sum(score_by_symbol[g] for g in observed_genes if g in score_by_symbol)
    all_scores = scores.to_numpy()
    if len(all_scores) < sample_size:
        return {"observed_score": observed, "p_value": math.nan, "null_mean": math.nan, "null_sd": math.nan}
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = rng.choice(all_scores, size=sample_size, replace=False).sum()
    p = (1 + np.sum(null >= observed)) / (n_perm + 1)
    return {
        "observed_score": observed,
        "p_value": p,
        "null_mean": float(null.mean()),
        "null_sd": float(null.std(ddof=1)),
        "null_q95": float(np.quantile(null, 0.95)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--human", type=pathlib.Path, default=DEFAULT_HUMAN)
    parser.add_argument("--translation", type=pathlib.Path, default=DEFAULT_TRANSLATION)
    parser.add_argument("--hgnc", type=pathlib.Path, default=DEFAULT_HGNC)
    parser.add_argument("--background", type=pathlib.Path, default=DEFAULT_BACKGROUND)
    parser.add_argument("--module-output", type=pathlib.Path, default=DEFAULT_MODULE_OUT)
    parser.add_argument("--background-output", type=pathlib.Path, default=DEFAULT_BACKGROUND_OUT)
    parser.add_argument("--permutation-output", type=pathlib.Path, default=DEFAULT_PERM_OUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    parser.add_argument("--background-size", type=int, default=300)
    parser.add_argument("--open-targets-top-n", type=int, default=10)
    parser.add_argument("--n-permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260611)
    parser.add_argument("--refresh-background", action="store_true")
    args = parser.parse_args()

    human = pd.read_csv(args.human, sep="\t")
    translation = pd.read_csv(args.translation, sep="\t")
    maint = maintenance_evidence(human, translation)
    background = build_background(args, set(maint["human_gene_symbol"].astype(str)))

    module_rows = []
    for module in sorted(maint["maintenance_module"].unique()):
        group = maint[maint["maintenance_module"] == module]
        other = maint[maint["maintenance_module"] != module]
        for flag in BINARY_FLAGS:
            module_rows.append(fisher_row(group, other, flag, module, "other_maintenance_genes"))
    module_df = pd.DataFrame(module_rows)
    module_df["p_bh"] = bh_adjust(module_df["p_value"].tolist())
    args.module_output.parent.mkdir(parents=True, exist_ok=True)
    module_df.to_csv(args.module_output, sep="\t", index=False)

    bg_rows = []
    comparisons = [
        ("maintenance_41", maint, "background_HGNC_sample", background),
        ("transposon_5", maint[maint["human_gene_symbol"].isin(TRANSPOSON_GENES)], "background_HGNC_sample", background),
    ]
    for label_a, group_a, label_b, group_b in comparisons:
        for flag in BINARY_FLAGS:
            bg_rows.append(fisher_row(group_a, group_b, flag, label_a, label_b))
        u = mannwhitneyu(
            group_a["max_open_targets_score"],
            group_b["max_open_targets_score"],
            alternative="greater",
        )
        bg_rows.append(
            {
                "test": "mannwhitneyu_greater",
                "comparison_a": label_a,
                "comparison_b": label_b,
                "feature": "max_open_targets_score",
                "a_hits": "",
                "a_n": len(group_a),
                "b_hits": "",
                "b_n": len(group_b),
                "odds_ratio": "",
                "p_value": float(u.pvalue),
                "a_mean": float(group_a["max_open_targets_score"].mean()),
                "b_mean": float(group_b["max_open_targets_score"].mean()),
            }
        )
    bg_df = pd.DataFrame(bg_rows)
    bg_df["p_bh"] = bh_adjust(bg_df["p_value"].tolist())
    args.background_output.parent.mkdir(parents=True, exist_ok=True)
    bg_df.to_csv(args.background_output, sep="\t", index=False)

    perm_rows = []
    for universe_name, universe in [
        ("maintenance_41", maint),
        ("maintenance_plus_background", pd.concat([maint, background], ignore_index=True, sort=False)),
    ]:
        lead = permutation_test(universe, LEAD_GENES, 3, args.n_permutations, args.seed)
        perm_rows.append({"test": "lead_TRIM28_SETDB1_MOV10_vs_random3", "universe": universe_name, "sample_size": 3, **lead})
        focal = permutation_test(universe, TRANSPOSON_GENES, 5, args.n_permutations, args.seed + 1)
        perm_rows.append({"test": "transposon_5_vs_random5", "universe": universe_name, "sample_size": 5, **focal})
    perm_df = pd.DataFrame(perm_rows)
    perm_df["p_bh"] = bh_adjust(perm_df["p_value"].tolist())
    args.permutation_output.parent.mkdir(parents=True, exist_ok=True)
    perm_df.to_csv(args.permutation_output, sep="\t", index=False)

    sig_module = module_df[module_df["p_bh"] < 0.05].sort_values("p_bh")
    sig_bg = bg_df[bg_df["p_bh"] < 0.05].sort_values("p_bh")
    sig_perm = perm_df[perm_df["p_bh"] < 0.05].sort_values("p_bh")
    report_lines = [
        "# Week 5 Statistical Support Report",
        "",
        "Added statistical support for the Week 5 human translation layer.",
        "",
        "## Inputs",
        "",
        f"- Maintenance genes: {len(maint)}",
        f"- HGNC protein-coding background sample: {len(background)}",
        f"- Background source: HGNC complete set, deterministic random sample, seed={args.seed}",
        f"- Permutations per test: {args.n_permutations}",
        "",
        "## Outputs",
        "",
        f"- Module enrichment: `{args.module_output}`",
        f"- Background comparison: `{args.background_output}`",
        f"- Permutation tests: `{args.permutation_output}`",
        f"- Background evidence: `{args.background}`",
        "",
        "## Significant Module Enrichments (BH < 0.05)",
        "",
    ]
    if sig_module.empty:
        report_lines.append("- None")
    else:
        for _, row in sig_module.iterrows():
            report_lines.append(
                f"- {row['comparison_a']} / {row['feature']}: {row['a_hits']}/{row['a_n']} vs "
                f"{row['b_hits']}/{row['b_n']}, OR={row['odds_ratio']:.3g}, p={row['p_value']:.3g}, BH={row['p_bh']:.3g}"
            )
    report_lines.extend(["", "## Significant Background Comparisons (BH < 0.05)", ""])
    if sig_bg.empty:
        report_lines.append("- None")
    else:
        for _, row in sig_bg.iterrows():
            extra = ""
            if row["test"] == "mannwhitneyu_greater":
                extra = f", mean_a={row.get('a_mean', ''):.3g}, mean_b={row.get('b_mean', ''):.3g}"
            report_lines.append(
                f"- {row['comparison_a']} / {row['feature']}: p={row['p_value']:.3g}, BH={row['p_bh']:.3g}{extra}"
            )
    report_lines.extend(["", "## Permutation Tests", ""])
    for _, row in perm_df.iterrows():
        report_lines.append(
            f"- {row['test']} in {row['universe']}: observed={row['observed_score']:.3f}, "
            f"null_mean={row['null_mean']:.3f}, null_q95={row['null_q95']:.3f}, "
            f"p={row['p_value']:.4g}, BH={row['p_bh']:.4g}"
        )
    report_lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "These tests ask whether the Week 5 human translation layer is stronger than expected under module-level, maintenance-gene, and random human protein-coding backgrounds. They remain exploratory because the background sample is finite and Open Targets scores are disease-evidence scores rather than ageing-causality measures.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.module_output}, {args.background_output}, {args.permutation_output}, and {args.report}")


if __name__ == "__main__":
    main()
