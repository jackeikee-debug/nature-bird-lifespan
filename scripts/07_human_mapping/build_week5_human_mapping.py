"""Build Week 5 human ageing and disease relevance mapping."""

from __future__ import annotations

import argparse
import csv
import io
import json
import pathlib
import re
import time
import zipfile
from collections import Counter
from typing import Any

import pandas as pd
import requests


DEFAULT_GENE_SET = pathlib.Path("data/processed/maintenance_gene_sets.tsv")
DEFAULT_RAW_DIR = pathlib.Path("data/raw/hagr")
DEFAULT_INTERIM_DIR = pathlib.Path("data/interim/human_mapping")
DEFAULT_OUTPUT = pathlib.Path("data/processed/human_mapping.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/week5_human_mapping_report.md")

HAGR_URLS = {
    "genage_human": "https://genomics.senescence.info/genes/human_genes.zip",
    "longevity_map": "https://genomics.senescence.info/longevity/longevity_genes.zip",
    "cellage": "https://genomics.senescence.info/cells/cellAge.zip",
}

MYGENE_URL = "https://mygene.info/v3/query"

CATEGORY_PATTERNS = {
    "genome_instability_evidence": [
        r"dna repair",
        r"double.strand",
        r"single.strand",
        r"genome instability",
        r"genomic instability",
        r"chromosome instability",
        r"telomere",
        r"damage response",
        r"homologous recombination",
        r"non.homologous",
    ],
    "cancer_evidence": [
        r"cancer",
        r"carcinoma",
        r"tumou?r",
        r"neoplasm",
        r"oncogene",
        r"tumou?r suppressor",
        r"leukemia",
        r"lymphoma",
    ],
    "inflammation_evidence": [
        r"inflammation",
        r"inflammatory",
        r"cytokine",
        r"nf.?kappa",
        r"interleukin",
        r"tnf",
        r"inflammasome",
        r"immune",
    ],
    "transposon_or_repeat_evidence": [
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
    "proteostasis_evidence": [
        r"proteostasis",
        r"chaperone",
        r"heat shock",
        r"ubiquitin",
        r"proteasome",
        r"protein folding",
    ],
    "mitochondrial_evidence": [
        r"mitochond",
        r"mitophagy",
        r"oxidative phosphorylation",
        r"reactive oxygen",
        r"parkinson",
    ],
    "autophagy_evidence": [
        r"autophagy",
        r"autophagosome",
        r"lysosome",
        r"mitophagy",
    ],
    "senescence_evidence": [
        r"senescence",
        r"senescent",
        r"cell cycle arrest",
        r"sasp",
    ],
}

TRANSPOSON_GENES = {"PIWIL1", "PIWIL2", "MOV10", "TRIM28", "SETDB1"}


def download(url: str, path: pathlib.Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    response = requests.get(url, timeout=90)
    response.raise_for_status()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(response.content)


def read_zip_tables(path: pathlib.Path) -> dict[str, pd.DataFrame]:
    tables: dict[str, pd.DataFrame] = {}
    with zipfile.ZipFile(path) as zf:
        for name in zf.namelist():
            if name.endswith("/") or not re.search(r"\.(csv|tsv|txt)$", name, re.I):
                continue
            raw = zf.read(name)
            text = raw.decode("utf-8", errors="replace")
            sample = text[:4096]
            delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
            try:
                df = pd.read_csv(io.StringIO(text), sep=delimiter)
            except Exception:
                df = pd.read_csv(io.StringIO(text), sep=None, engine="python")
            tables[name] = df
    return tables


def find_symbol_column(df: pd.DataFrame) -> str | None:
    candidates = [
        "symbol",
        "gene symbol",
        "gene_symbol",
        "hgnc symbol",
        "hgnc_symbol",
        "Gene Symbol",
        "Symbol",
        "HGNC symbol",
        "Gene(s)",
    ]
    lower_map = {str(col).strip().lower(): col for col in df.columns}
    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    for col in df.columns:
        col_l = str(col).lower()
        if "symbol" in col_l or col_l in {"gene", "name", "gene(s)"}:
            return col
    return None


def collect_symbols(tables: dict[str, pd.DataFrame]) -> tuple[set[str], dict[str, str]]:
    symbols: set[str] = set()
    evidence_text: dict[str, str] = {}
    for _, df in tables.items():
        sym_col = find_symbol_column(df)
        if sym_col is None:
            continue
        for _, row in df.iterrows():
            raw_symbol = str(row.get(sym_col, "")).strip().upper()
            if not raw_symbol or raw_symbol == "NAN":
                continue
            row_symbols = [
                part.strip()
                for part in re.split(r"[,;/\s]+", raw_symbol)
                if part.strip() and part.strip() != "NAN"
            ]
            row_text = " | ".join(
                str(row.get(col, ""))
                for col in df.columns
                if pd.notna(row.get(col, ""))
            )
            for symbol in row_symbols:
                symbols.add(symbol)
                if symbol not in evidence_text or len(row_text) > len(evidence_text[symbol]):
                    evidence_text[symbol] = row_text
    return symbols, evidence_text


def query_mygene(symbols: list[str]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    fields = "symbol,name,summary,entrezgene,ensembl.gene,uniprot.Swiss-Prot,go.BP,go.CC,go.MF,pathway,phenotype,refseq"
    for symbol in symbols:
        response = requests.get(
            MYGENE_URL,
            params={
                "q": f"symbol:{symbol}",
                "species": "human",
                "fields": fields,
                "size": 5,
            },
            timeout=60,
        )
        response.raise_for_status()
        hits = response.json().get("hits", [])
        exact = [
            hit for hit in hits
            if str(hit.get("symbol", "")).upper() == symbol
        ]
        if exact:
            out[symbol] = exact[0]
        elif hits:
            out[symbol] = hits[0]
        else:
            out[symbol] = {}
        time.sleep(0.05)
    return out


def flatten_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, (int, float)):
        return str(value)
    if isinstance(value, list):
        return " | ".join(flatten_text(v) for v in value)
    if isinstance(value, dict):
        return " | ".join(flatten_text(v) for v in value.values())
    return str(value)


def evidence_flag(text: str, patterns: list[str]) -> str:
    text_l = text.lower()
    matches = [pat for pat in patterns if re.search(pat, text_l)]
    if matches:
        return "keyword_supported"
    return "not_detected"


def priority(row: dict[str, str]) -> str:
    gene = row["human_gene_symbol"]
    flags = [
        row["genage_human_evidence"] == "yes",
        row["cellage_evidence"] == "yes",
        row["longevitymap_evidence"] == "yes",
        row["transposon_or_repeat_evidence"] == "keyword_supported",
        row["genome_instability_evidence"] == "keyword_supported",
        row["senescence_evidence"] == "keyword_supported",
    ]
    if gene in TRANSPOSON_GENES and sum(flags) >= 2:
        return "high_transposon_translation_priority"
    if row["maintenance_module"] == "transposon_suppression":
        return "transposon_module_follow_up"
    if sum(flags) >= 3:
        return "high_general_maintenance_priority"
    if sum(flags) >= 1:
        return "supporting_maintenance_context"
    return "low_current_human_mapping_support"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-set", type=pathlib.Path, default=DEFAULT_GENE_SET)
    parser.add_argument("--raw-dir", type=pathlib.Path, default=DEFAULT_RAW_DIR)
    parser.add_argument("--interim-dir", type=pathlib.Path, default=DEFAULT_INTERIM_DIR)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    args.raw_dir.mkdir(parents=True, exist_ok=True)
    args.interim_dir.mkdir(parents=True, exist_ok=True)

    raw_paths = {}
    for name, url in HAGR_URLS.items():
        path = args.raw_dir / f"{name}.zip"
        download(url, path)
        raw_paths[name] = path

    source_symbols: dict[str, set[str]] = {}
    source_text: dict[str, dict[str, str]] = {}
    source_shapes = []
    for name, path in raw_paths.items():
        tables = read_zip_tables(path)
        symbols, text = collect_symbols(tables)
        source_symbols[name] = symbols
        source_text[name] = text
        for table_name, df in tables.items():
            df.to_csv(args.interim_dir / f"{name}__{pathlib.Path(table_name).name}.tsv", sep="\t", index=False)
            source_shapes.append((name, table_name, len(df), len(df.columns), find_symbol_column(df) or ""))

    genes = pd.read_csv(args.gene_set, sep="\t")
    symbols = genes["human_gene_symbol"].dropna().astype(str).str.upper().tolist()
    mygene = query_mygene(symbols)
    time.sleep(0.2)

    rows: list[dict[str, str]] = []
    for _, gene_row in genes.iterrows():
        symbol = str(gene_row["human_gene_symbol"]).upper()
        hit = mygene.get(symbol, {})
        text_parts = [
            flatten_text(hit.get("name")),
            flatten_text(hit.get("summary")),
            flatten_text(hit.get("go")),
            flatten_text(hit.get("pathway")),
            flatten_text(hit.get("phenotype")),
            source_text.get("genage_human", {}).get(symbol, ""),
            source_text.get("longevity_map", {}).get(symbol, ""),
            source_text.get("cellage", {}).get(symbol, ""),
        ]
        evidence_text_blob = " | ".join(x for x in text_parts if x)
        row = {
            "human_gene_symbol": symbol,
            "maintenance_module": str(gene_row["maintenance_module"]),
            "entrezgene": str(hit.get("entrezgene", "")),
            "human_gene_name": str(hit.get("name", "")),
            "mygene_symbol": str(hit.get("symbol", "")),
            "genage_human_evidence": "yes" if symbol in source_symbols.get("genage_human", set()) else "no",
            "longevitymap_evidence": "yes" if symbol in source_symbols.get("longevity_map", set()) else "no",
            "cellage_evidence": "yes" if symbol in source_symbols.get("cellage", set()) else "no",
            "human_aging_evidence": "",
            "disease_association_source": "HAGR_GenAge;HAGR_LongevityMap;HAGR_CellAge;MyGene_summary_keywords",
            "drug_or_targetability_note": "not_assessed_week5_v1",
            "evidence_text_short": evidence_text_blob[:900].replace("\n", " "),
        }
        for category, patterns in CATEGORY_PATTERNS.items():
            row[category] = evidence_flag(evidence_text_blob, patterns)
        aging_sources = []
        if row["genage_human_evidence"] == "yes":
            aging_sources.append("GenAge_human")
        if row["longevitymap_evidence"] == "yes":
            aging_sources.append("LongevityMap")
        if row["cellage_evidence"] == "yes":
            aging_sources.append("CellAge")
        if row["senescence_evidence"] == "keyword_supported":
            aging_sources.append("senescence_keyword")
        row["human_aging_evidence"] = ";".join(aging_sources) if aging_sources else "not_detected_in_week5_v1_sources"
        row["interpretation_priority"] = priority(row)
        rows.append(row)

    output = pd.DataFrame(rows)
    output = output[
        [
            "human_gene_symbol",
            "maintenance_module",
            "entrezgene",
            "human_gene_name",
            "mygene_symbol",
            "human_aging_evidence",
            "genage_human_evidence",
            "longevitymap_evidence",
            "cellage_evidence",
            "genome_instability_evidence",
            "cancer_evidence",
            "inflammation_evidence",
            "transposon_or_repeat_evidence",
            "proteostasis_evidence",
            "mitochondrial_evidence",
            "autophagy_evidence",
            "senescence_evidence",
            "disease_association_source",
            "drug_or_targetability_note",
            "interpretation_priority",
            "evidence_text_short",
        ]
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)

    counts = Counter(output["interpretation_priority"])
    module_counts = output.groupby("maintenance_module")[
        ["genage_human_evidence", "longevitymap_evidence", "cellage_evidence"]
    ].apply(lambda df: pd.Series({col: int((df[col] == "yes").sum()) for col in df.columns}))
    source_shape_lines = [
        f"- {name} / {table_name}: rows={nrow}, columns={ncol}, symbol_column={sym_col}"
        for name, table_name, nrow, ncol, sym_col in source_shapes
    ]
    trans_rows = output[output["maintenance_module"] == "transposon_suppression"]
    report_lines = [
        "# Week 5 Human Mapping Report",
        "",
        "Built the first reproducible human ageing/disease relevance map for the 41 maintenance seed genes.",
        "",
        "## Sources",
        "",
        "- HAGR GenAge human genes: official zipped tab-delimited dataset.",
        "- HAGR LongevityMap: official zipped longevity gene dataset.",
        "- HAGR CellAge: official zipped cellular senescence gene dataset.",
        "- MyGene.info: human gene name, summary, GO/pathway/phenotype keyword context.",
        "",
        "Downloaded/interim tables:",
        *source_shape_lines,
        "",
        "## Output",
        "",
        f"- Rows: {len(output)}",
        f"- Genes: {output['human_gene_symbol'].nunique()}",
        f"- Output: `{args.output}`",
        "",
        "## Interpretation Priorities",
        "",
        *[f"- {key}: {value}" for key, value in sorted(counts.items())],
        "",
        "## Source Hits by Module",
        "",
        module_counts.to_markdown(),
        "",
        "## Transposon Module Snapshot",
        "",
    ]
    for _, row in trans_rows.iterrows():
        report_lines.append(
            f"- {row['human_gene_symbol']}: aging={row['human_aging_evidence']}; "
            f"repeat={row['transposon_or_repeat_evidence']}; "
            f"senescence={row['senescence_evidence']}; priority={row['interpretation_priority']}."
        )
    report_lines.extend(
        [
            "",
            "## Caveats",
            "",
            "- This is a first-pass computational map. Keyword evidence is useful for triage but should be replaced or supplemented with curated disease/drug databases before manuscript claims.",
            "- Transposon-module translation should focus on pathway plausibility, not direct claims that these exact genes explain human longevity.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    provenance = {
        "hagr_urls": HAGR_URLS,
        "mygene_url": MYGENE_URL,
        "output": str(args.output),
        "rows": len(output),
        "source_shapes": source_shapes,
    }
    (args.interim_dir / "week5_human_mapping_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    print(f"Wrote {args.output} and {args.report}")


if __name__ == "__main__":
    main()
