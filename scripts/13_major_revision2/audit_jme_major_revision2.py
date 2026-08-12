#!/usr/bin/env python
"""Audit scientific and structural consistency of the JME revision package."""

from __future__ import annotations

import argparse
import csv
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def docx_text(path: Path) -> str:
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    return "\n".join(node.text or "" for node in root.iter(f"{{{W}}}t"))


def has_line_numbers(path: Path) -> bool:
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    return root.find(f".//{{{W}}}lnNumType") is not None


def require(condition: bool, message: str, failures: list[str]) -> None:
    if condition:
        print(f"PASS\t{message}")
    else:
        print(f"FAIL\t{message}")
        failures.append(message)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "package",
        nargs="?",
        default="submission/journal_of_molecular_evolution_major_revision4",
    )
    args = parser.parse_args()
    package = Path(args.package)
    failures: list[str] = []

    required = [
        "01_Manuscript_with_Figures.docx",
        "04_Supplementary_Information.docx",
        "07_Point_by_Point_Response.pdf",
        "Reproducibility_Guide.md",
        "file_manifest.tsv",
    ]
    for relative in required:
        path = package / relative
        require(path.is_file() and path.stat().st_size > 0, f"required file: {relative}", failures)

    for path in sorted(package.glob("*.docx")):
        try:
            with zipfile.ZipFile(path) as archive:
                bad = archive.testzip()
            require(bad is None, f"valid DOCX archive: {path.name}", failures)
        except zipfile.BadZipFile:
            require(False, f"valid DOCX archive: {path.name}", failures)

    manuscript = package / "01_Manuscript_with_Figures.docx"
    supplement = package / "04_Supplementary_Information.docx"
    manuscript_text = docx_text(manuscript)
    supplement_text = docx_text(supplement)

    require(has_line_numbers(manuscript), "continuous line numbering is encoded", failures)
    for phrase in (
        "We test three hypotheses",
        "First, loss of observed gene rows",
        "Second, random observation loss",
        "Third, assembly-level completeness",
        "Generative AI statement",
        "OpenAI ChatGPT and Codex",
        "BUSCO",
        "phenotype study effort",
        "500 replicates",
    ):
        require(phrase.lower() in manuscript_text.lower(), f"manuscript contains: {phrase}", failures)

    for stale in ("UFBoot did not converge", "inference boundary"):
        require(stale.lower() not in manuscript_text.lower(), f"stale wording absent: {stale}", failures)
    require(
        "one common axis rather than six independent confirmations" in manuscript_text.lower(),
        "module concordance is interpreted as one common axis",
        failures,
    )

    require("Supplementary Figure 9" in supplement_text, "supplement includes Figures S1-S9", failures)
    require("Sequence classification" in supplement_text, "sequence follow-up uses direct terminology", failures)

    numbered = sorted((package / "Supplementary_Data").glob("Supplementary_Data_*_*.tsv"))
    numbers = sorted(int(re.search(r"Supplementary_Data_(\d+)_", p.name).group(1)) for p in numbered)
    require(numbers == list(range(1, 47)), "Supplementary Data are numbered 1-46", failures)
    require(all(p.stat().st_size > 0 for p in numbered), "all Supplementary Data files are nonempty", failures)

    main_png = sorted((package / "Figures/Main").glob("Figure_*.png"))
    main_pdf = sorted((package / "Figures/Main").glob("Figure_*.pdf"))
    require(len(main_png) == 4 and len(main_pdf) == 4, "four main figures in PNG and PDF", failures)

    model_path = package / "Supplementary_Data/Supplementary_Data_43_effort_BUSCO_PGLS.tsv"
    table1_path = package / "Source_Data/Figure_2_module_effects.tsv"
    if model_path.exists() and table1_path.exists():
        models = read_tsv(model_path)
        base = {
            row["predictor"]: float(row["estimate_per_predictor_sd"])
            for row in models
            if row["model"] == "base" and row["predictor"] != "shared_module_pc1"
        }
        source = read_tsv(table1_path)
        source_values = {row["predictor"]: float(row["estimate_per_predictor_sd"]) for row in source}
        same = base.keys() == source_values.keys() and all(
            abs(base[key] - source_values[key]) < 1e-10 for key in base
        )
        require(same, "Figure 2 source data match base PGLS estimates", failures)
    else:
        require(False, "Figure 2 source/model tables present", failures)

    analysis = package / "Supplementary_Data/Supplementary_Data_42_study_effort_and_BUSCO_covariates.tsv"
    if analysis.exists():
        rows = read_tsv(analysis)
        require(len(rows) == 68, "revision analysis table contains 68 species", failures)
        require("busco_lineage" in rows[0], "BUSCO lineage is retained", failures)
    else:
        require(False, "revision analysis table present", failures)

    manifest = read_tsv(package / "file_manifest.tsv")
    first_column = next(iter(manifest[0]))
    manifest_paths = {row[first_column] for row in manifest}
    require("07_Point_by_Point_Response.pdf" in manifest_paths, "response PDF is in manifest", failures)

    print(f"\nAudit failures: {len(failures)}")
    if failures:
        for item in failures:
            print(f"- {item}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
