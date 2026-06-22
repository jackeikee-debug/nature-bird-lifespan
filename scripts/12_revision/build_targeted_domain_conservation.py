#!/usr/bin/env python3
"""Project human Pfam boundaries onto cross-species alignments and score conservation."""

from pathlib import Path

import numpy as np
import pandas as pd


GENES = ["DNMT1", "DNMT3A", "DNMT3B", "HELLS", "MBD2", "MBD3", "MORC3", "SAMHD1", "SETDB2", "UHRF1"]
WORK = Path("data/interim/protein_conservation")
INTERPRO = Path("results/tables/human_reference_10genes_interpro_pfam.tsv")
ROWS = Path("results/tables/targeted_protein_conservation_rows.tsv")
ROW_OUTPUT = Path("results/tables/targeted_domain_conservation_rows.tsv")
INTERVAL_OUTPUT = Path("results/tables/targeted_domain_conservation_intervals.tsv")
GENE_OUTPUT = Path("results/tables/targeted_domain_conservation_by_gene.tsv")
REPORT = Path("results/reports/targeted_domain_conservation_report.md")

INTERPRO_COLUMNS = [
    "record_id", "md5", "protein_length", "analysis", "signature_accession",
    "signature_description", "domain_start", "domain_end", "score", "status",
    "run_date", "interpro_accession", "interpro_description", "go_terms", "pathways",
]


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header = None
    sequence = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(sequence)
            header = line[1:].split()[0]
            sequence = []
        else:
            sequence.append(line)
    if header is not None:
        records[header] = "".join(sequence)
    return records


def raw_position_to_alignment(reference_aligned: str) -> dict[int, int]:
    mapping = {}
    raw_position = 0
    for column, amino_acid in enumerate(reference_aligned):
        if amino_acid != "-":
            raw_position += 1
            mapping[raw_position] = column
    return mapping


def region_metrics(reference: str, target: str, columns: set[int]) -> dict[str, float]:
    reference_columns = [column for column in sorted(columns) if reference[column] != "-"]
    paired = [(reference[column], target[column]) for column in reference_columns if target[column] != "-"]
    coverage = len(paired) / len(reference_columns) if reference_columns else np.nan
    identity = sum(a == b for a, b in paired) / len(paired) if paired else np.nan
    return {
        "reference_residues": len(reference_columns),
        "paired_residues": len(paired),
        "coverage": coverage,
        "identity": identity,
        "identity_coverage_product": identity * coverage if np.isfinite(identity) and np.isfinite(coverage) else np.nan,
    }


def main() -> None:
    interpro = pd.read_csv(INTERPRO, sep="\t", names=INTERPRO_COLUMNS)
    interpro["human_gene_symbol"] = interpro["record_id"].str.replace("REF_", "", regex=False)
    interpro["domain_instance"] = interpro.groupby(["human_gene_symbol", "signature_accession"]).cumcount() + 1
    rows = pd.read_csv(ROWS, sep="\t")
    metric_rows = []
    interval_rows = []

    for gene in GENES:
        aligned = read_fasta(WORK / f"{gene}.aligned.faa")
        reference = aligned["REF_Homo_sapiens"]
        position_map = raw_position_to_alignment(reference)
        domains = interpro.loc[interpro["human_gene_symbol"].eq(gene)].copy()
        if domains.empty:
            raise RuntimeError(f"No Pfam domains found for {gene}")
        domain_columns: set[int] = set()
        interval_columns: dict[tuple[str, int], set[int]] = {}
        for _, domain in domains.iterrows():
            raw_positions = range(int(domain["domain_start"]), int(domain["domain_end"]) + 1)
            columns = {position_map[position] for position in raw_positions if position in position_map}
            key = (domain["signature_accession"], int(domain["domain_instance"]))
            interval_columns[key] = columns
            domain_columns.update(columns)
        all_reference_columns = {column for column, amino_acid in enumerate(reference) if amino_acid != "-"}
        nondomain_columns = all_reference_columns - domain_columns

        gene_rows = rows.loc[rows["human_gene_symbol"].eq(gene)].copy()
        for _, row in gene_rows.iterrows():
            record_id = str(row.get("alignment_record_id", ""))
            available = bool(row["sequence_available"]) and record_id in aligned
            base = {
                "scientific_name": row["scientific_name"],
                "opentree_tip_label": row.get("opentree_tip_label", ""),
                "clade": row["clade"],
                "flight_status": row["flight_status"],
                "human_gene_symbol": gene,
                "alignment_record_id": record_id if available else "",
                "selected_accession": row.get("selected_accession", ""),
                "sequence_source": row.get("sequence_source", ""),
                "sequence_available": available,
                "pfam_domain_instances": len(domains),
                "pfam_domain_signatures": domains["signature_accession"].nunique(),
                "human_domain_reference_residues": len(domain_columns),
                "human_nondomain_reference_residues": len(nondomain_columns),
            }
            if not available:
                base.update(
                    {
                        "domain_reference_coverage": np.nan,
                        "domain_aligned_identity": np.nan,
                        "domain_identity_coverage_product": np.nan,
                        "nondomain_reference_coverage": np.nan,
                        "nondomain_aligned_identity": np.nan,
                        "nondomain_identity_coverage_product": np.nan,
                        "domain_minus_nondomain_identity": np.nan,
                        "domain_minus_nondomain_product": np.nan,
                    }
                )
                metric_rows.append(base)
                continue

            target = aligned[record_id]
            domain_metric = region_metrics(reference, target, domain_columns)
            nondomain_metric = region_metrics(reference, target, nondomain_columns)
            base.update(
                {
                    "domain_reference_coverage": domain_metric["coverage"],
                    "domain_aligned_identity": domain_metric["identity"],
                    "domain_identity_coverage_product": domain_metric["identity_coverage_product"],
                    "nondomain_reference_coverage": nondomain_metric["coverage"],
                    "nondomain_aligned_identity": nondomain_metric["identity"],
                    "nondomain_identity_coverage_product": nondomain_metric["identity_coverage_product"],
                    "domain_minus_nondomain_identity": domain_metric["identity"] - nondomain_metric["identity"],
                    "domain_minus_nondomain_product": domain_metric["identity_coverage_product"] - nondomain_metric["identity_coverage_product"],
                }
            )
            metric_rows.append(base)

            for _, domain in domains.iterrows():
                key = (domain["signature_accession"], int(domain["domain_instance"]))
                interval_metric = region_metrics(reference, target, interval_columns[key])
                interval_rows.append(
                    {
                        "scientific_name": row["scientific_name"],
                        "clade": row["clade"],
                        "human_gene_symbol": gene,
                        "selected_accession": row.get("selected_accession", ""),
                        "signature_accession": domain["signature_accession"],
                        "signature_description": domain["signature_description"],
                        "interpro_accession": domain["interpro_accession"],
                        "interpro_description": domain["interpro_description"],
                        "domain_instance": int(domain["domain_instance"]),
                        "human_domain_start": int(domain["domain_start"]),
                        "human_domain_end": int(domain["domain_end"]),
                        "human_domain_length": interval_metric["reference_residues"],
                        "paired_residues": interval_metric["paired_residues"],
                        "domain_reference_coverage": interval_metric["coverage"],
                        "domain_aligned_identity": interval_metric["identity"],
                        "domain_identity_coverage_product": interval_metric["identity_coverage_product"],
                    }
                )

    metrics = pd.DataFrame(metric_rows).sort_values(["human_gene_symbol", "scientific_name"])
    intervals = pd.DataFrame(interval_rows).sort_values(
        ["human_gene_symbol", "signature_accession", "domain_instance", "scientific_name"]
    )
    summaries = (
        metrics.loc[metrics["sequence_available"]]
        .groupby("human_gene_symbol", as_index=False)
        .agg(
            species_with_sequence=("scientific_name", "nunique"),
            pfam_domain_instances=("pfam_domain_instances", "first"),
            human_domain_reference_residues=("human_domain_reference_residues", "first"),
            median_domain_coverage=("domain_reference_coverage", "median"),
            median_domain_identity=("domain_aligned_identity", "median"),
            median_domain_product=("domain_identity_coverage_product", "median"),
            median_nondomain_identity=("nondomain_aligned_identity", "median"),
            median_domain_minus_nondomain_identity=("domain_minus_nondomain_identity", "median"),
        )
    )
    ROW_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(ROW_OUTPUT, sep="\t", index=False, na_rep="")
    intervals.to_csv(INTERVAL_OUTPUT, sep="\t", index=False, na_rep="")
    summaries.to_csv(GENE_OUTPUT, sep="\t", index=False, na_rep="")

    report = [
        "# Targeted Pfam-Domain Conservation",
        "",
        f"Pfam identified {len(interpro)} domain instances across all {len(GENES)} reviewed human reference proteins. Domain coordinates were projected onto the existing per-gene MAFFT alignments; no target-species domain annotation was assumed.",
        "",
        f"- Gene-species rows: {len(metrics)}.",
        f"- Rows with aligned protein sequence: {int(metrics['sequence_available'].sum())}.",
        f"- Domain interval comparison rows: {len(intervals)}.",
        f"- Rows with at least 50% human-domain coverage: {int((metrics['domain_reference_coverage'] >= 0.5).sum())}.",
        "",
        "Domain conservation is a sequence-level proxy. It does not measure protein activity, expression, chromatin state, or repeat burden. Human-reference boundaries also make the metric best suited to sensitivity analysis rather than claims of lineage-specific domain gain or loss.",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Wrote {ROW_OUTPUT}, {INTERVAL_OUTPUT}, {GENE_OUTPUT}, and {REPORT}")


if __name__ == "__main__":
    main()
