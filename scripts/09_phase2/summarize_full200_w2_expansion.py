"""Combine and summarize full200 NCBI Gene expansion batches."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def read_many(paths: list[pathlib.Path]) -> pd.DataFrame:
    frames = []
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path, sep="\t")
        frame["source_file"] = path.name
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-results", nargs="+", type=pathlib.Path, required=True)
    parser.add_argument("--batch-qc", nargs="+", type=pathlib.Path, required=True)
    parser.add_argument("--combined-output", type=pathlib.Path, required=True)
    parser.add_argument("--gene-summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--batch-qc-summary-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    rows = read_many(args.batch_results)
    qc = read_many(args.batch_qc)

    args.combined_output.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(args.combined_output, sep="\t", index=False)

    gene_summary = (
        rows.groupby(["human_gene_symbol", "maintenance_module_v2", "submodule_v2"], as_index=False)
        .agg(
            rows=("scientific_name", "size"),
            candidate_found=("ncbi_pilot_status", lambda x: int((x == "candidate_found").sum())),
            no_candidate=("ncbi_pilot_status", lambda x: int((x == "no_candidate").sum())),
            query_error=("ncbi_pilot_status", lambda x: int((x == "query_error").sum())),
            medium_confidence=("ncbi_symbol_confidence", lambda x: int((x == "medium").sum())),
            low_confidence=("ncbi_symbol_confidence", lambda x: int((x == "low").sum())),
        )
        .assign(candidate_fraction=lambda d: d["candidate_found"] / d["rows"])
        .sort_values(["maintenance_module_v2", "candidate_fraction", "human_gene_symbol"])
    )
    args.gene_summary_output.parent.mkdir(parents=True, exist_ok=True)
    gene_summary.to_csv(args.gene_summary_output, sep="\t", index=False)

    qc = qc.sort_values("batch_id")
    args.batch_qc_summary_output.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(args.batch_qc_summary_output, sep="\t", index=False)

    clade_summary = (
        rows.groupby(["maintenance_module_v2", "clade", "ncbi_pilot_status"], as_index=False)
        .size()
        .pivot(index=["maintenance_module_v2", "clade"], columns="ncbi_pilot_status", values="size")
        .fillna(0)
        .astype(int)
        .reset_index()
    )

    total = len(rows)
    candidate_found = int((rows["ncbi_pilot_status"] == "candidate_found").sum())
    no_candidate = int((rows["ncbi_pilot_status"] == "no_candidate").sum())
    query_error = int((rows["ncbi_pilot_status"] == "query_error").sum())
    lines = [
        "# Full200 NCBI Gene Expansion Summary",
        "",
        f"Batches combined: {len(args.batch_results)}",
        f"Genes queried: {rows['human_gene_symbol'].nunique()}",
        f"Rows queried: {total}",
        f"Candidate found: {candidate_found} ({candidate_found / total:.3f})",
        f"No candidate: {no_candidate}",
        f"Query errors: {query_error}",
        "",
        "## Batch QC",
        "",
    ]
    for _, row in qc.iterrows():
        lines.append(
            f"- {row['batch_id']}: candidate_found={int(row['candidate_found'])}/{int(row['rows'])} "
            f"({row['candidate_fraction']:.3f}), query_error={int(row['query_error'])}, "
            f"continue={row['continue_recommended']}"
        )
    lines.extend(["", "## Lowest-Coverage Genes", ""])
    for _, row in gene_summary.head(12).iterrows():
        lines.append(
            f"- {row['human_gene_symbol']} ({row['maintenance_module_v2']}): "
            f"{int(row['candidate_found'])}/{int(row['rows'])}, fraction={row['candidate_fraction']:.3f}"
        )
    lines.extend(["", "## Clade Summary", ""])
    for _, row in clade_summary.iterrows():
        parts = [f"{col}={row[col]}" for col in clade_summary.columns if col not in {"maintenance_module_v2", "clade"}]
        lines.append(f"- {row['maintenance_module_v2']} / {row['clade']}: " + ", ".join(parts))
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Full-species NCBI Gene expansion is complete at the first-pass symbol/taxid level. Candidate rows can enter sensitivity mapping; no-candidate rows remain protected from absence claims and should be triaged for protein/external orthology follow-up if they influence module scores.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.report}")


if __name__ == "__main__":
    main()
