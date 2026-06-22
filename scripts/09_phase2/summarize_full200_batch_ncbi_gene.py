"""Summarize a full200 expansion NCBI Gene batch."""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, required=True)
    parser.add_argument("--gene-summary", type=pathlib.Path, required=True)
    parser.add_argument("--qc-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    rows = pd.read_csv(args.input, sep="\t")
    genes = pd.read_csv(args.gene_summary, sep="\t")
    total = len(rows)
    candidate_found = int((rows["ncbi_pilot_status"] == "candidate_found").sum())
    no_candidate = int((rows["ncbi_pilot_status"] == "no_candidate").sum())
    query_error = int((rows["ncbi_pilot_status"] == "query_error").sum())
    low_conf = int(
        (rows["ncbi_pilot_status"].eq("candidate_found") & rows["ncbi_symbol_confidence"].ne("medium")).sum()
    )
    query_error_fraction = query_error / total if total else 0.0
    candidate_fraction = candidate_found / total if total else 0.0
    continue_recommended = query_error_fraction <= 0.05

    clade = (
        rows.groupby(["clade", "ncbi_pilot_status"], as_index=False)
        .size()
        .pivot(index="clade", columns="ncbi_pilot_status", values="size")
        .fillna(0)
        .astype(int)
        .reset_index()
    )

    qc = pd.DataFrame(
        [
            {
                "batch_id": rows["batch_id"].iloc[0] if "batch_id" in rows.columns and len(rows) else "",
                "rows": total,
                "genes": rows["human_gene_symbol"].nunique(),
                "candidate_found": candidate_found,
                "no_candidate": no_candidate,
                "query_error": query_error,
                "candidate_fraction": candidate_fraction,
                "query_error_fraction": query_error_fraction,
                "non_medium_candidate_rows": low_conf,
                "continue_recommended": continue_recommended,
                "next_step": "continue_W2_batches" if continue_recommended else "pause_and_fix_query_errors",
            }
        ]
    )
    args.qc_output.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(args.qc_output, sep="\t", index=False)

    lines = [
        "# Full200 Sensitivity Batch NCBI Gene QC Report",
        "",
        f"Batch: {qc['batch_id'].iloc[0]}",
        f"Rows queried: {total}",
        f"Genes queried: {rows['human_gene_symbol'].nunique()}",
        f"Candidate found: {candidate_found} ({candidate_fraction:.3f})",
        f"No candidate: {no_candidate}",
        f"Query errors: {query_error} ({query_error_fraction:.3f})",
        f"Non-medium candidate rows: {low_conf}",
        f"Continue recommended: {continue_recommended}",
        "",
        "## Gene Summary",
        "",
    ]
    for _, row in genes.sort_values("candidate_fraction").iterrows():
        lines.append(
            f"- {row['human_gene_symbol']}: {int(row['candidate_found'])}/{int(row['rows'])} "
            f"candidate_found, fraction={row['candidate_fraction']:.3f}"
        )
    lines.extend(["", "## Clade Counts", ""])
    for _, row in clade.iterrows():
        parts = [f"{col}={row[col]}" for col in clade.columns if col != "clade"]
        lines.append(f"- {row['clade']}: " + ", ".join(parts))
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This NCBI Gene batch is a first-pass symbol/taxid expansion. Candidate rows support sensitivity mapping, but no-candidate rows are not biological absences without external sequence or orthology evidence.",
        ]
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.report}")


if __name__ == "__main__":
    main()
