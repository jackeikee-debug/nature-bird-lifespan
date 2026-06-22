"""Run selected full200 expansion batches end-to-end."""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import pandas as pd


def slug_for_batch(batch_id: str) -> str:
    body = batch_id.split("__", 1)[1] if "__" in batch_id else batch_id
    body = body.replace("-", "_")
    body = body.replace("__", "_")
    return f"phase2_full200_{body}"


def run(cmd: list[str]) -> None:
    print("RUN", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=pathlib.Path, required=True)
    parser.add_argument("--batch-summary", type=pathlib.Path, required=True)
    parser.add_argument("--wave", default="W3_standard_mapping_full_species")
    parser.add_argument("--modules", nargs="*", default=[])
    parser.add_argument("--max-batches", type=int, default=0)
    parser.add_argument("--sleep", type=float, default=0.34)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--manifest-output", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()

    summary = pd.read_csv(args.batch_summary, sep="\t")
    subset = summary[summary["expansion_wave"].eq(args.wave)].copy()
    if args.modules:
        subset = subset[subset["maintenance_module_v2"].isin(args.modules)].copy()
    subset = subset.sort_values(["maintenance_module_v2", "batch_id"])
    if args.max_batches and args.max_batches > 0:
        subset = subset.head(args.max_batches)
    if subset.empty:
        raise SystemExit("No batches selected.")

    manifest_rows = []
    for _, batch in subset.iterrows():
        batch_id = batch["batch_id"]
        slug = slug_for_batch(batch_id)
        queue_out = pathlib.Path("data/processed") / f"{slug}_queue.tsv"
        queue_report = pathlib.Path("results/reports") / f"{slug}_queue_report.md"
        results_out = pathlib.Path("data/processed") / f"{slug}_ncbi_gene.tsv"
        gene_summary = pathlib.Path("results/tables") / f"{slug}_ncbi_gene_summary.tsv"
        ncbi_report = pathlib.Path("results/reports") / f"{slug}_ncbi_gene_report.md"
        qc_out = pathlib.Path("results/tables") / f"{slug}_qc.tsv"
        qc_report = pathlib.Path("results/reports") / f"{slug}_qc_report.md"
        cache_dir = pathlib.Path("data/interim/phase2") / f"ncbi_gene_{slug}_cache"

        if not queue_out.exists():
            run(
                [
                    args.python,
                    "scripts/09_phase2/select_full200_expansion_batch.py",
                    "--queue",
                    str(args.queue),
                    "--batch-id",
                    batch_id,
                    "--output",
                    str(queue_out),
                    "--report",
                    str(queue_report),
                ]
            )
        else:
            print(f"SKIP existing {queue_out}", flush=True)

        if not results_out.exists():
            run(
                [
                    args.python,
                    "scripts/09_phase2/run_p2_2_ncbi_gene_pilot.py",
                    "--input",
                    str(queue_out),
                    "--output",
                    str(results_out),
                    "--summary-output",
                    str(gene_summary),
                    "--report",
                    str(ncbi_report),
                    "--cache-dir",
                    str(cache_dir),
                    "--sleep",
                    str(args.sleep),
                ]
            )
        else:
            print(f"SKIP existing {results_out}", flush=True)

        if not qc_out.exists():
            run(
                [
                    args.python,
                    "scripts/09_phase2/summarize_full200_batch_ncbi_gene.py",
                    "--input",
                    str(results_out),
                    "--gene-summary",
                    str(gene_summary),
                    "--qc-output",
                    str(qc_out),
                    "--report",
                    str(qc_report),
                ]
            )
        else:
            print(f"SKIP existing {qc_out}", flush=True)

        manifest_rows.append(
            {
                "batch_id": batch_id,
                "slug": slug,
                "maintenance_module_v2": batch["maintenance_module_v2"],
                "genes": int(batch["genes"]),
                "species_rows_needed": int(batch["species_rows_needed"]),
                "queue": str(queue_out),
                "results": str(results_out),
                "gene_summary": str(gene_summary),
                "qc": str(qc_out),
                "report": str(qc_report),
            }
        )

    manifest = pd.DataFrame(manifest_rows)
    args.manifest_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.manifest_output, sep="\t", index=False)

    lines = [
        "# Full200 Expansion Batch Runner Report",
        "",
        f"Wave: {args.wave}",
        f"Batches selected: {len(manifest)}",
        "",
        "## Batches",
        "",
    ]
    for _, row in manifest.iterrows():
        lines.append(
            f"- {row['batch_id']}: module={row['maintenance_module_v2']}, genes={row['genes']}, rows={row['species_rows_needed']}"
        )
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.manifest_output} and {args.report}")


if __name__ == "__main__":
    main()
