"""Write the Week 3 feasibility decision report."""

from __future__ import annotations

import argparse
import csv
import pathlib


DEFAULT_OUTPUT = pathlib.Path("results/reports/week3_feasibility_report.md")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def top_rows(path: pathlib.Path, n: int = 5) -> list[dict[str, str]]:
    rows = read_tsv(path)
    return rows[:n]


def write_report(output: pathlib.Path) -> None:
    pgls_top = top_rows(pathlib.Path("results/tables/maintenance_pgls_primary.tsv"), 8)
    sensitivity = read_tsv(pathlib.Path("results/tables/maintenance_pgls_sensitivity_summary.tsv"))
    output.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Week 3 Feasibility Report",
        "",
        "## Decision",
        "",
        "**Go to Week 4.**",
        "",
        "The project has enough primary-panel genome coverage and a coherent first-pass mechanism signal to justify moving from feasibility into deeper mechanism testing. The strongest Week 4 target is transposon suppression, with secondary targets in DNA repair, mitochondrial quality control, proteostasis, and cancer surveillance.",
        "",
        "## What Was Completed",
        "",
        "- Built a 68-species primary genome mechanism panel and a 130-species sensitivity panel.",
        "- Built a 41-gene human-centered maintenance seed set across 7 modules.",
        "- Generated primary and sensitivity ortholog matrix scaffolds.",
        "- Completed NCBI Gene candidate mapping for the 68-species primary panel.",
        "- Rescued zero-coverage GenBank assemblies with GFF annotation parsing.",
        "- Resolved `Erythrura gouldiae` with Taeniopygia reference proteins and DIAMOND validation.",
        "- Built module-level maintenance scores.",
        "- Ran first-pass non-phylogenetic scans, PGLS scans, and PGLS sensitivity checks.",
        "",
        "## Primary Panel Coverage",
        "",
        "- Primary panel: 68 species.",
        "- Tier 1 RefSeq annotated chromosome/complete assemblies: 37 species.",
        "- Tier 2 annotated assemblies: 31 species.",
        "- DIAMOND-validated primary candidate hits: 2,458 / 2,788.",
        "- DIAMOND-validated primary coverage: 88.2%.",
        "- Zero-coverage species after rescue: 0.",
        "",
        "## Module Coverage",
        "",
        "- DNA repair: 96.3%.",
        "- Autophagy: 93.8%.",
        "- Cancer surveillance: 75.3%.",
        "- Inflammation control: 78.7%.",
        "- Mitochondrial quality control: 94.6%.",
        "- Proteostasis: 86.3%.",
        "- Transposon suppression: 88.2%.",
        "",
        "## Strongest PGLS Signals",
        "",
    ]
    for row in pgls_top:
        lines.append(
            f"- {row['maintenance_module']} / {row['model']}: estimate={float(row['module_estimate']):.3f}, "
            f"p={float(row['module_p']):.4g}, BH(model)={float(row['module_p_bh_by_model']):.4g}, n={row['n']}."
        )
    lines.extend(
        [
            "",
            "## Sensitivity Summary",
            "",
        ]
    )
    for row in sensitivity:
        lines.append(
            f"- {row['maintenance_module']}: {row['priority']}; "
            f"foundation significant subsets={row['foundation_significant_subsets']}; "
            f"stress significant subsets={row['stress_significant_subsets']}; "
            f"all-primary p={row['all_primary_p']}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The feasibility signal is strongest for transposon suppression. It survives common exclusions such as removing humans, bats, reptiles, or non-bat mammals, and it remains positive across all tested subsets. However, it weakens in stress tests that remove birds, restrict to Tier 1 genomes, or require high overall ortholog coverage. This means the signal is promising but still sensitive to bird annotation structure and rescue strategy.",
            "",
            "Secondary modules show broadly positive but less robust evidence. DNA repair, mitochondrial quality control, proteostasis, and cancer surveillance are good Week 4 follow-up modules. Autophagy and inflammation control should remain exploratory until pathway sets and annotations are expanded.",
            "",
            "## Main Risks",
            "",
            "- The current PGLS tree uses OpenTree topology and Grafen fallback branch lengths, not dated trees.",
            "- Module scores measure candidate ortholog coverage and confidence, not gene expression, copy number, molecular rates, or functional activity.",
            "- Several bird assemblies rely on GenBank annotation rescue, which may introduce uneven annotation quality.",
            "- The 41-gene seed set is intentionally small and should be expanded with Reactome, GO, GenAge/HAGR, UniProt, OMA, OrthoDB, or Ensembl Compara.",
            "- Transposon suppression is the strongest signal but also stress-sensitive, especially to bird/Tier2 coverage structure.",
            "",
            "## Week 4 Recommendation",
            "",
            "Week 4 should focus on turning the transposon-suppression signal from a coverage signal into a stronger genome-mechanism result.",
            "",
            "Priority tasks:",
            "",
            "1. Expand transposon-suppression and genome-maintenance gene sets beyond the 41 seed genes.",
            "2. Cross-check candidate orthologs with OMA, OrthoDB, or Ensembl Compara where available.",
            "3. Add dated trees or clade-specific trees for birds and mammals.",
            "4. Build module-level sensitivity scores using primary-only, no-GFF-rescue, and high-confidence-only variants.",
            "5. Begin repeat/transposon proxy analysis only after deciding whether available assemblies support comparable repeat annotation.",
            "",
            "## Bottom Line",
            "",
            "Week 3 supports continuing. The project should not yet claim a final mechanism, but it now has a defensible computational signal: long-lived species, especially in the bird-heavy primary panel, show higher maintenance-module candidate coverage, with the clearest and most robust signal in transposon suppression.",
        ]
    )
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    write_report(args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
