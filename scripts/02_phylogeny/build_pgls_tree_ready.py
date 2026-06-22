"""Join PGLS phenotype rows with OpenTree TNRS matches."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import Counter


DEFAULT_PGLS = pathlib.Path("data/processed/pgls_species.tsv")
DEFAULT_TNRS = pathlib.Path("data/processed/opentree_tnrs_matches.tsv")
DEFAULT_OUTPUT = pathlib.Path("data/processed/pgls_tree_ready.tsv")
DEFAULT_REVIEW = pathlib.Path("data/processed/pgls_tree_manual_review.tsv")
DEFAULT_REPORT = pathlib.Path("results/reports/opentree_tnrs_report.md")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def tree_label(name: str) -> str:
    return name.strip().replace(" ", "_")


def review_reason(row: dict[str, str]) -> str:
    reasons = []
    if row.get("pgls_tree_status") != "ready":
        reasons.append("unmatched")
    if row.get("opentree_rank") != "species":
        reasons.append(f"rank_{row.get('opentree_rank') or 'missing'}")
    if row.get("opentree_is_synonym") == "True":
        reasons.append("opentree_synonym")
    if row.get("opentree_is_approximate") == "True":
        reasons.append("opentree_approximate")
    risk_flags = row.get("risk_flags", "")
    for flag in ["anage_alias_or_subspecies", "trinomial_or_subspecies", "anage_quality_questionable"]:
        if flag in risk_flags:
            reasons.append(flag)
    return ";".join(reasons) if reasons else "none"


def build_rows(pgls_rows: list[dict[str, str]], tnrs_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    tnrs_by_search = {row["tree_search_name"]: row for row in tnrs_rows}
    output = []
    for row in pgls_rows:
        match = tnrs_by_search.get(row["tree_search_name"], {})
        final_name = match.get("matched_name") or row["tree_search_name"]
        joined = dict(row)
        joined.update(
            {
                "ott_id": match.get("ott_id", ""),
                "opentree_matched_name": final_name,
                "opentree_unique_name": match.get("unique_name", ""),
                "opentree_rank": match.get("rank", ""),
                "opentree_score": match.get("score", ""),
                "opentree_is_synonym": match.get("is_synonym", ""),
                "opentree_is_approximate": match.get("is_approximate_match", ""),
                "final_tree_label": tree_label(final_name),
                "pgls_tree_status": "ready" if match.get("match_status") == "matched" else "unmatched",
            }
        )
        joined["manual_review_reason"] = review_reason(joined)
        output.append(joined)
    output.sort(key=lambda item: (item["clade"], item["final_tree_label"]))
    return output


def write_report(report: pathlib.Path, rows: list[dict[str, str]]) -> None:
    status = Counter(row["pgls_tree_status"] for row in rows)
    ranks = Counter(row["opentree_rank"] for row in rows)
    synonyms = [row for row in rows if row["opentree_is_synonym"] == "True"]
    manual = [row for row in rows if row["manual_review_reason"] != "none"]
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "\n".join(
            [
                "# OpenTree TNRS Report",
                "",
                f"Rows: {len(rows)}",
                "",
                "## Match Status",
                *[f"- {key}: {value}" for key, value in sorted(status.items())],
                "",
                "## OpenTree Rank",
                *[f"- {key}: {value}" for key, value in sorted(ranks.items())],
                "",
                "## Synonym Matches",
                *[
                    f"- {row['tree_search_name']} -> {row['opentree_matched_name']} ({row['opentree_rank']})"
                    for row in synonyms
                ],
                "",
                "## Manual Review",
                f"Rows flagged for manual review before final tree pruning: {len(manual)}",
                "",
                "Main review reasons are held in `data/processed/pgls_tree_manual_review.tsv`.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgls", type=pathlib.Path, default=DEFAULT_PGLS)
    parser.add_argument("--tnrs", type=pathlib.Path, default=DEFAULT_TNRS)
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--review-output", type=pathlib.Path, default=DEFAULT_REVIEW)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    rows = build_rows(read_tsv(args.pgls), read_tsv(args.tnrs))
    fields = list(rows[0].keys())
    write_tsv(args.output, rows, fields)
    review_rows = [row for row in rows if row["manual_review_reason"] != "none"]
    write_tsv(args.review_output, review_rows, fields)
    write_report(args.report, rows)
    print(f"Wrote {args.output}, {args.review_output}, and {args.report}")


if __name__ == "__main__":
    main()
