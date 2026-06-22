"""Fetch an OpenTree induced subtree for PGLS-ready species."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re
import urllib.error
import urllib.request


DEFAULT_INPUT = pathlib.Path("data/processed/pgls_tree_ready.tsv")
DEFAULT_TREE = pathlib.Path("data/processed/phylogeny_inputs/opentree_induced_subtree.tre")
DEFAULT_TRAITS = pathlib.Path("data/processed/pgls_trait_table.tsv")
DEFAULT_AUDIT = pathlib.Path("data/processed/opentree_subtree_audit.tsv")
DEFAULT_RESPONSE = pathlib.Path("data/processed/phylogeny_inputs/opentree_induced_subtree_response.json")
DEFAULT_REPORT = pathlib.Path("results/reports/opentree_subtree_report.md")
DEFAULT_ENDPOINT = "https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree"


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def post_subtree(endpoint: str, ott_ids: list[int], timeout: int) -> dict:
    request = urllib.request.Request(
        endpoint,
        data=json.dumps({"ott_ids": ott_ids}).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "User-Agent": "bird-lifespan-pgls-prep/0.1",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"OpenTree induced subtree HTTP {exc.code}: {body}") from exc


def extract_tip_labels(newick: str) -> list[str]:
    # OpenTree labels are emitted before ':branch_length', ',' or ')'.
    labels = []
    token = []
    in_quote = False
    for char in newick:
        if char == "'":
            in_quote = not in_quote
            continue
        if not in_quote and char in "(),:;":
            text = "".join(token).strip()
            if text and not re.fullmatch(r"[0-9. Ee+-]+", text):
                labels.append(text)
            token = []
        else:
            token.append(char)
    return labels


def label_ott_id(label: str) -> str:
    match = re.search(r"_ott(\d+)$", label)
    return match.group(1) if match else ""


def build_outputs(
    rows: list[dict[str, str]],
    response: dict,
) -> tuple[str, list[dict[str, str]], list[dict[str, str]]]:
    newick = response.get("newick", "")
    if not newick:
        raise RuntimeError(f"OpenTree response did not include newick. Keys: {sorted(response)}")

    tip_labels = extract_tip_labels(newick)
    tip_by_ott = {label_ott_id(label): label for label in tip_labels if label_ott_id(label)}
    not_in_tree = {str(item) for item in response.get("ott_ids_not_in_tree", []) or []}
    node_ids_not_in_tree = {str(item) for item in response.get("node_ids_not_in_tree", []) or []}
    broken = {str(key).replace("ott", ""): value for key, value in (response.get("broken", {}) or {}).items()}

    trait_rows = []
    audit_rows = []
    for row in rows:
        ott_id = row["ott_id"]
        tip_label = tip_by_ott.get(ott_id, "")
        status = "in_tree" if tip_label else "not_in_tree"
        broken_node = broken.get(ott_id, "")
        if broken_node:
            status = "broken_to_internal_node"
        elif ott_id in not_in_tree or ott_id in node_ids_not_in_tree:
            status = "reported_not_in_tree"
        audit_rows.append(
            {
                "scientific_name": row["scientific_name"],
                "ott_id": ott_id,
                "opentree_matched_name": row["opentree_matched_name"],
                "final_tree_label": row["final_tree_label"],
                "opentree_tip_label": tip_label,
                "subtree_status": status,
                "broken_node": broken_node,
                "manual_review_reason": row["manual_review_reason"],
            }
        )
        if tip_label:
            trait = dict(row)
            trait["opentree_tip_label"] = tip_label
            trait_rows.append(trait)

    return newick, trait_rows, audit_rows


def write_report(
    report: pathlib.Path,
    requested_n: int,
    trait_rows: list[dict[str, str]],
    audit_rows: list[dict[str, str]],
    response: dict,
) -> None:
    missing = [row for row in audit_rows if row["subtree_status"] != "in_tree"]
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "\n".join(
            [
                "# OpenTree Subtree Report",
                "",
                f"Requested OTT IDs: {requested_n}",
                f"Tips recovered into trait table: {len(trait_rows)}",
                f"Missing from induced subtree: {len(missing)}",
                "",
                "## Response Metadata",
                f"- ott_ids_not_in_tree: {response.get('ott_ids_not_in_tree', [])}",
                f"- node_ids_not_in_tree: {response.get('node_ids_not_in_tree', [])}",
                f"- broken taxa: {response.get('broken', {})}",
                "",
                "## Outputs",
                "- `data/processed/phylogeny_inputs/opentree_induced_subtree.tre`",
                "- `data/processed/pgls_trait_table.tsv`",
                "- `data/processed/opentree_subtree_audit.tsv`",
                "",
                "## Next Step",
                "Run a first PGLS model on `pgls_trait_table.tsv` and the induced subtree, then compare against OLS on the same species set.",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, default=DEFAULT_INPUT)
    parser.add_argument("--tree-output", type=pathlib.Path, default=DEFAULT_TREE)
    parser.add_argument("--traits-output", type=pathlib.Path, default=DEFAULT_TRAITS)
    parser.add_argument("--audit-output", type=pathlib.Path, default=DEFAULT_AUDIT)
    parser.add_argument("--response-output", type=pathlib.Path, default=DEFAULT_RESPONSE)
    parser.add_argument("--report", type=pathlib.Path, default=DEFAULT_REPORT)
    parser.add_argument("--endpoint", default=DEFAULT_ENDPOINT)
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()

    rows = read_tsv(args.input)
    ott_ids = [int(row["ott_id"]) for row in rows if row.get("ott_id", "").isdigit()]
    response = post_subtree(args.endpoint, ott_ids, args.timeout)
    args.response_output.parent.mkdir(parents=True, exist_ok=True)
    args.response_output.write_text(json.dumps(response, indent=2, sort_keys=True), encoding="utf-8")

    newick, trait_rows, audit_rows = build_outputs(rows, response)
    args.tree_output.parent.mkdir(parents=True, exist_ok=True)
    args.tree_output.write_text(newick + "\n", encoding="utf-8")
    write_tsv(args.traits_output, trait_rows, list(trait_rows[0].keys()))
    write_tsv(
        args.audit_output,
        audit_rows,
        [
            "scientific_name",
            "ott_id",
            "opentree_matched_name",
            "final_tree_label",
            "opentree_tip_label",
            "subtree_status",
            "broken_node",
            "manual_review_reason",
        ],
    )
    write_report(args.report, len(ott_ids), trait_rows, audit_rows, response)
    print(f"Wrote {args.tree_output}, {args.traits_output}, and {args.audit_output}")


if __name__ == "__main__":
    main()
