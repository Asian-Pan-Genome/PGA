#!/usr/bin/env python3
"""Refine PGA flanking-tree clusters with SH-enriched monophyletic clades."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Dict, List, Tuple

import pandas as pd
from ete3 import Tree
from sklearn.metrics import adjusted_rand_score

FINAL_METHOD = "hdbscan"
FINAL_PURITY = 0.80
FINAL_MIN_SIZE = 10
PURITY_GRID = [0.60, 0.70, 0.80, 0.90, 0.95]
MIN_SIZE_GRID = [5, 10, 20]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--ut-tree", required=True, type=Path)
    parser.add_argument("--dt-tree", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def read_tree(path: Path) -> Tree:
    try:
        return Tree(str(path), format=1)
    except Exception:
        return Tree(str(path))


def clean(value) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text in {"", "NA", "nan", "None", "-1"} else text


def majority(labels: List[str]) -> Tuple[str, int, float]:
    counts = Counter(label if label else "Unknown" for label in labels)
    label, count = sorted(counts.items(), key=lambda x: (-x[1], x[0]))[0]
    return label, count, count / len(labels)


def refine(tree: Tree, df: pd.DataFrame, tree_name: str, method: str, purity: float, min_size: int):
    records = df.to_dict("records")
    by_tip = {row["tree_tip"]: row for row in records}
    valid = set(by_tip)
    raw = {tip: clean(row.get("cluster_id")) for tip, row in by_tip.items()}
    sh = {tip: clean(row.get("final_hap_label")) or "Unknown" for tip, row in by_tip.items()}

    candidates = []
    for node in tree.traverse("preorder"):
        if node.is_leaf():
            continue
        tips = [leaf.name for leaf in node.iter_leaves() if leaf.name in valid]
        if len(tips) < min_size:
            continue
        raw_clusters = {raw[tip] for tip in tips if raw[tip]}
        if len(raw_clusters) != 1:
            continue
        major_sh, major_count, group_purity = majority([sh[tip] for tip in tips])
        if group_purity < purity:
            continue
        candidates.append(
            {
                "tips": tips,
                "tip_set": set(tips),
                "parent_cluster": next(iter(raw_clusters)),
                "major_SH": major_sh,
                "major_count": major_count,
                "purity": group_purity,
                "size": len(tips),
            }
        )

    candidates.sort(key=lambda x: (-x["size"], -x["purity"], sorted(x["tips"])[0]))
    used = set()
    assignments = {}
    groups = []
    group_no = 0

    for candidate in candidates:
        if candidate["tip_set"] & used:
            continue
        group_no += 1
        group_id = f"{tree_name}_{method}_R{group_no:03d}"
        groups.append(
            {
                "tree": tree_name,
                "method": method,
                "refined_group": group_id,
                "parent_cluster": candidate["parent_cluster"],
                "major_SH": candidate["major_SH"],
                "group_size": candidate["size"],
                "major_SH_count": candidate["major_count"],
                "purity": candidate["purity"],
                "purity_threshold": purity,
                "min_size": min_size,
            }
        )
        for tip in candidate["tips"]:
            row = by_tip[tip]
            hap = row["haplotype"]
            assignments[hap] = {
                **row,
                "refined_group": group_id,
                "parent_cluster": candidate["parent_cluster"],
                "major_SH": candidate["major_SH"],
                "member_SH": sh[tip],
                "group_size": candidate["size"],
                "purity": candidate["purity"],
                "purity_threshold": purity,
                "min_size": min_size,
            }
        used |= candidate["tip_set"]

    return assignments, groups


def pair_labels(ut: Dict[str, dict], dt: Dict[str, dict]) -> Dict[str, str]:
    out = {}
    for hap in set(ut) | set(dt):
        if hap not in ut or hap not in dt:
            out[hap] = "NA"
        else:
            out[hap] = f"{ut[hap]['major_SH']}->{dt[hap]['major_SH']}"
    return out


def stability(current: Dict[str, str], neighbors: List[Dict[str, str]]) -> float:
    scores = []
    for other in neighbors:
        haps = sorted(set(current) | set(other))
        if len(haps) < 2:
            continue
        scores.append(
            adjusted_rand_score(
                [current.get(hap, "NA") for hap in haps],
                [other.get(hap, "NA") for hap in haps],
            )
        )
    return sum(scores) / len(scores) if scores else 1.0


def write_tsv(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_raw(raw_dir: Path):
    result = {}
    for path in sorted(raw_dir.glob("cluster_assignments.*.tsv")):
        method = path.name.split(".")[1]
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        result[method] = df
    if not result:
        raise FileNotFoundError(f"No cluster_assignments.*.tsv in {raw_dir}")
    return result


def method_comparison(raw_by_method, trees):
    rows = []
    final_results = {}
    for method, df in raw_by_method.items():
        method_result = {}
        for tree_name in ["UT", "DT"]:
            tree_df = df[df["tree"] == tree_name].copy()
            assignments, groups = refine(
                trees[tree_name], tree_df, tree_name, method, FINAL_PURITY, FINAL_MIN_SIZE
            )
            method_result[tree_name] = (assignments, groups)
        ut, ut_groups = method_result["UT"]
        dt, dt_groups = method_result["DT"]
        shared = set(ut) & set(dt)
        purities = [g["purity"] for g in ut_groups + dt_groups]
        rows.append(
            {
                "method": method,
                "purity_threshold": FINAL_PURITY,
                "min_size": FINAL_MIN_SIZE,
                "UT_assigned": len(ut),
                "DT_assigned": len(dt),
                "UT_refined_groups": len(ut_groups),
                "DT_refined_groups": len(dt_groups),
                "shared_UT_DT_assigned": len(shared),
                "median_purity": median(purities) if purities else 0.0,
            }
        )
        final_results[method] = method_result
    return rows, final_results


def hdbscan_sensitivity(df: pd.DataFrame, trees):
    scan = {}
    for purity in PURITY_GRID:
        for min_size in MIN_SIZE_GRID:
            results = {}
            for tree_name in ["UT", "DT"]:
                tree_df = df[df["tree"] == tree_name].copy()
                results[tree_name] = refine(
                    trees[tree_name], tree_df, tree_name, FINAL_METHOD, purity, min_size
                )
            scan[(purity, min_size)] = results

    rows = []
    for purity in PURITY_GRID:
        for min_size in MIN_SIZE_GRID:
            current = scan[(purity, min_size)]
            current_pairs = pair_labels(current["UT"][0], current["DT"][0])
            neighbors = []
            pi = PURITY_GRID.index(purity)
            mi = MIN_SIZE_GRID.index(min_size)
            for pidx, midx in [(pi - 1, mi), (pi + 1, mi), (pi, mi - 1), (pi, mi + 1)]:
                if 0 <= pidx < len(PURITY_GRID) and 0 <= midx < len(MIN_SIZE_GRID):
                    neighbor = scan[(PURITY_GRID[pidx], MIN_SIZE_GRID[midx])]
                    neighbors.append(pair_labels(neighbor["UT"][0], neighbor["DT"][0]))

            ut, ut_groups = current["UT"]
            dt, dt_groups = current["DT"]
            rows.append(
                {
                    "method": FINAL_METHOD,
                    "purity_threshold": purity,
                    "min_size": min_size,
                    "UT_assigned": len(ut),
                    "DT_assigned": len(dt),
                    "UT_refined_groups": len(ut_groups),
                    "DT_refined_groups": len(dt_groups),
                    "shared_UT_DT_assigned": len(set(ut) & set(dt)),
                    "pair_stability": stability(current_pairs, neighbors),
                    "selected": int(math.isclose(purity, FINAL_PURITY) and min_size == FINAL_MIN_SIZE),
                }
            )
    return rows, scan[(FINAL_PURITY, FINAL_MIN_SIZE)]


def comparator_candidates(focal: str, ut: Dict[str, dict], dt: Dict[str, dict]):
    u = ut[focal]
    d = dt[focal]
    left = []
    right = []

    for hap in sorted(set(ut) & set(dt)):
        if hap == focal:
            continue
        if (
            ut[hap]["refined_group"] == u["refined_group"]
            and clean(ut[hap]["member_SH"]) == clean(u["major_SH"])
            and dt[hap]["refined_group"] != d["refined_group"]
        ):
            left.append(hap)
        if (
            dt[hap]["refined_group"] == d["refined_group"]
            and clean(dt[hap]["member_SH"]) == clean(d["major_SH"])
            and ut[hap]["refined_group"] != u["refined_group"]
        ):
            right.append(hap)
    return left, right


def candidate_rows(ut: Dict[str, dict], dt: Dict[str, dict]) -> List[dict]:
    rows = []
    for hap in sorted(set(ut) & set(dt)):
        observed = clean(ut[hap].get("final_hap_label")) or clean(ut[hap].get("member_SH"))
        if observed == clean(ut[hap]["major_SH"]) or observed == clean(dt[hap]["major_SH"]):
            continue
        left, right = comparator_candidates(hap, ut, dt)
        rows.append(
            {
                "haplotype": hap,
                "observed_SH": observed,
                "UT_refined_group": ut[hap]["refined_group"],
                "UT_major_SH": ut[hap]["major_SH"],
                "DT_refined_group": dt[hap]["refined_group"],
                "DT_major_SH": dt[hap]["major_SH"],
                "left_context_candidates": ",".join(left),
                "right_context_candidates": ",".join(right),
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    trees = {"UT": read_tree(args.ut_tree), "DT": read_tree(args.dt_tree)}
    raw_by_method = load_raw(args.raw_dir)

    comparison, fixed_results = method_comparison(raw_by_method, trees)
    write_tsv(args.outdir / "method_comparison.SH_refined.tsv", comparison)

    if FINAL_METHOD not in raw_by_method:
        raise ValueError("HDBSCAN cluster assignments were not found in raw-dir")

    sensitivity_rows, final = hdbscan_sensitivity(raw_by_method[FINAL_METHOD], trees)
    write_tsv(args.outdir / "parameter_scan.hdbscan.SH_refined.tsv", sensitivity_rows)

    ut, ut_groups = final["UT"]
    dt, dt_groups = final["DT"]
    assignment_rows = []
    for tree_name, assignments in [("UT", ut), ("DT", dt)]:
        for hap, row in assignments.items():
            assignment_rows.append({"tree": tree_name, **row})
    write_tsv(args.outdir / "refined_cluster_assignments.hdbscan.tsv", assignment_rows)
    write_tsv(args.outdir / "refined_cluster_summary.hdbscan.tsv", ut_groups + dt_groups)
    write_tsv(args.outdir / "ancestral_NAHR_candidates.hdbscan.tsv", candidate_rows(ut, dt))

    print(f"Final method: {FINAL_METHOD}")
    print(f"Final refinement: purity >= {FINAL_PURITY}, min_size >= {FINAL_MIN_SIZE}")
    print(f"Wrote outputs to {args.outdir}")


if __name__ == "__main__":
    main()
