#!/usr/bin/env python3
"""Benchmark tree-only clustering of PGA upstream/downstream flanking trees."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean
from typing import Dict, List

import hdbscan
import numpy as np
import pandas as pd
from ete3 import Tree
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score

OUTGROUP = "chimpanzee"
HDBSCAN_SIZES = [5, 10, 20, 40]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ut-tree", required=True, type=Path)
    parser.add_argument("--dt-tree", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def normalize_hap(value: str) -> str:
    text = str(value).strip().replace("HG002v1.1", "HG002")
    if text in {"GRCh38", "CHM13v2", "CN1v1"}:
        return f"{text}.hap0"
    if not text.endswith((".hap0", ".hap1", ".hap2")):
        parts = text.rsplit(".", 1)
        if len(parts) == 2 and parts[1] in {"1", "2"}:
            text = f"{parts[0]}.hap{parts[1]}"
    return text


def annotation_table(path: Path) -> Dict[str, dict]:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if {"Sample", "Hap"}.issubset(df.columns):
        df["tree_id"] = [normalize_hap(f"{s}.{h}") for s, h in zip(df["Sample"], df["Hap"])]
    elif "sample_hap" in df.columns:
        df["tree_id"] = df["sample_hap"].map(normalize_hap)
    else:
        raise ValueError("Annotation table must contain Sample/Hap or sample_hap")
    return {row["tree_id"]: row for row in df.to_dict("records")}


def read_tree(path: Path) -> Tree:
    try:
        return Tree(str(path), format=1)
    except Exception:
        return Tree(str(path))


def patristic_matrix(tree: Tree, tips: List[str]) -> np.ndarray:
    nodes = list(tree.traverse("preorder"))
    node_id = {node: i for i, node in enumerate(nodes)}
    adjacency = [[] for _ in nodes]
    for node in nodes:
        if node.up is None:
            continue
        a, b = node_id[node], node_id[node.up]
        w = float(getattr(node, "dist", 0.0) or 0.0)
        adjacency[a].append((b, w))
        adjacency[b].append((a, w))

    leaf_id = {leaf.name: node_id[leaf] for leaf in tree.iter_leaves()}
    matrix = np.zeros((len(tips), len(tips)), dtype=float)
    for i, tip in enumerate(tips):
        dist = np.full(len(nodes), np.inf)
        start = leaf_id[tip]
        dist[start] = 0.0
        stack = [(start, -1)]
        while stack:
            current, parent = stack.pop()
            for nxt, weight in adjacency[current]:
                if nxt == parent:
                    continue
                dist[nxt] = dist[current] + weight
                stack.append((nxt, current))
        for j, other in enumerate(tips):
            matrix[i, j] = dist[leaf_id[other]]
    return matrix


def relabel(raw, prefix: str) -> List[str]:
    counts = Counter(x for x in raw if int(x) >= 0)
    ordered = sorted(counts, key=lambda x: (-counts[x], x))
    mapping = {label: f"{prefix}{i+1:04d}" for i, label in enumerate(ordered)}
    return [mapping.get(x, "NA") for x in raw]


def metrics(tree: Tree, tips: List[str], matrix: np.ndarray, labels: List[str]) -> dict:
    assigned = [i for i, label in enumerate(labels) if label != "NA"]
    clusters = defaultdict(list)
    for i in assigned:
        clusters[labels[i]].append(i)

    mono = []
    tip_set = set(tips)
    for indices in clusters.values():
        names = {tips[i] for i in indices}
        if len(names) == 1:
            mono.append(1)
        else:
            mrca = tree.get_common_ancestor(list(names))
            descendants = {x.name for x in mrca.iter_leaves() if x.name in tip_set}
            mono.append(int(descendants == names))

    silhouette = np.nan
    if len(assigned) >= 3 and 1 < len(clusters) < len(assigned):
        sub = matrix[np.ix_(assigned, assigned)]
        silhouette = silhouette_score(sub, [labels[i] for i in assigned], metric="precomputed")

    small = sum(len(v) for v in clusters.values() if len(v) < 5)
    coverage = len(assigned) / len(tips) if tips else 0.0
    return {
        "coverage": coverage,
        "cluster_count": len(clusters),
        "small_cluster_fraction": small / len(assigned) if assigned else 1.0,
        "monophyletic_cluster_fraction": mean(mono) if mono else 0.0,
        "silhouette": silhouette,
    }


def score(row: dict) -> float:
    sil = row["silhouette"] if np.isfinite(row["silhouette"]) else 0.0
    return (
        0.35 * row["coverage"]
        + 0.35 * row["monophyletic_cluster_fraction"]
        + 0.15 * max(sil, 0.0)
        - 0.15 * row["small_cluster_fraction"]
        - 0.001 * row["cluster_count"]
    )


def hdbscan_results(tree_name, tree, tips, matrix):
    out = []
    for size in HDBSCAN_SIZES:
        raw = hdbscan.HDBSCAN(
            metric="precomputed",
            min_cluster_size=size,
            min_samples=None,
            cluster_selection_method="eom",
        ).fit_predict(matrix)
        labels = relabel(raw, "HD")
        row = {"tree": tree_name, "method": "hdbscan", "parameter": f"min_cluster_size={size}"}
        row.update(metrics(tree, tips, matrix, labels))
        row["score"] = score(row)
        out.append((row, labels))
    return out


def agglomerative_results(tree_name, tree, tips, matrix):
    positive = matrix[np.triu_indices_from(matrix, 1)]
    positive = positive[positive > 0]
    if len(positive) == 0:
        return []
    thresholds = sorted(set(float(np.quantile(positive, q)) for q in [0.01, 0.02, 0.05, 0.1, 0.2]))
    out = []
    for threshold in thresholds:
        model = AgglomerativeClustering(
            n_clusters=None,
            metric="precomputed",
            linkage="average",
            distance_threshold=threshold,
        )
        labels = relabel(model.fit_predict(matrix), "AG")
        row = {"tree": tree_name, "method": "agglomerative", "parameter": f"distance={threshold:g}"}
        row.update(metrics(tree, tips, matrix, labels))
        row["score"] = score(row)
        out.append((row, labels))
    return out


def treecluster_results(tree_name, tree, tips, matrix, tree_path: Path):
    executable = shutil.which("TreeCluster.py") or shutil.which("TreeCluster")
    if executable is None:
        return []

    positive = matrix[np.triu_indices_from(matrix, 1)]
    positive = positive[positive > 0]
    thresholds = sorted(set(float(np.quantile(positive, q)) for q in [0.01, 0.02, 0.05, 0.1, 0.2]))
    out = []
    for method in ["max_clade", "avg_clade"]:
        for threshold in thresholds:
            with tempfile.NamedTemporaryFile(suffix=".tsv") as tmp:
                launcher = [sys.executable, executable] if executable.endswith(".py") else [executable]
                cmd = launcher + ["-i", str(tree_path), "-t", str(threshold), "-m", method, "-o", tmp.name]
                subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                df = pd.read_csv(tmp.name, sep="\t", dtype=str)
            name_col = "SequenceName" if "SequenceName" in df.columns else "Sequence"
            cluster_col = "ClusterNumber" if "ClusterNumber" in df.columns else "Cluster"
            mapping = dict(zip(df[name_col], df[cluster_col]))
            raw = [int(mapping.get(tip, -1)) for tip in tips]
            labels = relabel(raw, "TC")
            row = {"tree": tree_name, "method": "treecluster", "parameter": f"{method};distance={threshold:g}"}
            row.update(metrics(tree, tips, matrix, labels))
            row["score"] = score(row)
            out.append((row, labels))
    return out


def write_tsv(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    annotations = annotation_table(args.annotation)

    all_rows = []
    assignments = {}
    for tree_name, tree_path in [("UT", args.ut_tree), ("DT", args.dt_tree)]:
        tree = read_tree(tree_path)
        tips = [leaf.name for leaf in tree.iter_leaves() if leaf.name != OUTGROUP]
        matrix = patristic_matrix(tree, tips)

        results = []
        results.extend(hdbscan_results(tree_name, tree, tips, matrix))
        results.extend(agglomerative_results(tree_name, tree, tips, matrix))
        results.extend(treecluster_results(tree_name, tree, tips, matrix, tree_path))
        all_rows.extend(row for row, _ in results)

        by_method = defaultdict(list)
        for row, labels in results:
            by_method[row["method"]].append((row, labels))

        for method, candidates in by_method.items():
            if method == "hdbscan":
                chosen = next(x for x in candidates if x[0]["parameter"] == "min_cluster_size=5")
            else:
                chosen = max(candidates, key=lambda x: x[0]["score"])
            row, labels = chosen
            assignments[(tree_name, method)] = (row, tips, labels)

    write_tsv(args.outdir / "parameter_scan.tsv", all_rows)

    comparison = []
    for (tree_name, method), (row, _, _) in sorted(assignments.items()):
        comparison.append(row)
    write_tsv(args.outdir / "method_comparison.tsv", comparison)

    methods = sorted({method for _, method in assignments})
    for method in methods:
        rows = []
        for tree_name in ["UT", "DT"]:
            if (tree_name, method) not in assignments:
                continue
            selected, tips, labels = assignments[(tree_name, method)]
            for tip, label in zip(tips, labels):
                hap = normalize_hap(tip)
                ann = annotations.get(hap, {})
                rows.append(
                    {
                        "tree": tree_name,
                        "method": method,
                        "parameter": selected["parameter"],
                        "tree_tip": tip,
                        "haplotype": hap,
                        "cluster_id": label,
                        "final_hap_label": ann.get("final_hap_label", ""),
                        "CN": ann.get("CN", ""),
                        "Superpop": ann.get("Superpop", ""),
                        "New_superpop": ann.get("New_superpop", ""),
                        "Pop": ann.get("Pop", ""),
                    }
                )
        write_tsv(args.outdir / f"cluster_assignments.{method}.tsv", rows)

    print(f"Wrote clustering benchmark to {args.outdir}")


if __name__ == "__main__":
    main()
