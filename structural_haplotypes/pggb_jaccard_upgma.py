#!/usr/bin/env python3
"""Build the structural-haplotype UPGMA tree from PGGB/ODGI Jaccard distances."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform

SAMPLE_ALIASES = {
    "YAO.Mat": "YAO.hap2",
    "YAO.Pat": "YAO.hap1",
    "CHM13v2.hap0": "CHM13v2",
    "HG002v1.1.hap1": "HG002.hap1",
    "HG002v1.1.hap2": "HG002.hap2",
    "CN1v1.hap0": "CN1v1",
    "GRCh38.hap0": "GRCh38",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--distance",
        required=True,
        type=Path,
        help="TSV from odgi similarity -d with group.a, group.b and jaccard.distance.",
    )
    parser.add_argument(
        "--haplotypes",
        required=True,
        type=Path,
        help="Structural-haplotype TSV containing Sample and Hap columns.",
    )
    parser.add_argument("--output", required=True, type=Path, help="Output Newick tree.")
    return parser.parse_args()


def normalize_sample_hap(sample: str, hap: str) -> str:
    sample_hap = f"{sample}.{hap}"
    sample_hap = SAMPLE_ALIASES.get(sample_hap, sample_hap)
    if sample_hap.startswith("apr"):
        sample_hap = sample_hap.replace("hap", "")
    return sample_hap


def odgi_key(value: str) -> str:
    """Return the Sample#Hap key from an ODGI path/group label."""
    head = str(value).split(":", 1)[0]
    fields = head.split("#")
    if len(fields) < 2:
        raise ValueError(f"Cannot parse Sample#Hap from ODGI group label: {value}")
    return "#".join(fields[:2])


def scipy_tree_to_newick(node, parent_distance: float, leaf_names: np.ndarray) -> str:
    branch = max(0.0, parent_distance - node.dist)
    if node.is_leaf():
        return f"{leaf_names[node.id]}:{branch:.6g}"

    left = scipy_tree_to_newick(node.get_left(), node.dist, leaf_names)
    right = scipy_tree_to_newick(node.get_right(), node.dist, leaf_names)
    return f"({left},{right}):{branch:.6g}"


def main() -> None:
    args = parse_args()

    haplotypes = pd.read_csv(args.haplotypes, sep="\t", dtype=str)
    required_hap = {"Sample", "Hap"}
    missing = required_hap - set(haplotypes.columns)
    if missing:
        raise ValueError(f"Missing columns in haplotype table: {sorted(missing)}")

    key_to_label = {}
    for sample, hap in zip(haplotypes["Sample"], haplotypes["Hap"]):
        key = f"{sample}#{hap}"
        label = normalize_sample_hap(str(sample), str(hap))
        if key in key_to_label and key_to_label[key] != label:
            raise ValueError(f"Conflicting label for {key}")
        key_to_label[key] = label

    labels = list(key_to_label.values())
    if len(labels) != len(set(labels)):
        duplicated = pd.Series(labels)[pd.Series(labels).duplicated()].unique().tolist()
        raise ValueError(f"Duplicated normalized haplotype identifiers: {duplicated[:10]}")

    distances = pd.read_csv(args.distance, sep="\t")
    required_dist = {"group.a", "group.b", "jaccard.distance"}
    missing = required_dist - set(distances.columns)
    if missing:
        raise ValueError(f"Missing columns in ODGI distance table: {sorted(missing)}")

    matrix = pd.DataFrame(np.nan, index=labels, columns=labels, dtype=float)
    np.fill_diagonal(matrix.values, 0.0)

    for _, row in distances.iterrows():
        key_a = odgi_key(row["group.a"])
        key_b = odgi_key(row["group.b"])
        if key_a not in key_to_label or key_b not in key_to_label:
            continue

        sample_a = key_to_label[key_a]
        sample_b = key_to_label[key_b]
        distance = float(row["jaccard.distance"])
        if distance < 0:
            raise ValueError(f"Negative Jaccard distance for {sample_a}, {sample_b}: {distance}")

        previous = matrix.at[sample_a, sample_b]
        if np.isfinite(previous) and not np.isclose(previous, distance):
            raise ValueError(
                f"Conflicting distances for {sample_a}, {sample_b}: {previous} vs {distance}"
            )
        matrix.at[sample_a, sample_b] = distance
        matrix.at[sample_b, sample_a] = distance

    values = matrix.to_numpy(dtype=float)
    missing_pairs = np.argwhere(np.triu(np.isnan(values), k=1))
    if len(missing_pairs):
        examples = [f"{labels[i]} <> {labels[j]}" for i, j in missing_pairs[:10]]
        raise ValueError(
            f"Distance matrix is incomplete: {len(missing_pairs)} pair(s) missing. "
            f"Examples: {', '.join(examples)}"
        )

    if not np.allclose(values, values.T):
        raise ValueError("Jaccard distance matrix is not symmetric")

    condensed = squareform(values, checks=True)
    linkage_matrix = linkage(condensed, method="average")
    root = to_tree(linkage_matrix, rd=False)
    newick = scipy_tree_to_newick(root, root.dist, np.asarray(labels, dtype=object)) + ";\n"

    args.output.write_text(newick, encoding="utf-8")
    print(f"Haplotypes: {len(labels)}")
    print(f"Wrote UPGMA tree: {args.output}")


if __name__ == "__main__":
    main()
