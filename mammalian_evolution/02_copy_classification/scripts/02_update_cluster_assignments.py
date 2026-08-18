#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Update cluster_assignments.tsv: convert numeric cluster_best values to labels.
Default mapping:
  3 -> PGA
  1 -> PAG
  2,4 -> PGA-like

Input columns expected from cluster_assignments.tsv:
  id, header, gene, species, cluster_best, cluster_at_ngenes
"""
import argparse
import csv
import os
import sys


def parse_map(map_items):
    """Parse mapping items like 3=PGA 1=PAG 2=PGA-like 4=PGA-like."""
    m = {}
    for item in map_items:
        if "=" not in item:
            raise ValueError(f"Bad --map item: {item!r}; expected format like 3=PGA")
        k, v = item.split("=", 1)
        k = k.strip()
        v = v.strip()
        if not k or not v:
            raise ValueError(f"Bad --map item: {item!r}; empty key or value")
        m[k] = v
    return m


def main():
    ap = argparse.ArgumentParser(
        description="Replace cluster_best numeric codes in cluster_assignments.tsv with PGA/PAG/PGA-like labels."
    )
    ap.add_argument("-i", "--input", required=True, help="Input cluster_assignments.tsv")
    ap.add_argument("-o", "--output", default=None, help="Output TSV. Default: <input_dir>/cluster_assignments.updated.tsv")
    ap.add_argument(
        "--map",
        nargs="+",
        default=["3=PGA", "1=PAG", "2=PGA-like", "4=PGA-like"],
        help="Cluster mapping, e.g. --map 3=PGA 1=PAG 2=PGA-like 4=PGA-like",
    )
    ap.add_argument(
        "--cluster-col",
        default="cluster_best",
        help="Column to update. Default: cluster_best. If absent, column 5 is used.",
    )
    ap.add_argument(
        "--strict",
        action="store_true",
        help="Error if a cluster_best value is not present in --map. Default: keep original value.",
    )
    args = ap.parse_args()

    if args.output is None:
        args.output = os.path.join(os.path.dirname(os.path.abspath(args.input)), "cluster_assignments.updated.tsv")

    cluster_map = parse_map(args.map)

    n_total = 0
    n_changed = 0
    counts = {}

    with open(args.input, newline="") as fin, open(args.output, "w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")

        try:
            header = next(reader)
        except StopIteration:
            sys.exit(f"ERROR: empty input file: {args.input}")

        if args.cluster_col in header:
            cidx = header.index(args.cluster_col)
        else:
            if len(header) < 5:
                sys.exit("ERROR: input has fewer than 5 columns and cluster_best column was not found.")
            cidx = 4
            print(
                f"WARNING: column {args.cluster_col!r} not found; using column 5 ({header[cidx]!r}) instead.",
                file=sys.stderr,
            )

        writer.writerow(header)
        for row_no, row in enumerate(reader, start=2):
            if not row or all(x == "" for x in row):
                continue
            if len(row) <= cidx:
                sys.exit(f"ERROR: line {row_no} has too few columns: {row}")
            old = row[cidx]
            if old in cluster_map:
                new = cluster_map[old]
                row[cidx] = new
                if new != old:
                    n_changed += 1
            else:
                if args.strict:
                    sys.exit(f"ERROR: line {row_no}: cluster value {old!r} is not in mapping {cluster_map}")
                new = old
            counts[new] = counts.get(new, 0) + 1
            n_total += 1
            writer.writerow(row)

    print(f"Wrote: {args.output}", file=sys.stderr)
    print(f"Rows processed: {n_total}", file=sys.stderr)
    print(f"Rows changed:   {n_changed}", file=sys.stderr)
    print("Label counts:", file=sys.stderr)
    for k in sorted(counts):
        print(f"  {k}\t{counts[k]}", file=sys.stderr)


if __name__ == "__main__":
    main()
