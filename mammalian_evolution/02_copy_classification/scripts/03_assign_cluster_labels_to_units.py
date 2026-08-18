#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filter all.VPS37C-VWCE.locus_annotations.bed by complete-unit representative IDs,
replace BED name/header column with cluster_best label from cluster_assignments.updated.tsv,
sort by chrom/start/end, then add per-label suffixes:
  PGA_1, PGA_2, PAG_1, PGA-like_1, ...

Default assumptions:
  - BED name/header is column 4.
  - ID list entries match the BED name/header column exactly.
  - cluster_assignments.updated.tsv contains columns: id, header, cluster_best.
  - cluster_best is already updated to labels: PGA/PAG/PGA-like.
"""
import argparse
import csv
import re
import sys
from collections import defaultdict


def chrom_key(chrom):
    """Lexicographic-like key with embedded numbers handled naturally within the chromosome string."""
    parts = re.split(r"(\d+)", chrom)
    return tuple(int(p) if p.isdigit() else p for p in parts)


def load_ids(path):
    ids = set()
    with open(path) as f:
        for line in f:
            x = line.strip()
            if x and not x.startswith("#"):
                ids.add(x)
    if not ids:
        sys.exit(f"ERROR: no IDs loaded from {path}")
    return ids


def load_assignment_map(path, label_col="cluster_best"):
    """
    Return a mapping from both id and header to cluster_best label.
    This makes the script work whether the BED name column contains 'id' or 'header'.
    """
    m = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            sys.exit(f"ERROR: empty assignment file: {path}")
        required = {"id", "header", label_col}
        missing = required - set(reader.fieldnames)
        if missing:
            sys.exit(f"ERROR: missing columns in {path}: {', '.join(sorted(missing))}")

        for row_no, row in enumerate(reader, start=2):
            label = row[label_col].strip()
            if not label:
                sys.exit(f"ERROR: line {row_no} has empty {label_col}")
            for key_col in ("id", "header"):
                key = row.get(key_col, "").strip()
                if key:
                    # If duplicated keys have different labels, stop rather than silently overwrite.
                    if key in m and m[key] != label:
                        sys.exit(
                            f"ERROR: duplicated key with conflicting labels in {path}: {key!r}: {m[key]!r} vs {label!r}"
                        )
                    m[key] = label
    if not m:
        sys.exit(f"ERROR: no assignment records loaded from {path}")
    return m


def parse_int(x, line_no, col_name):
    try:
        return int(x)
    except ValueError:
        sys.exit(f"ERROR: line {line_no}: BED {col_name} is not an integer: {x!r}")


def main():
    ap = argparse.ArgumentParser(
        description="Create complete-unit BED with cluster_best labels and ordered suffixes."
    )
    ap.add_argument("--assign", required=True, help="cluster_assignments.updated.tsv")
    ap.add_argument("--ids", required=True, help="toga.PGA_like.local.v3.complete_unit_representative_ids.txt")
    ap.add_argument("--bed", required=True, help="all.VPS37C-VWCE.locus_annotations.bed")
    ap.add_argument("-o", "--output", required=True, help="Output BED")
    ap.add_argument(
        "--bed-name-col",
        type=int,
        default=4,
        help="1-based BED column containing representative id/header to filter and replace. Default: 4",
    )
    ap.add_argument(
        "--label-col",
        default="cluster_best",
        help="Label column in updated assignment table. Default: cluster_best",
    )
    ap.add_argument(
        "--filter-any-column",
        action="store_true",
        help="Mimic grep more broadly: keep a BED line if any tab-delimited field exactly matches an ID. Replacement still uses --bed-name-col.",
    )
    ap.add_argument(
        "--keep-original-name",
        action="store_true",
        help="Append original BED name/header as the last column for traceability. Default: do not append.",
    )
    args = ap.parse_args()

    name_idx = args.bed_name_col - 1
    if name_idx < 0:
        sys.exit("ERROR: --bed-name-col must be >= 1")

    ids = load_ids(args.ids)
    assign_map = load_assignment_map(args.assign, args.label_col)

    records = []
    n_bed = 0
    n_keep = 0
    n_missing_assignment = 0

    with open(args.bed) as f:
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            row = line.split("\t")
            n_bed += 1
            if len(row) <= max(2, name_idx):
                sys.exit(f"ERROR: line {line_no} has too few BED columns for --bed-name-col {args.bed_name_col}: {line}")

            original_name = row[name_idx]
            if args.filter_any_column:
                keep = any(field in ids for field in row)
            else:
                keep = original_name in ids
            if not keep:
                continue

            if original_name not in assign_map:
                n_missing_assignment += 1
                sys.exit(
                    f"ERROR: BED name/header {original_name!r} from line {line_no} is not found in id/header columns of {args.assign}. "
                    "Check whether --bed-name-col is correct, or whether the ID list and assignment table use different naming."
                )

            label = assign_map[original_name]
            start = parse_int(row[1], line_no, "start/column2")
            end = parse_int(row[2], line_no, "end/column3")
            records.append((chrom_key(row[0]), start, end, row[0], start, end, label, original_name, row))
            n_keep += 1

    records.sort(key=lambda x: (x[0], x[1], x[2], x[7]))

    label_counts = defaultdict(int)
    with open(args.output, "w") as out:
        for _, _, _, _chrom, _start, _end, label, original_name, row in records:
            label_counts[label] += 1
            row[name_idx] = f"{label}_{label_counts[label]}"
            if args.keep_original_name:
                row = row + [original_name]
            out.write("\t".join(row) + "\n")

    print(f"Wrote: {args.output}", file=sys.stderr)
    print(f"BED records read: {n_bed}", file=sys.stderr)
    print(f"BED records kept: {n_keep}", file=sys.stderr)
    print("Output label counts:", file=sys.stderr)
    for label in sorted(label_counts):
        print(f"  {label}\t{label_counts[label]}", file=sys.stderr)


if __name__ == "__main__":
    main()
