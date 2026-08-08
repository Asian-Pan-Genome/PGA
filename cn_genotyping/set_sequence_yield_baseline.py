#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd


DEFAULT_GENOME_SIZE = 3.1e9
BWA_RE = re.compile(r"read\s+\d+\s+sequences\s+\((\d+)\s+bp\)")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Replace Baseline_depth using total input sequence yield."
    )
    parser.add_argument("--depth", required=True, help="Depth TSV produced by the PGA extraction step.")
    parser.add_argument("--out", required=True, help="Output depth TSV.")
    parser.add_argument("--genome-size", type=float, default=DEFAULT_GENOME_SIZE)

    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--bwa-log",
        help="stderr log from the initial BWA mapping; all '[M::process] read ... (N bp)' entries are summed.",
    )
    source.add_argument(
        "--metadata",
        help="TSV containing per-sample total sequenced bases.",
    )

    parser.add_argument("--sample", help="Sample ID for --metadata mode. Defaults to the depth-file prefix.")
    parser.add_argument("--sample-column", default="run_accession")
    parser.add_argument("--base-count-column", default="base_count")
    return parser.parse_args()


def bases_from_bwa_log(path):
    total = 0
    with open(path) as handle:
        for line in handle:
            match = BWA_RE.search(line)
            if match:
                total += int(match.group(1))
    if total == 0:
        raise ValueError(f"No BWA '[M::process] read ... (N bp)' entries found in {path}")
    return total


def bases_from_metadata(path, sample, sample_column, base_count_column):
    table = pd.read_csv(path, sep="\t")
    for column in (sample_column, base_count_column):
        if column not in table.columns:
            raise ValueError(f"Missing column '{column}' in {path}")

    match = table.loc[table[sample_column].astype(str) == str(sample), base_count_column]
    if len(match) != 1:
        raise ValueError(f"Expected one metadata row for sample '{sample}', found {len(match)}")
    return float(match.iloc[0])


def main():
    args = parse_args()

    if args.bwa_log:
        base_count = bases_from_bwa_log(args.bwa_log)
    else:
        sample = args.sample or os.path.basename(args.depth).split(".")[0]
        base_count = bases_from_metadata(
            args.metadata,
            sample,
            args.sample_column,
            args.base_count_column,
        )

    depth = pd.read_csv(args.depth, sep="\t")
    if "Baseline_depth" not in depth.columns:
        raise ValueError(f"{args.depth} does not contain a Baseline_depth column")

    depth["Baseline_depth"] = base_count / args.genome_size
    depth.to_csv(args.out, sep="\t", index=False)
    print(f"Baseline depth: {base_count / args.genome_size:.6f}")
    print(f"Output: {args.out}")


if __name__ == "__main__":
    main()
