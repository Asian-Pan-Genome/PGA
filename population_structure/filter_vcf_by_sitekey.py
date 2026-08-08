#!/usr/bin/env python3

import argparse
import gzip
import sys


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_sites(path):
    sites = set()
    with open(path) as fin:
        for line in fin:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"Malformed site key: {line.rstrip()}")
            sites.add(tuple(fields[:4]))
    return sites


def main():
    parser = argparse.ArgumentParser(
        description="Keep VCF records matching exact CHROM/POS/REF/ALT site keys."
    )
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--sites", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    sites = load_sites(args.sites)
    kept = 0
    total = 0

    with open_text(args.vcf) as fin, open(args.out, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            total += 1
            key = (fields[0], fields[1], fields[3], fields[4])
            if key in sites:
                fout.write(line)
                kept += 1

    print(f"Kept {kept} / {total} variants from {args.vcf}", file=sys.stderr)


if __name__ == "__main__":
    main()
