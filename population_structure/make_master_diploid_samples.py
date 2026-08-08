#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(
        description="Keep samples represented by both .1 and .2 assembly haplotypes."
    )
    parser.add_argument("--master-haps", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    hap_by_sample = defaultdict(set)
    with open(args.master_haps) as fin:
        for line in fin:
            hap = line.strip()
            if not hap:
                continue
            match = re.match(r"^(.+)\.([12])$", hap)
            if match:
                hap_by_sample[match.group(1)].add(match.group(2))

    with open(args.out, "w") as fout:
        for sample in sorted(hap_by_sample):
            if hap_by_sample[sample] == {"1", "2"}:
                fout.write(sample + "\n")


if __name__ == "__main__":
    main()
