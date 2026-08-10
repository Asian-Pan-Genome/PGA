#!/usr/bin/env python3
"""Convert a three-haplotype NAHR alignment to a focal-referenced variant table."""

from __future__ import annotations

import argparse
from collections import OrderedDict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", required=True, type=Path, help="Aligned FASTA from MAFFT.")
    parser.add_argument("--focal", required=True, help="Focal candidate haplotype ID.")
    parser.add_argument("--left", required=True, help="Left-context comparator haplotype ID.")
    parser.add_argument("--right", required=True, help="Right-context comparator haplotype ID.")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV.")
    return parser.parse_args()


def read_fasta(path: Path) -> OrderedDict[str, str]:
    records: OrderedDict[str, str] = OrderedDict()
    current = None
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in records:
                    raise ValueError(f"Duplicate FASTA ID: {current}")
                records[current] = ""
            else:
                if current is None:
                    raise ValueError("Sequence found before FASTA header")
                records[current] += line.upper()

    lengths = {len(seq) for seq in records.values()}
    if len(lengths) != 1:
        raise ValueError("Alignment sequences have unequal lengths")
    return records


def main() -> None:
    args = parse_args()
    alignment = read_fasta(args.alignment)

    required = [args.left, args.focal, args.right]
    missing = [name for name in required if name not in alignment]
    if missing:
        raise ValueError(f"Missing sequence(s) from alignment: {missing}")

    left = alignment[args.left]
    focal = alignment[args.focal]
    right = alignment[args.right]

    focal_pos = 0
    rows = []
    for i, (lbase, fbase, rbase) in enumerate(zip(left, focal, right), start=1):
        if fbase == "-":
            continue
        focal_pos += 1

        if len({lbase, fbase, rbase}) == 1:
            continue

        if lbase != rbase and fbase == lbase and fbase != rbase:
            state = "left"
        elif lbase != rbase and fbase == rbase and fbase != lbase:
            state = "right"
        elif fbase == lbase == rbase:
            state = "both"
        else:
            state = "other"

        rows.append((args.focal, focal_pos, i, fbase, lbase, rbase, state))

    with args.output.open("w") as out:
        out.write("focal_haplotype\tposition\talignment_position\tfocal_allele\tleft_allele\tright_allele\tfocal_matches\n")
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")

    print(f"Variable sites: {len(rows)}")
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
