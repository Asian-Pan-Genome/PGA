#!/usr/bin/env python3

import argparse

import numpy as np
import pandas as pd
import pysam


PARALOGS = ("PGA34A", "PGA34B", "PGA5")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a paralog-informative site manifest from the pseudo-PGA VCF."
    )
    parser.add_argument("--vcf", required=True, help="Pseudo-reference VCF.")
    parser.add_argument("--out", required=True, help="Output manifest TSV.")
    parser.add_argument(
        "--min-frequency",
        type=float,
        default=0.70,
        help="Minimum within-paralog predominant-allele frequency (default: 0.70).",
    )
    return parser.parse_args()


def normalize_specific(value):
    if value is None:
        return []
    if isinstance(value, str):
        values = [value]
    else:
        values = list(value)

    groups = []
    for value in values:
        for block in str(value).split("|"):
            groups.extend(x for x in block.split(",") if x)
    return groups


def main():
    args = parse_args()
    rows = []

    with pysam.VariantFile(args.vcf) as vcf:
        missing = [name for name in PARALOGS if name not in vcf.header.samples]
        if missing:
            raise ValueError(
                "The VCF must use the current paralog labels PGA34A, PGA34B and PGA5. "
                f"Missing sample columns: {', '.join(missing)}"
            )

        for record in vcf:
            specific = [x for x in normalize_specific(record.info.get("SPECIFIC")) if x in PARALOGS]
            if not specific:
                continue

            gfs = []
            for paralog in specific:
                gf = record.samples[paralog].get("GF")
                if gf is not None:
                    gfs.append(float(gf))
            if not gfs:
                continue

            specificity_weight = float(np.mean(gfs) / len(gfs))
            retained = [
                paralog
                for paralog in specific
                if record.samples[paralog].get("GF") is not None
                and float(record.samples[paralog]["GF"]) >= args.min_frequency
            ]
            if not retained:
                continue

            row = {
                "CHROM": record.chrom,
                "POS": record.pos,
                "REF": record.ref,
                "ALT": record.alts[0],
                "SPECIFIC_COUNT": len(retained),
                "SPECIFIC": ",".join(retained),
                "GF_SPECIFICITY": specificity_weight,
            }

            for paralog in PARALOGS:
                gt = record.samples[paralog].get("GT")
                gf = record.samples[paralog].get("GF")
                gt0 = "." if not gt or gt[0] is None else str(gt[0])
                gf_value = np.nan if gf is None else float(gf)
                row[paralog] = f"{gt0}:{gf_value}"

            rows.append(row)

    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    print(f"Manifest: {args.out} ({len(rows)} sites)")


if __name__ == "__main__":
    main()
