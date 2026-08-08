#!/usr/bin/env python3

import argparse
import gzip
import os
import re
import sys


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_samples(path):
    with open(path) as fin:
        return [line.strip() for line in fin if line.strip()]


def read_failed_samples(path):
    """Collapse failed Sample.1/Sample.2 haplotypes to diploid sample IDs."""
    failed = set()
    if not path or not os.path.exists(path):
        return failed

    with open(path) as fin:
        for line in fin:
            value = line.strip()
            if not value:
                continue
            match = re.match(r"^(.+)\.([12])$", value)
            failed.add(match.group(1) if match else value)
    return failed


def normalize_gt(sample_field):
    """Return phased biallelic diploid GT or .|. for missing/invalid GT."""
    if sample_field is None:
        return ".|."

    gt = sample_field.split(":", 1)[0]
    if gt in {".", "./.", ".|."} or "|" not in gt:
        return ".|."

    alleles = gt.split("|")
    if len(alleles) != 2 or any(a not in {"0", "1"} for a in alleles):
        return ".|."
    return "|".join(alleles)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Standardize one Minigraph-Cactus regional VCF to a fixed diploid "
            "assembly sample set and mask region-specific failed samples."
        )
    )
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--samples", required=True, help="Master diploid sample list")
    parser.add_argument("--failed-list", required=True, help="Region-specific failed.list")
    parser.add_argument("--expected-chrom", required=True, help="Expected contig, e.g. region1")
    parser.add_argument("--out", required=True, help="Output uncompressed VCF")
    args = parser.parse_args()

    master_samples = read_samples(args.samples)
    failed_samples = read_failed_samples(args.failed_list)

    n_variants = 0
    n_wrong_chrom = 0
    n_masked = 0
    saw_header = False

    with open_text(args.vcf) as fin, open(args.out, "w") as fout:
        sample_to_idx = {}

        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue

            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                sample_to_idx = {s: i for i, s in enumerate(fields[9:])}
                fout.write("\t".join(fields[:9] + master_samples) + "\n")
                saw_header = True
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            if fields[0] != args.expected_chrom:
                n_wrong_chrom += 1
                continue

            fixed = fields[:9]
            fixed[8] = "GT"
            input_gts = fields[9:]
            out_gts = []

            for sample in master_samples:
                if sample in failed_samples or sample not in sample_to_idx:
                    out_gts.append(".|.")
                    n_masked += 1
                    continue

                gt = normalize_gt(input_gts[sample_to_idx[sample]])
                out_gts.append(gt)
                if gt == ".|.":
                    n_masked += 1

            fout.write("\t".join(fixed + out_gts) + "\n")
            n_variants += 1

    if not saw_header:
        raise ValueError(f"VCF header not found: {args.vcf}")

    print(f"input_vcf={args.vcf}", file=sys.stderr)
    print(f"expected_chrom={args.expected_chrom}", file=sys.stderr)
    print(f"n_master_samples={len(master_samples)}", file=sys.stderr)
    print(f"n_region_failed_diploid_samples={len(failed_samples)}", file=sys.stderr)
    print(f"n_variants_written={n_variants}", file=sys.stderr)
    print(f"n_wrong_chrom_skipped={n_wrong_chrom}", file=sys.stderr)
    print(f"n_masked_gt_fields={n_masked}", file=sys.stderr)


if __name__ == "__main__":
    main()
