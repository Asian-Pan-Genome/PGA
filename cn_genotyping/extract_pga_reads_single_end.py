#!/usr/bin/env python3

import argparse
import gzip

import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract PGA-overlapping reads and depth features from a pre-aligned single-end BAM/CRAM."
    )
    parser.add_argument("--alignment", required=True, help="Input BAM or CRAM.")
    parser.add_argument("--reference", required=True, help="Whole-genome reference FASTA.")
    parser.add_argument("--pga-bed", required=True, help="BED file with PGA3, PGA4 and PGA5 in column 4.")
    parser.add_argument("--control-bed", required=True, help="BED file containing depth-control intervals.")
    parser.add_argument("--out-prefix", required=True, help="Output prefix.")
    return parser.parse_args()


def read_regions(path, named=False):
    regions = [] if not named else {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if named:
                if len(fields) < 4:
                    raise ValueError("PGA BED requires four columns: chrom, start, end, gene.")
                regions[fields[3]] = (chrom, start, end)
            else:
                regions.append((chrom, start, end))
    return regions


def overlap(start1, end1, start2, end2):
    return max(0, min(end1, end2) - max(start1, start2))


def overlap_with_regions(read, regions):
    if read.reference_end is None:
        return 0
    total = 0
    chrom = read.reference_name
    for region_chrom, start, end in regions:
        if chrom == region_chrom:
            total += overlap(read.reference_start, read.reference_end, start, end)
    return total


def main():
    args = parse_args()
    pga = read_regions(args.pga_bed, named=True)
    controls = read_regions(args.control_bed)

    required = {"PGA3", "PGA4", "PGA5"}
    missing = required - set(pga)
    if missing:
        raise ValueError(f"PGA BED is missing: {', '.join(sorted(missing))}")

    chrom = pga["PGA3"][0]
    if not all(pga[g][0] == chrom for g in required):
        raise ValueError("PGA3, PGA4 and PGA5 must occur on the same contig.")

    targets = {
        "Baseline": controls,
        "PGA": [(chrom, pga["PGA3"][1], pga["PGA5"][2])],
        "PGA34": [pga["PGA3"], pga["PGA4"]],
        "PGA5": [pga["PGA5"]],
    }
    lengths = {
        name: sum(end - start for _, start, end in regions)
        for name, regions in targets.items()
    }
    covered = {name: 0 for name in targets}

    mode = "rc" if args.alignment.lower().endswith(".cram") else "rb"
    kwargs = {"reference_filename": args.reference} if mode == "rc" else {}

    fastq_path = f"{args.out_prefix}.fq.gz"
    with pysam.AlignmentFile(args.alignment, mode, **kwargs) as alignment, \
         gzip.open(fastq_path, "wt") as fastq:
        for read in alignment.fetch(until_eof=True):
            if read.is_unmapped or read.is_supplementary or read.is_duplicate:
                continue

            for name, regions in targets.items():
                covered[name] += overlap_with_regions(read, regions)

            if overlap_with_regions(read, list(pga.values())) > 0:
                sequence = read.query_sequence
                quality = read.qual
                if not sequence or quality is None:
                    continue
                fastq.write(f"@{read.query_name}\n{sequence}\n+\n{quality}\n")

    depth = {
        name: covered[name] / lengths[name] if lengths[name] else 0.0
        for name in targets
    }

    depth_file = f"{args.out_prefix}.depth.tsv"
    with open(depth_file, "w") as handle:
        handle.write("REF_Copy\tBaseline_depth\tPGA_depth\tPGA34_depth\tPGA5_depth\n")
        handle.write(
            f"{len(pga)}\t{depth['Baseline']:.6f}\t{depth['PGA']:.6f}\t"
            f"{depth['PGA34']:.6f}\t{depth['PGA5']:.6f}\n"
        )

    print(f"Reads: {fastq_path}")
    print(f"Depth features: {depth_file}")


if __name__ == "__main__":
    main()
