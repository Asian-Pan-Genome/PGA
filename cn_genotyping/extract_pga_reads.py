#!/usr/bin/env python3

import argparse
import os
import subprocess
import tempfile

import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract PGA-overlapping read pairs and calculate depth features from a BAM/CRAM alignment."
    )
    parser.add_argument("--alignment", required=True, help="Input BAM or CRAM.")
    parser.add_argument("--reference", required=True, help="Whole-genome reference FASTA used by the alignment.")
    parser.add_argument("--pga-bed", required=True, help="BED file with PGA3, PGA4 and PGA5 in column 4.")
    parser.add_argument("--control-bed", required=True, help="BED file containing depth-control intervals.")
    parser.add_argument("--out-prefix", required=True, help="Output prefix.")
    parser.add_argument("--threads", type=int, default=16)
    return parser.parse_args()


def read_pga_bed(path):
    regions = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 4:
                raise ValueError("PGA BED requires four columns: chrom, start, end, gene.")
            chrom, start, end, gene = fields[:4]
            regions[gene] = (chrom, int(start), int(end))

    required = {"PGA3", "PGA4", "PGA5"}
    missing = required - set(regions)
    if missing:
        raise ValueError(f"PGA BED is missing: {', '.join(sorted(missing))}")
    return regions


def read_bed(path):
    regions = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end = line.rstrip().split("\t")[:3]
            regions.append((chrom, int(start), int(end)))
    if not regions:
        raise ValueError(f"No intervals found in {path}")
    return regions


def mean_pileup_depth(alignment, reference, regions):
    depth_sum = 0
    length_sum = 0
    for chrom, start, end in regions:
        for column in alignment.pileup(
            chrom,
            start,
            end,
            truncate=True,
            stepper="samtools",
            fastafile=reference,
            min_base_quality=0,
            ignore_overlaps=True,
            ignore_orphans=False,
        ):
            depth_sum += column.nsegments
        length_sum += end - start
    return depth_sum / length_sum if length_sum else 0.0


def alignment_mode(path):
    return "rc" if path.lower().endswith(".cram") else "rb"


def main():
    args = parse_args()
    pga = read_pga_bed(args.pga_bed)
    controls = read_bed(args.control_bed)

    mode = alignment_mode(args.alignment)
    open_kwargs = {"reference_filename": args.reference} if mode == "rc" else {}

    with tempfile.TemporaryDirectory(prefix="pga_reads_") as tmpdir:
        candidate_bam = os.path.join(tmpdir, "candidate.bam")
        name_sorted_bam = os.path.join(tmpdir, "candidate.name_sorted.bam")

        with pysam.AlignmentFile(args.alignment, mode, **open_kwargs) as alignment, \
             pysam.FastaFile(args.reference) as reference, \
             pysam.AlignmentFile(candidate_bam, "wb", template=alignment) as out_bam:

            seen = set()
            for chrom, start, end in pga.values():
                for read in alignment.fetch(chrom, start, end):
                    if read.is_unmapped or read.is_supplementary or read.is_duplicate:
                        continue
                    key = (read.query_name, read.flag, read.reference_id, read.reference_start)
                    if key not in seen:
                        out_bam.write(read)
                        seen.add(key)

            baseline_depth = mean_pileup_depth(alignment, reference, controls)

            pga_chrom = pga["PGA3"][0]
            if not all(pga[g][0] == pga_chrom for g in ("PGA3", "PGA4", "PGA5")):
                raise ValueError("PGA3, PGA4 and PGA5 must occur on the same contig.")

            pga_depth = mean_pileup_depth(
                alignment,
                reference,
                [(pga_chrom, pga["PGA3"][1], pga["PGA5"][2])],
            )
            pga34_depth = mean_pileup_depth(
                alignment,
                reference,
                [pga["PGA3"], pga["PGA4"]],
            )
            pga5_depth = mean_pileup_depth(alignment, reference, [pga["PGA5"]])

        subprocess.run(
            [
                "samtools", "sort", "-n", "-@", str(args.threads), "-m", "5G",
                "-o", name_sorted_bam, candidate_bam,
            ],
            check=True,
        )
        subprocess.run(
            [
                "samtools", "fastq", "-@", str(args.threads), "-n",
                "-1", f"{args.out_prefix}.R1.fq.gz",
                "-2", f"{args.out_prefix}.R2.fq.gz",
                "-0", os.devnull,
                "-s", os.devnull,
                name_sorted_bam,
            ],
            check=True,
        )

    depth_file = f"{args.out_prefix}.depth.tsv"
    with open(depth_file, "w") as handle:
        handle.write("REF_Copy\tBaseline_depth\tPGA_depth\tPGA34_depth\tPGA5_depth\n")
        handle.write(
            f"{len(pga)}\t{baseline_depth:.6f}\t{pga_depth:.6f}\t{pga34_depth:.6f}\t{pga5_depth:.6f}\n"
        )

    print(f"Read pairs: {args.out_prefix}.R1.fq.gz, {args.out_prefix}.R2.fq.gz")
    print(f"Depth features: {depth_file}")


if __name__ == "__main__":
    main()
