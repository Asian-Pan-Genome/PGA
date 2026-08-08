#!/usr/bin/env python3

import argparse
import gzip
import sys
from collections import defaultdict


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def gt_dosage(sample_field):
    """Return diploid ALT dosage for a biallelic GT, ignoring phase."""
    gt = sample_field.split(":", 1)[0]
    if gt in {".", "./.", ".|."}:
        return None

    if "|" in gt:
        alleles = gt.split("|")
    elif "/" in gt:
        alleles = gt.split("/")
    else:
        return None

    if len(alleles) != 2 or any(a not in {"0", "1"} for a in alleles):
        return None
    return int(alleles[0]) + int(alleles[1])


def panel_missingness(genotypes, indices):
    called = sum(gt_dosage(genotypes[i]) is not None for i in indices)
    return 1.0 - called / len(indices)


def main():
    parser = argparse.ArgumentParser(
        description="Filter panel-biased SNPs using overlapping ASM_/KG_ individuals."
    )
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--out-vcf", required=True, help="Output uncompressed VCF")
    parser.add_argument("--min-overlap-called", type=int, default=10)
    parser.add_argument("--min-concordance", type=float, default=0.98)
    parser.add_argument("--max-missing-diff", type=float, default=0.05)
    parser.add_argument("--max-panel-missing", type=float, default=0.20)
    args = parser.parse_args()

    if args.out_vcf.endswith(".gz"):
        raise ValueError("--out-vcf must be an uncompressed .vcf")

    qc_path = args.out_prefix + ".site_qc.tsv"
    kept_path = args.out_prefix + ".kept.keys.tsv"
    rejected_path = args.out_prefix + ".rejected.keys.tsv"

    n_total = 0
    n_keep = 0
    reject_counter = defaultdict(int)

    with open_text(args.vcf) as fin, \
            open(args.out_vcf, "w") as fout, \
            open(qc_path, "w") as fqc, \
            open(kept_path, "w") as fkeep, \
            open(rejected_path, "w") as freject:

        asm_indices = []
        kg_indices = []
        overlap_pairs = []

        fqc.write(
            "CHROM\tPOS\tREF\tALT\tn_overlap_called\tn_overlap_concordant\t"
            "overlap_concordance\tasm_missing\tkg_missing\tmissing_diff\tkeep\tfail_reason\n"
        )

        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue

            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]

                asm_by_base = {}
                kg_by_base = {}
                for i, sample in enumerate(samples):
                    if sample.startswith("ASM_"):
                        base = sample[4:]
                        if base in asm_by_base:
                            raise ValueError(f"Duplicate ASM sample: {base}")
                        asm_by_base[base] = i
                        asm_indices.append(i)
                    elif sample.startswith("KG_"):
                        base = sample[3:]
                        if base in kg_by_base:
                            raise ValueError(f"Duplicate KG sample: {base}")
                        kg_by_base[base] = i
                        kg_indices.append(i)

                overlap = sorted(set(asm_by_base) & set(kg_by_base))
                overlap_pairs = [(asm_by_base[s], kg_by_base[s]) for s in overlap]

                if not asm_indices or not kg_indices or not overlap_pairs:
                    raise ValueError("Joint VCF must contain ASM_, KG_, and overlapping sample pairs")

                print(f"ASM samples: {len(asm_indices)}", file=sys.stderr)
                print(f"KG samples: {len(kg_indices)}", file=sys.stderr)
                print(f"Overlapping sample pairs: {len(overlap_pairs)}", file=sys.stderr)
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            genotypes = fields[9:]
            n_total += 1
            fail = []

            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                fail.append("not_biallelic_snp")

            n_called = 0
            n_concordant = 0
            for asm_i, kg_i in overlap_pairs:
                asm_gt = gt_dosage(genotypes[asm_i])
                kg_gt = gt_dosage(genotypes[kg_i])
                if asm_gt is None or kg_gt is None:
                    continue
                n_called += 1
                n_concordant += asm_gt == kg_gt

            concordance = n_concordant / n_called if n_called else float("nan")
            asm_missing = panel_missingness(genotypes, asm_indices)
            kg_missing = panel_missingness(genotypes, kg_indices)
            missing_diff = abs(asm_missing - kg_missing)

            if n_called < args.min_overlap_called:
                fail.append("low_overlap_called")
            if concordance != concordance or concordance < args.min_concordance:
                fail.append("low_concordance")
            if asm_missing > args.max_panel_missing:
                fail.append("high_asm_missing")
            if kg_missing > args.max_panel_missing:
                fail.append("high_kg_missing")
            if missing_diff > args.max_missing_diff:
                fail.append("high_panel_missing_diff")

            keep = not fail
            reason = "PASS" if keep else ",".join(fail)
            concordance_text = "NA" if concordance != concordance else f"{concordance:.6g}"
            fqc.write(
                f"{chrom}\t{pos}\t{ref}\t{alt}\t{n_called}\t{n_concordant}\t{concordance_text}\t"
                f"{asm_missing:.6g}\t{kg_missing:.6g}\t{missing_diff:.6g}\t"
                f"{int(keep)}\t{reason}\n"
            )

            key = f"{chrom}\t{pos}\t{ref}\t{alt}\n"
            if keep:
                fout.write(line)
                fkeep.write(key)
                n_keep += 1
            else:
                freject.write(key)
                for item in fail:
                    reject_counter[item] += 1

    print("Filtering finished.", file=sys.stderr)
    print(f"Total variants: {n_total}", file=sys.stderr)
    print(f"Kept variants: {n_keep}", file=sys.stderr)
    print(f"Rejected variants: {n_total - n_keep}", file=sys.stderr)
    for reason, count in sorted(reject_counter.items()):
        print(f"Rejected [{reason}]: {count}", file=sys.stderr)


if __name__ == "__main__":
    main()
