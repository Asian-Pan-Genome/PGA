#!/usr/bin/env python3

import argparse

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a gap-free pseudo-PGA reference from an aligned PGA34A template and a variable-site VCF."
    )
    parser.add_argument("--template", required=True, help="Aligned FASTA containing the PGA34A template sequence.")
    parser.add_argument("--vcf", required=True, help="Paralog-resolved variable-site VCF.")
    parser.add_argument("--out-prefix", required=True, help="Output prefix.")
    parser.add_argument(
        "--template-id",
        default=None,
        help="Template FASTA record ID. If omitted, the CHROM value of the first VCF record is used.",
    )
    parser.add_argument(
        "--reference-name",
        default=None,
        help="Name of the pseudo-reference contig. Defaults to the output-prefix basename.",
    )
    return parser.parse_args()


def load_template(path, template_id):
    records = {record.id: record for record in SeqIO.parse(path, "fasta")}
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    if template_id not in records:
        raise ValueError(f"Template ID '{template_id}' was not found in {path}")
    return list(str(records[template_id].seq))


def flip_biallelic_gt(gt):
    if gt is None:
        return gt
    flipped = []
    for allele in gt:
        if allele is None:
            flipped.append(None)
        elif allele == 0:
            flipped.append(1)
        elif allele == 1:
            flipped.append(0)
        else:
            flipped.append(allele)
    return tuple(flipped)


def main():
    args = parse_args()

    with pysam.VariantFile(args.vcf) as vcf_in:
        first_record = next(iter(vcf_in), None)
        if first_record is None:
            raise ValueError(f"No variant records found in {args.vcf}")
        inferred_template_id = first_record.chrom

    template_id = args.template_id or inferred_template_id
    template = load_template(args.template, template_id)
    reference_name = args.reference_name or args.out_prefix.rsplit("/", 1)[-1]

    pseudo_vcf = f"{args.out_prefix}.pseudo.vcf"
    pseudo_fasta = f"{args.out_prefix}.pseudo.fa"

    with pysam.VariantFile(args.vcf) as vcf_in:
        header = vcf_in.header.copy()
        if reference_name not in header.contigs:
            header.contigs.add(reference_name, length=len(template))

        with pysam.VariantFile(pseudo_vcf, "w", header=header) as vcf_out:
            for record in vcf_in:
                record.translate(header)

                if record.pos < 1 or record.pos > len(template):
                    raise ValueError(
                        f"VCF position {record.pos} is outside the template length ({len(template)})."
                    )

                if record.ref == "-":
                    if not record.alts or len(record.alts) != 1:
                        raise ValueError(f"Expected one ALT allele at {record.chrom}:{record.pos}")
                    replacement = record.alts[0]
                    if replacement == "-":
                        raise ValueError(f"Cannot replace a template gap with '-' at {record.chrom}:{record.pos}")

                    template[record.pos - 1] = replacement
                    record.ref = replacement
                    record.alts = ("-",)
                    for sample in record.samples:
                        record.samples[sample]["GT"] = flip_biallelic_gt(record.samples[sample].get("GT"))

                record.chrom = reference_name
                vcf_out.write(record)

    pseudo_sequence = "".join(template)
    if "-" in pseudo_sequence:
        raise ValueError(
            "The pseudo-reference still contains alignment gaps. Check that the variable-site VCF covers all template-gap positions."
        )

    SeqIO.write(
        [SeqRecord(Seq(pseudo_sequence), id=reference_name, description="")],
        pseudo_fasta,
        "fasta",
    )

    print(f"Pseudo-reference FASTA: {pseudo_fasta}")
    print(f"Coordinate-matched VCF: {pseudo_vcf}")


if __name__ == "__main__":
    main()
