#!/usr/bin/env python3

import argparse
from collections import defaultdict

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Convert OrthoFinder orthogroups into chromosome-local gene families "
            "using GENCODE gene coordinates."
        )
    )
    parser.add_argument("--orthogroups", required=True, help="OrthoFinder Orthogroups.tsv")
    parser.add_argument(
        "--unassigned",
        required=True,
        help="OrthoFinder Orthogroups_UnassignedGenes.tsv",
    )
    parser.add_argument("--gff3", required=True, help="GENCODE GFF3 used for coordinates")
    parser.add_argument(
        "--reference-column",
        required=True,
        help="Reference-proteome column in the OrthoFinder tables",
    )
    parser.add_argument(
        "--duplicate-column",
        required=True,
        help=(
            "Technical-duplicate reference-proteome column. Orthogroup membership "
            "must match --reference-column exactly."
        ),
    )
    parser.add_argument(
        "--max-extension",
        type=int,
        default=5_000_000,
        help="Maximum stepwise extension within a local family in bp (default: 5000000)",
    )
    parser.add_argument("--output", required=True, help="Output TSV")
    return parser.parse_args()


def parse_attributes(text):
    attributes = {}
    for field in text.split(";"):
        if "=" not in field:
            continue
        key, value = field.split("=", 1)
        attributes[key] = value
    return attributes


def load_gene_coordinates(gff3):
    coordinates = {}
    with open(gff3) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attributes = parse_attributes(fields[8])
            gene_name = attributes.get("gene_name")
            if gene_name is None:
                continue
            coordinates[gene_name] = (fields[0], int(fields[3]), int(fields[4]))
    return coordinates


def parse_gene_set(value):
    if pd.isna(value):
        return set()
    return {gene.strip() for gene in str(value).split(",") if gene.strip()}


def add_local_families(orthogroup, genes, coordinates, max_extension, records):
    genes = [gene for gene in genes if not gene.startswith("ENSG")]
    if not genes:
        return

    missing = sorted(gene for gene in genes if gene not in coordinates)
    if missing:
        raise KeyError(
            f"Missing GFF3 coordinates for {orthogroup}: {', '.join(missing)}"
        )

    by_chrom = defaultdict(list)
    for gene in genes:
        chrom, start, end = coordinates[gene]
        by_chrom[chrom].append((gene, start, end))

    for chrom, chrom_genes in by_chrom.items():
        chrom_genes.sort(key=lambda x: x[1])

        first_gene, start, end = chrom_genes[0]
        current_genes = [first_gene]

        for gene, gene_start, gene_end in chrom_genes[1:]:
            if gene_end <= end:
                current_genes.append(gene)
                continue

            if gene_end - end < max_extension:
                current_genes.append(gene)
                end = gene_end
                continue

            records.append((orthogroup, chrom, start, end, current_genes.copy()))
            start, end = gene_start, gene_end
            current_genes = [gene]

        records.append((orthogroup, chrom, start, end, current_genes.copy()))


def main():
    args = parse_args()
    coordinates = load_gene_coordinates(args.gff3)

    orthogroups = pd.read_csv(args.orthogroups, sep="\t", header=0, index_col=0)
    for column in (args.reference_column, args.duplicate_column):
        if column not in orthogroups.columns:
            raise KeyError(f"Column not found in {args.orthogroups}: {column}")

    records = []
    for orthogroup, row in orthogroups.iterrows():
        reference_genes = parse_gene_set(row[args.reference_column])
        duplicate_genes = parse_gene_set(row[args.duplicate_column])
        if reference_genes != duplicate_genes:
            raise ValueError(
                f"Technical-duplicate columns differ for {orthogroup}: "
                f"{sorted(reference_genes)} vs {sorted(duplicate_genes)}"
            )
        add_local_families(
            orthogroup,
            sorted(reference_genes),
            coordinates,
            args.max_extension,
            records,
        )

    unassigned = pd.read_csv(args.unassigned, sep="\t", header=0, index_col=0)
    if args.duplicate_column not in unassigned.columns:
        raise KeyError(f"Column not found in {args.unassigned}: {args.duplicate_column}")

    for orthogroup, row in unassigned.iterrows():
        genes = parse_gene_set(row[args.duplicate_column])
        add_local_families(
            orthogroup,
            sorted(genes),
            coordinates,
            args.max_extension,
            records,
        )

    output = pd.DataFrame(
        records,
        columns=["rowindex", "chrom", "start", "end", "genes"],
    )
    output["genes_number"] = output["genes"].str.len()
    output["genes"] = output["genes"].apply(",".join)
    output = output.sort_values(
        ["genes_number", "chrom", "start"],
        ascending=[False, True, True],
    )
    output = output.drop_duplicates(
        subset=["chrom", "start", "end", "genes"],
        keep="first",
    )
    output.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
