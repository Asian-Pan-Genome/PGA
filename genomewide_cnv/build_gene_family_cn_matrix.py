#!/usr/bin/env python3

import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate per-gene Liftoff copy-number counts into local gene-family "
            "copy numbers across haplotype assemblies."
        )
    )
    parser.add_argument("--families", required=True, help="Local gene-family TSV")
    parser.add_argument(
        "--manifest",
        required=True,
        help=(
            "Two-column TSV with header: sample_hap and gene_cn_file. Each gene CN "
            "file is a headerless two-column table: gene and CN."
        ),
    )
    parser.add_argument("--output", required=True, help="Output wide CN matrix")
    return parser.parse_args()


def load_families(path):
    families = pd.read_csv(path, sep="\t", header=0)
    required = {"genes"}
    missing = required - set(families.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")

    family_records = []
    gene_to_family = {}

    for _, row in families.iterrows():
        members = [gene for gene in str(row["genes"]).split(",") if gene]
        if not members:
            continue

        representative = members[0]
        family_name = ",".join(members)
        family_records.append((representative, family_name))

        for gene in members:
            if gene in gene_to_family and gene_to_family[gene] != representative:
                raise ValueError(
                    f"Gene {gene} occurs in multiple local families: "
                    f"{gene_to_family[gene]} and {representative}"
                )
            gene_to_family[gene] = representative

    family_df = pd.DataFrame(family_records, columns=["gene", "gene_family"])
    if family_df["gene"].duplicated().any():
        duplicated = family_df.loc[family_df["gene"].duplicated(), "gene"].tolist()
        raise ValueError(f"Duplicated representative genes: {duplicated[:10]}")

    return family_df, gene_to_family


def load_sample_cn(path, gene_to_family, family_order):
    cn = pd.read_csv(path, sep="\t", header=None, names=["gene", "CN"])
    cn["CN"] = pd.to_numeric(cn["CN"], errors="raise")

    cn["family"] = cn["gene"].map(gene_to_family)
    cn = cn.dropna(subset=["family"])

    family_cn = cn.groupby("family", sort=False)["CN"].sum()
    return family_cn.reindex(family_order, fill_value=0)


def main():
    args = parse_args()

    family_df, gene_to_family = load_families(args.families)
    family_order = family_df["gene"].tolist()

    manifest = pd.read_csv(args.manifest, sep="\t", header=0)
    required = {"sample_hap", "gene_cn_file"}
    missing = required - set(manifest.columns)
    if missing:
        raise ValueError(f"{args.manifest} missing columns: {sorted(missing)}")
    if manifest["sample_hap"].duplicated().any():
        raise ValueError("sample_hap values must be unique in the manifest")

    output = family_df.copy()
    for _, row in manifest.iterrows():
        sample_hap = str(row["sample_hap"])
        output[sample_hap] = load_sample_cn(
            str(row["gene_cn_file"]),
            gene_to_family,
            family_order,
        ).to_numpy()

    output.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
