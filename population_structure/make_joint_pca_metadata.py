#!/usr/bin/env python3

import argparse
import pandas as pd


POP_TO_NEW_SUPERPOP = {
    "GWD": "AFR-W", "MSL": "AFR-W", "ESN": "AFR-W", "GWW": "AFR-W",
    "GWJ": "AFR-W", "YRI": "AFR-W", "GWF": "AFR-W",
    "LWK": "AFR-E&S", "ASW": "AFR-Other", "ACB": "AFR-Other",
    "CLM": "AMR", "PEL": "AMR", "MXL": "AMR", "PUR": "AMR",
    "CHS": "EAS", "KHV": "EAS", "JPT": "EAS", "CHB": "EAS", "CDX": "EAS",
    "IBS": "EUR", "TSI": "EUR", "FIN": "EUR", "GBR": "EUR", "CEU": "EUR",
    "PJL": "CSA", "BEB": "SAS", "STU": "SAS", "ITU": "SAS", "GIH": "SAS",
}


def read_samples(path):
    with open(path) as fin:
        return [line.strip() for line in fin if line.strip()]


def require_columns(df, columns, label):
    missing = set(columns) - set(df.columns)
    if missing:
        raise ValueError(f"{label} missing columns: {sorted(missing)}")


def main():
    parser = argparse.ArgumentParser(description="Build metadata for joint assembly-1KGP PCA.")
    parser.add_argument("--joint-vcf-samples", required=True)
    parser.add_argument("--assembly-copies", required=True)
    parser.add_argument("--assembly-new-superpop", required=True)
    parser.add_argument("--kg-metadata", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    joint_samples = set(read_samples(args.joint_vcf_samples))

    copies = pd.read_csv(args.assembly_copies, sep="\t", dtype=str)
    new = pd.read_csv(args.assembly_new_superpop, sep="\t", dtype=str)
    require_columns(copies, ["Sample", "Superpop"], "assembly copy-number table")
    require_columns(new, ["Sample", "New_superpop"], "assembly population table")

    copies["Sample"] = copies["Sample"].replace({"HG002v1.1": "HG002"})
    new["Sample"] = new["Sample"].replace({"HG002v1.1": "HG002"})

    asm = copies[["Sample", "Superpop"]].drop_duplicates()
    asm = asm.merge(new[["Sample", "New_superpop"]].drop_duplicates(), on="Sample", how="left")
    asm["New_superpop"] = asm["New_superpop"].replace({"NA": pd.NA})
    asm["New_superpop"] = asm["New_superpop"].fillna("AFR-Other")

    asm_out = pd.DataFrame({
        "IID": "ASM_" + asm["Sample"].astype(str),
        "Basename": asm["Sample"],
        "Panel": "Assembly panel",
        "Population": asm["Superpop"],
        "Superpopulation": asm["New_superpop"],
    })

    kg = pd.read_csv(args.kg_metadata, sep=r"\s+", dtype=str)
    require_columns(kg, ["SampleID", "Population"], "1KGP metadata")
    kg["New_superpop"] = kg["Population"].map(POP_TO_NEW_SUPERPOP)

    unmapped = sorted(kg.loc[kg["New_superpop"].isna(), "Population"].dropna().unique())
    if unmapped:
        raise ValueError("Unmapped 1KGP population codes: " + ",".join(unmapped))

    kg_out = pd.DataFrame({
        "IID": "KG_" + kg["SampleID"].astype(str),
        "Basename": kg["SampleID"],
        "Panel": "1KGP panel",
        "Population": kg["Population"],
        "Superpopulation": kg["New_superpop"],
    })

    meta = pd.concat([asm_out, kg_out], ignore_index=True)
    meta = meta[meta["IID"].isin(joint_samples)].copy()
    meta.to_csv(args.out, sep="\t", index=False)

    print(f"metadata samples: {meta.shape[0]}")
    print(meta.groupby(["Panel", "Superpopulation"], dropna=False).size())


if __name__ == "__main__":
    main()
