#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


PC_FIELDS = [f"22009-0.{i}" for i in range(1, 11)]
QC_FIELDS = [
    "31-0.0",      # recorded sex
    "22001-0.0",   # genetic sex
    "30079-0.0",   # genomic ancestry
    "22006-0.0",   # genetic ethnic grouping
    "22019-0.0",   # sex chromosome aneuploidy
    "22020-0.0",   # used in genetic PCA
    "22021-0.0",   # genetic kinship
    "22027-0.0",   # heterozygosity/missingness outlier
]
AGE_FIELD = "21022-0.0"
ICD10_FIELD = "41270-0.0"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare the UK Biobank cohort for PGA34A molecular association and PheWAS."
    )
    parser.add_argument("--cn", required=True, help="UKB PGA CN prediction table.")
    parser.add_argument("--phenotypes", required=True, help="UKB phenotype table.")
    parser.add_argument(
        "--molecular-fields",
        required=True,
        help="One baseline UKB molecular/biochemical field per line.",
    )
    parser.add_argument("--out-prefix", required=True)
    return parser.parse_args()


def read_field_list(path):
    fields = [
        line.strip()
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if not fields:
        raise ValueError(f"No fields found in {path}")
    if len(fields) != len(set(fields)):
        raise ValueError(f"Duplicate fields found in {path}")
    return fields


def load_cn(path):
    cn = pd.read_csv(path, sep="\t")

    if "eid" not in cn.columns:
        if "Sample" not in cn.columns:
            raise ValueError("CN table must contain either 'eid' or 'Sample'.")
        cn["eid"] = cn["Sample"].astype(str).str.split("_").str[0].astype(int)

    if "Pred_PGA34A" not in cn.columns:
        raise ValueError("CN table is missing 'Pred_PGA34A'.")

    return cn[["eid", "Pred_PGA34A"]].copy()


def main():
    args = parse_args()
    molecular_fields = read_field_list(args.molecular_fields)

    required_fields = [
        "eid",
        AGE_FIELD,
        ICD10_FIELD,
        *QC_FIELDS,
        *PC_FIELDS,
        *molecular_fields,
    ]

    phenotype_header = pd.read_csv(args.phenotypes, sep="\t", nrows=0).columns
    missing = sorted(set(required_fields) - set(phenotype_header))
    if missing:
        raise ValueError(
            "Phenotype table is missing required fields: " + ", ".join(missing)
        )

    phenotypes = pd.read_csv(
        args.phenotypes,
        sep="\t",
        usecols=required_fields,
    )

    phenotypes = phenotypes[
        (phenotypes["31-0.0"] == phenotypes["22001-0.0"])
        & (phenotypes["30079-0.0"] == 5)
        & (phenotypes["22006-0.0"] == 1)
        & phenotypes["22019-0.0"].isna()
        & (phenotypes["22020-0.0"] == 1)
        & (phenotypes["22021-0.0"] != 10)
        & phenotypes["22027-0.0"].isna()
    ].copy()

    rename = {
        AGE_FIELD: "Age",
        "31-0.0": "Sex",
        **{f"22009-0.{i}": f"PC{i}" for i in range(1, 11)},
    }
    phenotypes.rename(columns=rename, inplace=True)

    cn = load_cn(args.cn)
    cohort = pd.merge(cn, phenotypes, on="eid", how="inner")

    covariates = ["Age", "Sex"] + [f"PC{i}" for i in range(1, 11)]
    cohort = cohort.dropna(subset=["Pred_PGA34A", *covariates]).copy()

    cohort_columns = [
        "eid",
        "Pred_PGA34A",
        *covariates,
        *molecular_fields,
    ]
    cohort[cohort_columns].to_csv(
        f"{args.out_prefix}.cohort.tsv",
        sep="\t",
        index=False,
    )

    icd = cohort[["eid", ICD10_FIELD]].dropna().copy()
    icd[ICD10_FIELD] = icd[ICD10_FIELD].astype(str).str.split("|")
    icd = icd.explode(ICD10_FIELD)
    icd.rename(columns={"eid": "id", ICD10_FIELD: "code"}, inplace=True)
    icd["code"] = icd["code"].astype(str).str.strip()
    icd = icd[icd["code"] != ""].copy()
    icd["vocabulary_id"] = "ICD10"
    icd["count"] = 1
    icd[["id", "vocabulary_id", "code", "count"]].to_csv(
        f"{args.out_prefix}.icd10.tsv",
        sep="\t",
        index=False,
    )

    print(f"Association cohort: {len(cohort)} participants")
    print(f"Participants with >=1 ICD-10 code: {icd['id'].nunique()}")
    print(f"ICD-10 entries: {len(icd)}")


if __name__ == "__main__":
    main()
