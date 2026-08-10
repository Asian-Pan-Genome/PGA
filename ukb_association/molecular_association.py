#!/usr/bin/env python3

import os

# Avoid BLAS/OpenMP oversubscription when multiprocessing is used.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import argparse
import multiprocessing as mp
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm, rankdata
from statsmodels.stats.multitest import multipletests


COVARIATES = ["Age", "Sex"] + [f"PC{i}" for i in range(1, 11)]
_WORKER_DF = None
_TARGET_CN = None
_MIN_N = None


def parse_args():
    parser = argparse.ArgumentParser(
        description="Test UK Biobank molecular traits for association with PGA34A copy number."
    )
    parser.add_argument("--cohort", required=True, help="Prepared UKB association cohort.")
    parser.add_argument("--proteomics", required=True, help="UKB proteomics table.")
    parser.add_argument(
        "--molecular-fields",
        required=True,
        help="One baseline UKB molecular/biochemical field per line.",
    )
    parser.add_argument("--cn-column", default="Pred_PGA34A")
    parser.add_argument("--min-n", type=int, default=100)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def read_field_list(path):
    fields = [
        line.strip()
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if not fields:
        raise ValueError(f"No fields found in {path}")
    return fields


def rint(series):
    ranks = rankdata(series, method="average")
    quantiles = (ranks - 0.5) / len(ranks)
    return norm.ppf(quantiles)


def init_worker(df, target_cn, min_n):
    global _WORKER_DF, _TARGET_CN, _MIN_N
    _WORKER_DF = df
    _TARGET_CN = target_cn
    _MIN_N = min_n


def run_association(trait):
    temp = _WORKER_DF[[_TARGET_CN, *COVARIATES, trait]].dropna()
    if len(temp) < _MIN_N:
        return None

    y = rint(temp[trait].to_numpy())
    x = temp[[_TARGET_CN, *COVARIATES]].astype(float)
    x = sm.add_constant(x, has_constant="add")
    model = sm.OLS(y, x).fit()

    return {
        "Biomarker": trait,
        "Beta": model.params[_TARGET_CN],
        "SE": model.bse[_TARGET_CN],
        "P_value": model.pvalues[_TARGET_CN],
        "N_samples": len(temp),
    }


def main():
    args = parse_args()
    molecular_fields = read_field_list(args.molecular_fields)

    cohort = pd.read_csv(args.cohort, sep="\t")
    required = {"eid", args.cn_column, *COVARIATES, *molecular_fields}
    missing = sorted(required - set(cohort.columns))
    if missing:
        raise ValueError("Cohort table is missing: " + ", ".join(missing))

    proteomics = pd.read_csv(args.proteomics, sep="\t")
    if "eid" not in proteomics.columns:
        raise ValueError("Proteomics table must contain 'eid'.")

    protein_fields = [col for col in proteomics.columns if col != "eid"]
    overlap = sorted(set(protein_fields) & set(molecular_fields))
    if overlap:
        raise ValueError("Traits occur in both inputs: " + ", ".join(overlap))

    df = pd.merge(cohort, proteomics, on="eid", how="left")
    traits = protein_fields + molecular_fields

    n_threads = max(1, args.threads)
    if n_threads == 1:
        init_worker(df, args.cn_column, args.min_n)
        results = [run_association(trait) for trait in traits]
    else:
        with mp.Pool(
            processes=n_threads,
            initializer=init_worker,
            initargs=(df, args.cn_column, args.min_n),
        ) as pool:
            results = pool.map(run_association, traits)

    results = pd.DataFrame([x for x in results if x is not None])
    if results.empty:
        raise RuntimeError("No traits passed the minimum sample-size filter.")

    results["Bonferroni_P"] = multipletests(
        results["P_value"],
        method="bonferroni",
    )[1]
    results["Significant"] = results["Bonferroni_P"] < 0.05
    results.sort_values("Bonferroni_P", inplace=True)
    results.to_csv(args.out, sep="\t", index=False)

    n_protein = results["Biomarker"].isin(protein_fields).sum()
    n_other = results["Biomarker"].isin(molecular_fields).sum()
    print(f"Tested traits: {len(results)}")
    print(f"Proteins: {n_protein}")
    print(f"Other molecular/biochemical traits: {n_other}")


if __name__ == "__main__":
    main()
