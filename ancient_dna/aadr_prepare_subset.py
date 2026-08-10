#!/usr/bin/env python3
"""Prepare a chromosome-wide ancient AADR subset for the PGA SNP-panel analysis.

The script retains ancient individuals, keeps one Genetic ID per Individual ID
(the record with the largest number of SNPs hit on the AADR 2M autosomal
targets), preserves the original AADR .ind sample order, writes matched metadata,
and generates an EIGENSOFT convertf parameter file for a single chromosome.

The downstream PGA analysis uses chromosome 11-wide genotypes because its
empirical background is sampled across chromosome 11.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def norm_col(value: str) -> str:
    return " ".join(str(value).replace("\xa0", " ").split()).lower()


def find_col(columns, *, startswith=None, contains_all=None, exact=None):
    matches = []
    for column in columns:
        normalized = norm_col(column)
        if exact is not None and normalized != exact.lower():
            continue
        if startswith is not None and not normalized.startswith(startswith.lower()):
            continue
        if contains_all is not None and not all(term.lower() in normalized for term in contains_all):
            continue
        matches.append(column)

    if len(matches) != 1:
        raise KeyError(f"Expected one matching column, found {len(matches)}: {matches}")
    return matches[0]


def to_numeric(series):
    return pd.to_numeric(series.replace({"..": np.nan, "": np.nan}), errors="coerce")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare a chromosome-wide ancient AADR subset for EIGENSTRAT analysis."
    )
    parser.add_argument("--prefix", required=True, help="AADR prefix, e.g. v66.2M.aadr.PUB")
    parser.add_argument("--anno", required=True, help="AADR annotation table (.anno)")
    parser.add_argument("--out-prefix", required=True, help="Output prefix")
    parser.add_argument("--chrom", default="11", help="Chromosome passed to convertf (default: 11)")
    return parser.parse_args()


def main():
    args = parse_args()

    prefix = Path(args.prefix)
    anno_path = Path(args.anno)
    ind_path = Path(f"{prefix}.ind")
    geno_path = Path(f"{prefix}.geno")
    snp_path = Path(f"{prefix}.snp")
    out = Path(args.out_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)

    for path in (anno_path, ind_path, geno_path, snp_path):
        if not path.exists():
            raise FileNotFoundError(path)

    anno = pd.read_csv(anno_path, sep="\t", dtype=str, low_memory=False)
    columns = list(anno.columns)

    gid_col = find_col(columns, startswith="genetic id")
    individual_col = find_col(columns, exact="individual id")
    lat_col = find_col(columns, exact="latitude")
    lon_col = find_col(columns, exact="longitude")
    date_col = find_col(columns, startswith="date mean in bp")
    snps_2m_col = find_col(
        columns,
        contains_all=["snps hit on autosomal targets", "enhance 2m"],
    )
    coverage_col = find_col(columns, startswith="mean coverage on 1.15m autosomal targets")

    anno["_date_bp"] = to_numeric(anno[date_col])
    anno["_snps_2m"] = to_numeric(anno[snps_2m_col])
    anno["_coverage"] = to_numeric(anno[coverage_col])

    ancient = anno.loc[anno["_date_bp"] > 0].copy()
    ancient = ancient.dropna(subset=[gid_col, individual_col])

    # AADR can contain multiple Genetic IDs for one archaeological individual.
    # Retain the representation with the most 2M autosomal target SNPs; use
    # coverage and Genetic ID only as deterministic tie-breakers.
    ancient["_snps_key"] = ancient["_snps_2m"].fillna(-1)
    ancient["_coverage_key"] = ancient["_coverage"].fillna(-1)
    selected = (
        ancient.sort_values(
            [individual_col, "_snps_key", "_coverage_key", gid_col],
            ascending=[True, False, False, True],
            kind="mergesort",
        )
        .drop_duplicates(subset=[individual_col], keep="first")
        .copy()
    )

    ind = pd.read_csv(
        ind_path,
        sep=r"\s+",
        header=None,
        names=["Genetic_ID", "Sex", "Group"],
        dtype=str,
    )

    keep = set(selected[gid_col])
    missing_in_ind = sorted(keep - set(ind["Genetic_ID"]))
    if missing_in_ind:
        raise RuntimeError(
            f"{len(missing_in_ind)} selected Genetic IDs are absent from {ind_path}; "
            f"first few: {missing_in_ind[:5]}"
        )

    # Match metadata to the sample order that convertf will retain.
    selected_order = ind.loc[ind["Genetic_ID"].isin(keep), "Genetic_ID"].tolist()
    selected = selected.set_index(gid_col).loc[selected_order].reset_index()

    metadata = selected[
        [individual_col, gid_col, lat_col, lon_col, "_date_bp"]
    ].copy()
    metadata.columns = ["Individual ID", "Genetic ID", "Latitude", "Longitude", "Date_BP"]
    metadata.to_csv(f"{out}.metadata.tsv", sep="\t", index=False)

    with open(f"{out}.samples.txt", "w") as handle:
        for sample in selected_order:
            handle.write(f"{sample}\n")

    # convertf removes individuals whose third .ind field is Ignore. Preserve the
    # original AADR order so genotype columns remain aligned with the source .geno.
    selected_ind = ind.copy()
    selected_ind.loc[~selected_ind["Genetic_ID"].isin(keep), "Group"] = "Ignore"
    selected_ind.to_csv(
        f"{out}.selected_full_order.ind",
        sep="\t",
        header=False,
        index=False,
    )

    # Do not set lopos/hipos: the downstream empirical background requires the
    # entire chromosome, not only the local PGA interval.
    convertf_par = f"""genotypename:    {geno_path}
snpname:         {snp_path}
indivname:       {out}.selected_full_order.ind
outputformat:    EIGENSTRAT
genotypeoutname: {out}.eigenstrat.geno
snpoutname:      {out}.eigenstrat.snp
indivoutname:    {out}.eigenstrat.ind
hashcheck:       NO
chrom:           {args.chrom}
"""
    with open(f"{out}.convertf.par", "w") as handle:
        handle.write(convertf_par)

    print(f"Ancient rows before Individual ID de-duplication: {len(ancient)}", file=sys.stderr)
    print(f"Selected unique ancient Individual IDs: {len(selected)}", file=sys.stderr)
    print(f"Wrote metadata: {out}.metadata.tsv", file=sys.stderr)
    print(f"Wrote convertf parameters: {out}.convertf.par", file=sys.stderr)


if __name__ == "__main__":
    main()
