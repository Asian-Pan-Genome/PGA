#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_POPULATIONS = ["AFR", "AMR", "EAS", "EUR", "SAS", "CSA"]


def load_association_file(path, population, chrom, p_threshold):
    df = pd.read_csv(path, sep=r"\s+")

    required = {"SNP", "BETA", "P"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")

    if "TEST" in df.columns:
        df = df[df["TEST"] == "ADD"].copy()

    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df["BETA"] = pd.to_numeric(df["BETA"], errors="coerce")
    df = df[df["P"] < p_threshold].copy()

    if df.empty:
        return None

    variant = df["SNP"].astype(str).str.split(":", n=3, expand=True)
    if variant.shape[1] != 4:
        raise ValueError(
            f"Expected SNP IDs in CHROM:POS:REF:ALT format in {path}"
        )

    df["CHR"] = str(chrom)
    df["pos"] = pd.to_numeric(variant[1], errors="raise").astype(int)
    df["ref"] = variant[2]
    df["alt"] = variant[3]
    df["Population"] = population

    return df[["Population", "CHR", "pos", "BETA", "P", "ref", "alt"]]


def main():
    parser = argparse.ArgumentParser(
        description="Merge genome-wide significant SNP associations with PGA34A copy number."
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing chr1/ through chr22/ association results.",
    )
    parser.add_argument(
        "--output",
        default="all_populations_aggregated.P5e-8.tsv",
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=5e-8,
    )
    parser.add_argument(
        "--populations",
        nargs="+",
        default=DEFAULT_POPULATIONS,
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    results = []

    for population in args.populations:
        for chrom in range(1, 23):
            path = (
                input_dir
                / f"chr{chrom}"
                / f"{population}.chr{chrom}.assoc_Final_Result.assoc.linear"
            )

            if not path.exists():
                raise FileNotFoundError(path)

            df = load_association_file(
                path=path,
                population=population,
                chrom=chrom,
                p_threshold=args.p_threshold,
            )
            if df is not None:
                results.append(df)

    if not results:
        raise RuntimeError(
            f"No associations with P < {args.p_threshold:g} were found."
        )

    merged = pd.concat(results, ignore_index=True)

    population_order = {pop: i for i, pop in enumerate(args.populations)}
    merged["_population_order"] = merged["Population"].map(population_order)
    merged["_chrom_order"] = pd.to_numeric(merged["CHR"], errors="raise")
    merged = merged.sort_values(
        ["_chrom_order", "pos", "ref", "alt", "_population_order"]
    )

    aggregated = (
        merged.groupby(["CHR", "pos", "ref", "alt"], sort=False, as_index=False)
        .agg(
            BETA=("BETA", lambda x: ";".join(map(str, x))),
            P=("P", lambda x: ";".join(map(str, x))),
            Population=("Population", lambda x: ";".join(x)),
        )
    )

    aggregated.to_csv(args.output, sep="\t", index=False)

    print(f"Significant population-variant associations: {len(merged):,}")
    print(f"Unique significant variants: {len(aggregated):,}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
