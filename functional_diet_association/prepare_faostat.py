#!/usr/bin/env python3

"""Prepare the 1961-1980 FAOSTAT plant-derived protein fraction table."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


REGION_RENAME = {
    "United States of America": "USA",
    "United Kingdom of Great Britain and Northern Ireland": "UK",
    "Russian Federation": "Russia",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate country-level plant-derived protein fractions from "
            "FAOSTAT Food Balances protein-supply data and average them over 1961-1980."
        )
    )
    parser.add_argument("--input", type=Path, required=True, help="FAOSTAT Food Balances CSV")
    parser.add_argument("--output", type=Path, required=True, help="Output CSV")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    df = pd.read_csv(args.input)
    required = {"Area", "Year", "Element", "Item", "Value"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"FAOSTAT input is missing required columns: {sorted(missing)}")

    df = df.loc[
        df["Element"].eq("Protein supply quantity (g/capita/day)")
        & df["Year"].between(1961, 1980)
    ].copy()
    if df.empty:
        raise ValueError("No 1961-1980 protein-supply records were found in the FAOSTAT input.")

    pivot = df.pivot_table(
        index=["Area", "Year"],
        columns="Item",
        values="Value",
        aggfunc="mean",
    ).reset_index()

    required_items = {"Grand Total", "Animal Products"}
    missing_items = required_items - set(pivot.columns)
    if missing_items:
        raise ValueError(f"FAOSTAT input is missing required Items: {sorted(missing_items)}")

    total_protein = pd.to_numeric(pivot["Grand Total"], errors="coerce")
    animal_protein = pd.to_numeric(pivot["Animal Products"], errors="coerce")
    plant_protein = total_protein - animal_protein

    pivot["plant_ratio"] = np.where(
        total_protein > 0,
        plant_protein / total_protein,
        np.nan,
    )

    baseline = (
        pivot.dropna(subset=["plant_ratio"])
        .groupby("Area", as_index=False)["plant_ratio"]
        .mean()
        .rename(columns={"Area": "region"})
    )
    baseline["region"] = baseline["region"].replace(REGION_RENAME)
    baseline = baseline.sort_values("region").reset_index(drop=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    baseline.to_csv(args.output, index=False)
    print(f"Wrote {len(baseline)} regions to {args.output}")


if __name__ == "__main__":
    main()
