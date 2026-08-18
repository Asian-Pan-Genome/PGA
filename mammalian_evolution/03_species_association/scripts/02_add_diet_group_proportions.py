#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


DIET_COLS = [
    "Diet-Inv", "Diet-Vend", "Diet-Vect", "Diet-Vfish", "Diet-Vunk",
    "Diet-Scav", "Diet-Fruit", "Diet-Nect", "Diet-Seed", "Diet-PlantO",
]

CN_ALIASES = {
    "PGA_CN_primary": ["PGA_CN_primary", "pga"],
    "PGA_like_CN_primary": ["PGA_like_CN_primary", "pga_like"],
    "PAG_CN_primary": ["PAG_CN_primary", "pag"],
    "PGA_all_CN": ["PGA_all_CN", "all"],
    "pga": ["pga", "PGA_CN_primary"],
    "pga_like": ["pga_like", "PGA_like_CN_primary"],
    "pag": ["pag", "PAG_CN_primary"],
    "all": ["all", "PGA_all_CN"],
}


def read_table(path: str) -> pd.DataFrame:
    if str(path).lower().endswith(".csv"):
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def resolve_cn_column(df: pd.DataFrame, requested: str) -> str:
    """Resolve standardized or raw updated CN column names."""
    if requested in df.columns:
        return requested

    for candidate in CN_ALIASES.get(requested, []):
        if candidate in df.columns:
            return candidate

    lower_to_actual = {str(c).lower(): c for c in df.columns}
    if requested.lower() in lower_to_actual:
        return lower_to_actual[requested.lower()]

    raise ValueError(
        f"Copy-number column '{requested}' was not found. Available supported columns include: "
        "PGA_CN_primary/pga, PGA_like_CN_primary/pga_like, "
        "PAG_CN_primary/pag and PGA_all_CN/all."
    )


def high_order_diet(diet_class, nectar_mode="drop") -> str:
    """
    High-order diet groups:
      Plant_dominant = Herbivore + Frugivore + Granivore
      Omnivore
      Carnivore
      Insectivore

    Nectarivore is dropped by default or included in Plant_dominant when
    --nectar-mode plant is selected.
    """
    if pd.isna(diet_class):
        return "Unknown"

    diet_class = str(diet_class).strip()

    if diet_class in ["Herbivore", "Frugivore", "Granivore"]:
        return "Plant_dominant"
    if diet_class == "Omnivore":
        return "Omnivore"
    if diet_class == "Carnivore":
        return "Carnivore"
    if diet_class == "Insectivore":
        return "Insectivore"
    if diet_class == "Nectarivore":
        return "Plant_dominant" if nectar_mode == "plant" else "Drop_Nectarivore"
    return "Unknown"


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Add high-order diet groups and food-ratio variables for PGLS. "
            "Supports updated CN columns pga, pga_like, pag and all, as well as "
            "their standardized aliases."
        )
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input table from 01_prepare_species_association_input.py")
    parser.add_argument("-o", "--output", required=True,
                        help="Output TSV with high-order diet and food-ratio variables")
    parser.add_argument(
        "--cn-col",
        default="PGA_CN_primary",
        help=(
            "Copy-number response column. Supported examples: PGA_CN_primary or pga; "
            "PGA_like_CN_primary or pga_like; PAG_CN_primary or pag; PGA_all_CN or all. "
            "Default: PGA_CN_primary"
        ),
    )
    parser.add_argument(
        "--nectar-mode",
        choices=["drop", "plant"],
        default="drop",
        help="Drop Nectarivore or include it in Plant_dominant. Default: drop",
    )
    args = parser.parse_args()

    df = read_table(args.input)
    cn_col = resolve_cn_column(df, args.cn_col)

    required = ["Diet_Class", "BodyMass-Value", cn_col] + DIET_COLS
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    for col in DIET_COLS + ["BodyMass-Value", cn_col]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    if (df[cn_col].dropna() < 0).any():
        raise ValueError(f"Negative copy numbers found in {cn_col}.")

    df["High_Order_Diet"] = df["Diet_Class"].apply(
        lambda value: high_order_diet(value, nectar_mode=args.nectar_mode)
    )
    df["Use_High_Order_Diet"] = ~df["High_Order_Diet"].isin(
        ["Unknown", "Drop_Nectarivore"]
    )

    # Base plant proportion: fruit + seed + other plant material.
    plant_pct = df["Diet-Fruit"] + df["Diet-Seed"] + df["Diet-PlantO"]
    if args.nectar_mode == "plant":
        # Keep the continuous diet ratio consistent with classifying Nectarivore
        # as Plant_dominant.
        plant_pct = plant_pct + df["Diet-Nect"]

    df["plant_dominant_pct"] = plant_pct
    df["carnivore_pct"] = (
        df["Diet-Vend"] + df["Diet-Vect"] + df["Diet-Vfish"] +
        df["Diet-Vunk"] + df["Diet-Scav"]
    )
    df["insectivore_pct"] = df["Diet-Inv"]
    df["nectarivore_pct"] = df["Diet-Nect"]

    df["plant_dominant_prop"] = df["plant_dominant_pct"] / 100
    df["carnivore_prop"] = df["carnivore_pct"] / 100
    df["insectivore_prop"] = df["insectivore_pct"] / 100
    df["nectarivore_prop"] = df["nectarivore_pct"] / 100
    df["diet_total_pct"] = df[DIET_COLS].sum(axis=1)

    df["CNV_source_col"] = cn_col
    df["CNV"] = df[cn_col]
    df["CNV_log1p"] = np.log1p(df["CNV"])
    df["log10_body_mass"] = np.where(
        df["BodyMass-Value"] > 0,
        np.log10(df["BodyMass-Value"]),
        np.nan,
    )

    if "tip_label" not in df.columns:
        if "Species_clean" in df.columns:
            df["tip_label"] = df["Species_clean"].astype(str).str.replace(
                " ", "_", regex=False
            )
        elif "Scientific" in df.columns:
            df["tip_label"] = df["Scientific"].astype(str).str.replace(
                " ", "_", regex=False
            )
        else:
            raise ValueError("No tip_label, Species_clean or Scientific column found.")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)

    print("Resolved CN column:", cn_col)
    print("Output:", args.output)
    print("\nHigh-order diet counts:")
    print(df["High_Order_Diet"].value_counts(dropna=False).to_string())
    print("\nSpecies used in high-order diet PGLS:")
    print(df["Use_High_Order_Diet"].value_counts(dropna=False).to_string())
    print("\nMissing selected CN values:", int(df["CNV"].isna().sum()))
    print("Missing body-mass values:", int(df["BodyMass-Value"].isna().sum()))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
