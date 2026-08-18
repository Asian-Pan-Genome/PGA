#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def read_table(path: str) -> pd.DataFrame:
    """Read comma- or tab-delimited input based on the filename suffix."""
    if str(path).lower().endswith(".csv"):
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def clean_species(value) -> str:
    """Standardize species names to 'Genus species'."""
    return str(value).strip().replace("_", " ")


def require_columns(df: pd.DataFrame, columns, table_name: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"{table_name} is missing required columns: {missing}")


def normalize_copy_table(cn: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize the updated assignment table.

    New input columns:
      assem, all, pga, pga_like, pag

    Standardized output aliases used by downstream scripts:
      assembly, PGA_all_CN, PGA_CN_primary,
      PGA_like_CN_primary, PAG_CN_primary

    The original columns are retained for traceability. Older tables containing
    assembly/PGA_CN_primary/PAG_CN_primary are also accepted when possible.
    """
    cn = cn.copy()

    if {"assem", "all", "pga", "pga_like", "pag"}.issubset(cn.columns):
        cn["assembly"] = cn["assem"].astype(str).str.strip()
        cn["PGA_all_CN"] = pd.to_numeric(cn["all"], errors="coerce")
        cn["PGA_CN_primary"] = pd.to_numeric(cn["pga"], errors="coerce")
        cn["PGA_like_CN_primary"] = pd.to_numeric(cn["pga_like"], errors="coerce")
        cn["PAG_CN_primary"] = pd.to_numeric(cn["pag"], errors="coerce")

        bad = cn[
            cn[["PGA_all_CN", "PGA_CN_primary", "PGA_like_CN_primary", "PAG_CN_primary"]]
            .isna()
            .any(axis=1)
        ]
        if not bad.empty:
            raise ValueError(
                "The updated copy-number table contains non-numeric values in "
                "all/pga/pga_like/pag. Example assemblies: "
                + ", ".join(bad["assembly"].head(5).astype(str))
            )

        expected_total = (
            cn["PGA_CN_primary"]
            + cn["PGA_like_CN_primary"]
            + cn["PAG_CN_primary"]
        )
        inconsistent = cn[cn["PGA_all_CN"] != expected_total]
        if not inconsistent.empty:
            raise ValueError(
                f"Found {len(inconsistent)} rows where all != pga + pga_like + pag. "
                "Please check the source table."
            )

    elif "assembly" in cn.columns and "PGA_CN_primary" in cn.columns:
        # Backward-compatible path for older source tables.
        cn["assembly"] = cn["assembly"].astype(str).str.strip()
        cn["PGA_CN_primary"] = pd.to_numeric(cn["PGA_CN_primary"], errors="coerce")

        if "PAG_CN_primary" in cn.columns:
            cn["PAG_CN_primary"] = pd.to_numeric(cn["PAG_CN_primary"], errors="coerce")
        else:
            cn["PAG_CN_primary"] = np.nan

        if "PGA_like_CN_primary" not in cn.columns:
            cn["PGA_like_CN_primary"] = np.nan
        if "PGA_all_CN" not in cn.columns:
            cn["PGA_all_CN"] = np.nan
    else:
        raise ValueError(
            "Unrecognized copy-number table. Expected updated columns "
            "assem/all/pga/pga_like/pag, or legacy columns including "
            "assembly/PGA_CN_primary."
        )

    if cn["assembly"].duplicated().any():
        dup = cn.loc[cn["assembly"].duplicated(keep=False), "assembly"].unique()
        raise ValueError(
            "Duplicated assembly identifiers in copy-number table: "
            + ", ".join(map(str, dup[:10]))
        )

    return cn


def safe_log1p(series: pd.Series, column_name: str) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    if (values.dropna() < 0).any():
        raise ValueError(f"Negative values found in {column_name}; log1p is undefined.")
    return np.log1p(values)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare a one-assembly-per-species PGA/diet/body-mass table for PGLS. "
            "The updated copy-number source must contain assem, all, pga, pga_like and pag."
        )
    )
    parser.add_argument("-a", "--assembly", required=True,
                        help="Assembly metadata TSV containing Directory Name and Species")
    parser.add_argument("-c", "--copy", required=True,
                        help="Copy-number TSV with columns: assem, all, pga, pga_like, pag")
    parser.add_argument("-d", "--diet", required=True,
                        help="Diet classification CSV/TSV containing Scientific and Diet_Class")
    parser.add_argument("-b", "--bodymass", required=True,
                        help="Body-mass TSV containing Scientific and BodyMass-Value")
    parser.add_argument("-o", "--out", required=True,
                        help="Output one-assembly-per-species TSV")
    parser.add_argument("--dup-out", default="duplicated_species.assembly_choice.tsv",
                        help="Output table documenting assembly selection for duplicated species")
    args = parser.parse_args()

    asm = read_table(args.assembly)
    cn = normalize_copy_table(read_table(args.copy))
    diet = read_table(args.diet)
    body = read_table(args.bodymass)

    require_columns(asm, ["Directory Name", "Species"], "assembly table")
    require_columns(diet, ["Scientific", "Diet_Class"], "diet table")
    require_columns(body, ["Scientific", "BodyMass-Value"], "body-mass table")

    asm = asm.copy()
    diet = diet.copy()
    body = body.copy()

    asm["Directory Name"] = asm["Directory Name"].astype(str).str.strip()
    asm["Species_clean"] = asm["Species"].apply(clean_species)
    diet["Scientific_clean"] = diet["Scientific"].apply(clean_species)
    body["Scientific_clean"] = body["Scientific"].apply(clean_species)
    body["BodyMass-Value"] = pd.to_numeric(body["BodyMass-Value"], errors="coerce")

    # The body-mass source may contain repeated records per species. Collapse
    # these to one species-level value using the median, preventing a
    # many-to-many merge and reducing sensitivity to duplicate source records.
    body_dup_rows = int(body.duplicated("Scientific_clean", keep=False).sum())
    body_dup_species = int(
        body.loc[body.duplicated("Scientific_clean", keep=False), "Scientific_clean"].nunique()
    )
    body = (
        body.groupby("Scientific_clean", as_index=False, dropna=False)["BodyMass-Value"]
        .median()
    )

    if diet["Scientific_clean"].duplicated().any():
        dup_names = diet.loc[
            diet["Scientific_clean"].duplicated(keep=False), "Scientific_clean"
        ].unique()
        raise ValueError(
            "Duplicated species in diet table after name cleaning: "
            + ", ".join(map(str, dup_names[:10]))
        )

    if asm["Directory Name"].duplicated().any():
        raise ValueError("Duplicated Directory Name values found in the assembly table.")

    assembly_copy = asm.merge(
        cn,
        left_on="Directory Name",
        right_on="assembly",
        how="inner",
        validate="one_to_one",
    )

    unmatched_assemblies = sorted(set(asm["Directory Name"]) - set(cn["assembly"]))
    unmatched_copy = sorted(set(cn["assembly"]) - set(asm["Directory Name"]))

    df = assembly_copy.merge(
        diet,
        left_on="Species_clean",
        right_on="Scientific_clean",
        how="inner",
        suffixes=("", "_diet"),
        validate="many_to_one",
    )

    df = df.merge(
        body[["Scientific_clean", "BodyMass-Value"]],
        on="Scientific_clean",
        how="left",
        validate="many_to_one",
    )

    numeric_cols = [
        "contig N50 (bp)",
        "scaffold N50 (bp)",
        "No. ancestral genes with intact ORF",
        "No. ancestral genes with missing sequences",
        "PGA_all_CN",
        "PGA_CN_primary",
        "PGA_like_CN_primary",
        "PAG_CN_primary",
        "BodyMass-Value",
        "Diet-Inv",
        "Diet-Vend",
        "Diet-Vect",
        "Diet-Vfish",
        "Diet-Vunk",
        "Diet-Scav",
        "Diet-Fruit",
        "Diet-Nect",
        "Diet-Seed",
        "Diet-PlantO",
    ]

    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    required_diet_ratio_cols = [
        "Diet-Inv", "Diet-Vend", "Diet-Vect", "Diet-Vfish", "Diet-Vunk",
        "Diet-Scav", "Diet-Fruit", "Diet-Nect", "Diet-Seed", "Diet-PlantO",
    ]
    require_columns(df, required_diet_ratio_cols, "merged diet table")

    df["plant_total_pct"] = (
        df["Diet-Fruit"] + df["Diet-Nect"] + df["Diet-Seed"] + df["Diet-PlantO"]
    )
    df["animal_total_pct"] = (
        df["Diet-Inv"] + df["Diet-Vend"] + df["Diet-Vect"] +
        df["Diet-Vfish"] + df["Diet-Vunk"] + df["Diet-Scav"]
    )
    df["plant_total_prop"] = df["plant_total_pct"] / 100
    df["plant_structural_prop"] = df["Diet-PlantO"] / 100

    df["log10_body_mass"] = np.where(
        df["BodyMass-Value"] > 0,
        np.log10(df["BodyMass-Value"]),
        np.nan,
    )

    df["PGA_all_log1p"] = safe_log1p(df["PGA_all_CN"], "PGA_all_CN")
    df["PGA_log1p"] = safe_log1p(df["PGA_CN_primary"], "PGA_CN_primary")
    df["PGA_like_log1p"] = safe_log1p(
        df["PGA_like_CN_primary"], "PGA_like_CN_primary"
    )
    df["PAG_log1p"] = safe_log1p(df["PAG_CN_primary"], "PAG_CN_primary")

    def diet_group(value: str) -> str:
        if value in ["Herbivore", "Frugivore", "Granivore", "Nectarivore"]:
            return "Plant_dominant"
        if value in ["Carnivore", "Insectivore"]:
            return "Animal_dominant"
        if value == "Omnivore":
            return "Omnivore"
        return "Other"

    df["Diet_Group"] = df["Diet_Class"].apply(diet_group)
    df["tip_label"] = df["Species_clean"].str.replace(" ", "_", regex=False)

    sort_cols = [
        "Species_clean",
        "contig N50 (bp)",
        "scaffold N50 (bp)",
        "No. ancestral genes with intact ORF",
        "No. ancestral genes with missing sequences",
    ]
    sort_ascending = [True, False, False, False, True]

    dup = df[df.duplicated("Species_clean", keep=False)].copy()
    dup = dup.sort_values(sort_cols, ascending=sort_ascending)
    Path(args.dup_out).parent.mkdir(parents=True, exist_ok=True)
    dup.to_csv(args.dup_out, sep="\t", index=False)

    rep = (
        df.sort_values(sort_cols, ascending=sort_ascending)
        .drop_duplicates("Species_clean", keep="first")
        .copy()
    )

    preferred_cols = [
        "Species_clean", "tip_label", "Directory Name", "assembly", "assem",
        "Assembly name", "NCBI accession", "Taxonomic Lineage",
        "contig N50 (bp)", "scaffold N50 (bp)",
        "No. ancestral genes with intact ORF",
        "No. ancestral genes with missing sequences",
        "PGA_all_CN", "PGA_all_log1p",
        "PGA_CN_primary", "PGA_log1p",
        "PGA_like_CN_primary", "PGA_like_log1p",
        "PAG_CN_primary", "PAG_log1p",
        "all", "pga", "pga_like", "pag",
        "Diet_Class", "Diet_Group",
        "Diet-Inv", "Diet-Vend", "Diet-Vect", "Diet-Vfish", "Diet-Vunk",
        "Diet-Scav", "Diet-Fruit", "Diet-Nect", "Diet-Seed", "Diet-PlantO",
        "plant_total_pct", "plant_total_prop", "plant_structural_prop",
        "animal_total_pct", "BodyMass-Value", "log10_body_mass",
    ]
    preferred_cols = [c for c in preferred_cols if c in rep.columns]
    other_cols = [c for c in rep.columns if c not in preferred_cols]
    rep = rep[preferred_cols + other_cols]

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    rep.to_csv(args.out, sep="\t", index=False)

    print("Assembly metadata rows:", len(asm))
    print("Copy-number rows:", len(cn))
    print("Assembly-copy matched rows:", len(assembly_copy))
    print("Assembly-copy-diet-body rows:", len(df))
    print("Final one-assembly-per-species rows:", len(rep))
    print("Species with missing body mass:", int(rep["BodyMass-Value"].isna().sum()))
    print(
        "Body-mass duplicate rows collapsed:", body_dup_rows,
        "across", body_dup_species, "species"
    )
    print("Unmatched assembly metadata rows:", len(unmatched_assemblies))
    print("Unmatched copy-number rows:", len(unmatched_copy))
    if unmatched_assemblies:
        print("Example assembly metadata IDs without CN:", ", ".join(unmatched_assemblies[:5]))
    if unmatched_copy:
        print("Example CN assembly IDs without metadata:", ", ".join(unmatched_copy[:5]))
    print("Output:", args.out)
    print("Duplicated-species assembly-choice table:", args.dup_out)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
