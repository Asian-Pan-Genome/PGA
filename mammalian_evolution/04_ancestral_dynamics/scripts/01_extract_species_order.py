#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


MAMMAL_ORDERS = set(
    """
Afrosoricida
Artiodactyla
Carnivora
Cetacea
Cetartiodactyla
Chiroptera
Cingulata
Dasyuromorphia
Dermoptera
Didelphimorphia
Diprotodontia
Eulipotyphla
Hyracoidea
Lagomorpha
Macroscelidea
Microbiotheria
Monotremata
Notoryctemorphia
Paucituberculata
Peramelemorphia
Perissodactyla
Pholidota
Pilosa
Primates
Proboscidea
Rodentia
Scandentia
Sirenia
Tubulidentata
""".split()
)


def infer_order(lineage):
    if pd.isna(lineage):
        return "Unknown"

    parts = [item.strip() for item in str(lineage).split(";")]
    for item in parts:
        if item in MAMMAL_ORDERS:
            return item

    anchors = [
        "Laurasiatheria",
        "Euarchontoglires",
        "Afrotheria",
        "Xenarthra",
        "Marsupialia",
        "Metatheria",
        "Eutheria",
    ]
    for anchor in anchors:
        if anchor in parts:
            index = parts.index(anchor)
            if index + 1 < len(parts):
                return parts[index + 1]

    return "Unknown"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract a species-to-order table from the prepared mammalian association table."
    )
    parser.add_argument("-i", "--input", required=True, help="Prepared species association TSV")
    parser.add_argument("-o", "--output", required=True, help="Species/order TSV")
    parser.add_argument(
        "--check-output",
        help="Optional audit TSV retaining species names and full taxonomic lineages",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    frame = pd.read_csv(args.input, sep="\t", dtype=str)

    required = {"Taxonomic Lineage"}
    if "tip_label" not in frame.columns:
        required.add("Species")
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Input is missing required columns: {missing}")

    if "tip_label" in frame.columns:
        species = frame["tip_label"].fillna("")
    else:
        species = frame["Species"].fillna("").str.replace(" ", "_", regex=False)

    orders = frame["Taxonomic Lineage"].apply(infer_order)
    result = pd.DataFrame({"species": species, "order": orders})
    result = result.drop_duplicates().sort_values(["order", "species"])

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, sep="\t", index=False)
    print(f"Written: {output}")

    if args.check_output:
        check_output = Path(args.check_output)
        check_output.parent.mkdir(parents=True, exist_ok=True)
        check = pd.DataFrame(
            {
                "species": species,
                "Species": frame["Species"] if "Species" in frame.columns else species,
                "order": orders,
                "Taxonomic Lineage": frame["Taxonomic Lineage"],
            }
        )
        check = check.drop_duplicates().sort_values(["order", "species"])
        check.to_csv(check_output, sep="\t", index=False)
        print(f"Written: {check_output}")

    print("\nOrder counts:")
    print(result["order"].value_counts().sort_index())
    n_unknown = int((result["order"] == "Unknown").sum())
    print(f"\nUnknown order rows: {n_unknown}")


if __name__ == "__main__":
    main()
