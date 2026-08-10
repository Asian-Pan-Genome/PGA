#!/usr/bin/env python3

import argparse
from collections import Counter

import pandas as pd


RAW_COLUMNS = ["contig", "start", "end", "bundle_info"]
OUTPUT_COLUMNS = [
    "contig",
    "start",
    "end",
    "bundle_info",
    "sample_hap",
    "bundle_id",
    "bundle_type",
    "bundle_path",
]

PGA34_RAW_BUNDLES = ["3", "7", "10", "11", "4"]
PGA5_RAW_BUNDLES = ["3", "6", "9", "13", "5"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert raw PGR-TK bundle decomposition into curated PGA principal bundles."
    )
    parser.add_argument("input_bed", help="Raw BED from pgr-pbundle-decomp.")
    parser.add_argument("output_bed", help="Curated principal-bundle BED.")
    return parser.parse_args()


def classify_gene_part(gene_part, sample_hap):
    counts = Counter(gene_part["bundle_id"])
    score34 = sum(counts[b] for b in PGA34_RAW_BUNDLES)
    score5 = sum(counts[b] for b in PGA5_RAW_BUNDLES)
    bundle_ids = set(gene_part["bundle_id"])

    if score34 > score5 and ("5" not in bundle_ids or counts["4"] == 2):
        return "PGA34"

    if score5 > score34 and ("4" not in bundle_ids or counts["5"] == 2):
        return "PGA5"

    if {"7", "5"}.issubset(bundle_ids) and "9" not in bundle_ids and "13" not in bundle_ids:
        return "PGA34_PGA5"

    if {"6", "4"}.issubset(bundle_ids) and "10" not in bundle_ids and "11" not in bundle_ids:
        return "PGA5_PGA34"

    raise ValueError(
        f"{sample_hap}: unexpected raw bundle combination in gene part: "
        f"{gene_part['bundle_id'].tolist()}"
    )


def add_record(records, group, start, end, bundle_id, bundle_info, bundle_type, sample_hap):
    records.append(
        {
            "contig": group.iloc[0]["contig"],
            "start": int(start),
            "end": int(end),
            "bundle_info": bundle_info,
            "sample_hap": sample_hap,
            "bundle_id": str(bundle_id),
            "bundle_type": bundle_type,
        }
    )


def curate_haplotype(group, sample_hap):
    group = group.reset_index(drop=True)
    records = []

    # Upstream boundary.
    upstream = group.loc[:3, "bundle_id"].tolist()
    case1 = ["1", "3"]
    case2 = ["1", "8", "3"]
    case3 = ["1", "2", "1", "3"]
    case4 = ["1", "2", "8", "3"]

    if upstream[:2] == case1 or upstream[:3] == case2 or upstream == case3:
        last_1 = len(upstream) - 1 - upstream[::-1].index("1")
        add_record(
            records,
            group,
            group["start"].min(),
            group.loc[last_1, "end"],
            "0",
            "0:52:0:0:51:U",
            "U",
            sample_hap,
        )
        gene_scan_start = last_1 + 1
    elif upstream == case4:
        add_record(
            records,
            group,
            group["start"].min(),
            group.loc[0, "end"],
            "0",
            "0:52:0:0:5:U",
            "U",
            sample_hap,
        )
        add_record(
            records,
            group,
            group.loc[1, "start"],
            group.loc[2, "end"],
            "5",
            "5:52:0:0:51:U",
            "U",
            sample_hap,
        )
        gene_scan_start = 3
    else:
        raise ValueError(
            f"{sample_hap}: unexpected upstream raw bundle pattern: {upstream}"
        )

    # Locate gene-sized raw bundle blocks.
    gene_region = group.loc[gene_scan_start:]
    start_indices = gene_region.index[gene_region["bundle_id"] == "3"].tolist()
    end_indices = []

    for row_index, row in gene_region.iterrows():
        if row_index + 1 >= len(group):
            continue
        next_bundle = group.loc[row_index + 1, "bundle_id"]
        if row["bundle_id"] in {"4", "5"} and next_bundle in {"12", "2", "0"}:
            end_indices.append(row_index)

    if len(start_indices) != len(end_indices):
        raise ValueError(
            f"{sample_hap}: unequal numbers of gene starts and ends: "
            f"{start_indices} vs {end_indices}"
        )

    gene_parts = [group.loc[start:end] for start, end in zip(start_indices, end_indices)]
    gene_classes = [classify_gene_part(part, sample_hap) for part in gene_parts]

    n_pga34 = gene_classes.count("PGA34")
    n_pga5 = gene_classes.count("PGA5")
    n_34_5 = gene_classes.count("PGA34_PGA5")
    n_5_34 = gene_classes.count("PGA5_PGA34")

    intergenic_parts = [
        group.loc[end_indices[i] + 1 : start_indices[i + 1] - 1]
        for i in range(len(start_indices) - 1)
    ]

    for i, (start_idx, end_idx, gene_class) in enumerate(
        zip(start_indices, end_indices, gene_classes)
    ):
        full_gene = group.loc[start_idx:end_idx]

        common_type = "R" if len(gene_parts) > 1 else "U"
        add_record(
            records,
            group,
            full_gene.iloc[0]["start"],
            full_gene.iloc[0]["end"],
            "2",
            f"2:34:0:1:33:{common_type}",
            common_type,
            sample_hap,
        )

        gene_part = group.loc[start_idx + 1 : end_idx]

        if gene_class == "PGA34":
            first_4 = gene_part.index[gene_part["bundle_id"] == "4"][0]

            type3 = "R" if n_pga34 + n_34_5 > 1 else "U"
            add_record(
                records,
                group,
                gene_part.loc[: first_4 - 1, "start"].min(),
                gene_part.loc[: first_4 - 1, "end"].max(),
                "3",
                f"3:94:0:0:93:{type3}",
                type3,
                sample_hap,
            )

            type4 = "R" if n_pga34 + n_34_5 + n_5_34 > 1 else "U"
            add_record(
                records,
                group,
                gene_part.loc[first_4:, "start"].min(),
                gene_part.loc[first_4:, "end"].max(),
                "4",
                f"4:19:0:0:18:{type4}",
                type4,
                sample_hap,
            )

        elif gene_class == "PGA5":
            first_5 = gene_part.index[gene_part["bundle_id"] == "5"][0]

            type6 = "R" if n_pga5 + n_5_34 > 1 else "U"
            add_record(
                records,
                group,
                gene_part.loc[: first_5 - 1, "start"].min(),
                gene_part.loc[: first_5 - 1, "end"].max(),
                "6",
                f"6:94:0:0:93:{type6}",
                type6,
                sample_hap,
            )

            type7 = "R" if n_pga5 + n_34_5 > 1 else "U"
            add_record(
                records,
                group,
                gene_part.loc[first_5:, "start"].min(),
                gene_part.loc[first_5:, "end"].max(),
                "7",
                f"7:19:0:0:18:{type7}",
                type7,
                sample_hap,
            )

        elif gene_class == "PGA34_PGA5":
            first_4 = gene_part.index[gene_part["bundle_id"] == "4"][0]

            type3 = "R" if n_34_5 + n_pga34 > 1 else "U"
            add_record(
                records,
                group,
                gene_part.loc[: first_4 - 1, "start"].min(),
                gene_part.loc[: first_4 - 1, "end"].max(),
                "3",
                f"3:94:0:0:93:{type3}",
                type3,
                sample_hap,
            )

            type4 = "R" if n_pga34 + n_34_5 + n_5_34 > 1 else "U"
            add_record(
                records,
                group,
                gene_part.loc[first_4, "start"],
                gene_part.loc[first_4, "end"],
                "4",
                f"4:19:0:0:18:{type4}",
                type4,
                sample_hap,
            )

            type7 = "R" if n_pga5 + n_34_5 > 1 else "U"
            add_record(
                records,
                group,
                gene_part.loc[first_4 + 1 :, "start"].min(),
                gene_part.loc[first_4 + 1 :, "end"].max(),
                "7",
                f"7:19:0:0:18:{type7}",
                type7,
                sample_hap,
            )

        elif gene_class == "PGA5_PGA34":
            gene_part = gene_part.reset_index(drop=True)
            first_5 = gene_part.index[gene_part["bundle_id"] == "5"][0]

            type6 = "R" if n_pga5 + n_5_34 > 1 else "U"
            add_record(
                records,
                group,
                gene_part.loc[: first_5 - 1, "start"].min(),
                gene_part.loc[: first_5 - 1, "end"].max(),
                "6",
                f"6:94:0:0:93:{type6}",
                type6,
                sample_hap,
            )

            type7_4 = "R" if n_pga34 + n_34_5 + n_5_34 > 1 else "U"
            add_record(
                records,
                group,
                gene_part.loc[first_5, "start"],
                gene_part.loc[first_5, "end"],
                "7",
                f"7:19:0:0:18:{type7_4}",
                type7_4,
                sample_hap,
            )
            add_record(
                records,
                group,
                gene_part.loc[first_5 + 1 :, "start"].min(),
                gene_part.loc[first_5 + 1 :, "end"].max(),
                "4",
                f"4:19:0:0:18:{type7_4}",
                type7_4,
                sample_hap,
            )

        if i < len(intergenic_parts):
            intergenic = intergenic_parts[i]
            type5 = "R" if len(intergenic_parts) > 1 else "U"
            add_record(
                records,
                group,
                intergenic["start"].min(),
                intergenic["end"].max(),
                "5",
                f"5:94:0:0:93:{type5}",
                type5,
                sample_hap,
            )

    # Downstream boundary.
    downstream = group.iloc[-2:]["bundle_id"].tolist()
    if len(downstream) == 2 and downstream[-2] == "12":
        start = group.iloc[-2:]["start"].min()
        end = group.iloc[-2:]["end"].max()
    else:
        start = group.iloc[-1]["start"]
        end = group.iloc[-1]["end"]

    add_record(
        records,
        group,
        start,
        end,
        "1",
        "1:61:0:0:60:U",
        "U",
        sample_hap,
    )

    return records


def main():
    args = parse_args()

    raw = pd.read_csv(
        args.input_bed,
        sep="\t",
        header=None,
        names=RAW_COLUMNS,
        dtype={"contig": str, "bundle_info": str},
    )
    raw["sample_hap"] = raw["contig"].str.split("::", n=1).str[0]
    raw["bundle_id"] = raw["bundle_info"].str.split(":").str[0]

    records = []
    for sample_hap, group in raw.groupby("sample_hap", sort=False):
        records.extend(curate_haplotype(group, sample_hap))

    out = pd.DataFrame(records)
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)

    for sample_hap, group in out.groupby("sample_hap", sort=False):
        bundle_path = "+".join(group["bundle_id"].astype(str))
        out.loc[group.index, "bundle_path"] = bundle_path

    out[OUTPUT_COLUMNS].to_csv(
        args.output_bed,
        sep="\t",
        index=False,
        header=False,
    )


if __name__ == "__main__":
    main()
