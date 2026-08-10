#!/usr/bin/env python3

import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a local PGA gene track for extracted haplotype sequences."
    )
    parser.add_argument("region_hit", help="Table describing extracted PGA regions.")
    parser.add_argument(
        "id_list",
        help="TSV containing Sample, Hap and GFF File columns.",
    )
    parser.add_argument(
        "paralog_table",
        help="TSV containing Sample, Hap and PGAs columns.",
    )
    parser.add_argument("output_bed", help="Output BED4 gene track.")
    return parser.parse_args()


def parse_source(src):
    parts = src.split(".")
    if len(parts) == 3:
        sample, hap = parts[:2]
        if sample.startswith("apr"):
            hap = f"hap{hap}"
        elif sample == "YAO":
            hap = "hap1" if hap == "Pat" else "hap2"
    elif len(parts) == 2:
        sample = parts[0]
        hap = "hap0"
    else:
        raise ValueError(f"Unexpected src format: {src}")
    return sample, hap


def parse_attributes(text):
    attributes = {}
    for item in text.rstrip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value
    return attributes


def main():
    args = parse_args()

    region_hit = pd.read_csv(args.region_hit, sep="\t")
    id_list = pd.read_csv(args.id_list, sep="\t")
    paralogs = pd.read_csv(args.paralog_table, sep="\t").dropna(subset=["PGAs"])

    required_id = {"Sample", "Hap", "GFF File"}
    if not required_id.issubset(id_list.columns):
        raise ValueError(f"id_list must contain: {sorted(required_id)}")

    required_paralog = {"Sample", "Hap", "PGAs"}
    if not required_paralog.issubset(paralogs.columns):
        raise ValueError(f"paralog_table must contain: {sorted(required_paralog)}")

    with open(args.output_bed, "w") as out:
        for _, region in region_hit.iterrows():
            sample, hap = parse_source(str(region["src"]))
            contig = str(region["ctg_name"]).split()[0]
            region_start = int(region["ctg_bgn"])
            region_end = int(region["ctg_end"])

            matches = id_list[
                (id_list["Sample"] == sample) & (id_list["Hap"] == hap)
            ]
            if len(matches) != 1:
                raise ValueError(
                    f"Expected one GFF for {sample}.{hap}, found {len(matches)}"
                )
            gff_file = matches.iloc[0]["GFF File"]

            paralog_match = paralogs[
                (paralogs["Sample"] == sample) & (paralogs["Hap"] == hap)
            ]
            gene_labels = None
            if not paralog_match.empty:
                if len(paralog_match) != 1:
                    raise ValueError(
                        f"Expected one paralog row for {sample}.{hap}, found {len(paralog_match)}"
                    )
                gene_labels = str(paralog_match.iloc[0]["PGAs"]).split("+")

            gene_index = 0
            with open(gff_file) as gff:
                for line in gff:
                    if line.startswith("#"):
                        continue

                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 9 or fields[2] != "gene":
                        continue

                    attributes = parse_attributes(fields[8])
                    gene_name = attributes.get("gene_name", "")
                    if not gene_name.startswith("PGA"):
                        continue

                    start = int(fields[3]) - 1
                    end = int(fields[4])

                    if start < region_start or end > region_end:
                        raise ValueError(
                            f"PGA gene outside extracted region for {contig}: "
                            f"{start}-{end} in {gff_file}"
                        )

                    if gene_labels is None:
                        label = "PGA"
                    else:
                        if gene_index >= len(gene_labels):
                            raise ValueError(
                                f"More PGA genes in GFF than labels in PGAs for {sample}.{hap}"
                            )
                        label = gene_labels[gene_index]

                    out.write(
                        f"{contig}\t{start - region_start}\t{end - region_start}\t{label}\n"
                    )
                    gene_index += 1

            if gene_labels is not None and gene_index != len(gene_labels):
                raise ValueError(
                    f"PGA gene count mismatch for {sample}.{hap}: "
                    f"GFF={gene_index}, PGAs={len(gene_labels)}"
                )


if __name__ == "__main__":
    main()
