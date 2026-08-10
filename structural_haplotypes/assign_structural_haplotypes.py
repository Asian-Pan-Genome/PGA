#!/usr/bin/env python3

import argparse

import pandas as pd


BUNDLE_COLUMNS = [
    "contig",
    "start",
    "end",
    "bundle_info",
    "sample_hap",
    "bundle_id",
    "bundle_type",
    "bundle_path",
]

SAMPLE_ALIASES = {
    "YAO.Mat": "YAO.hap2",
    "YAO.Pat": "YAO.hap1",
    "CHM13v2.hap0": "CHM13v2",
    "HG002v1.1.hap1": "HG002.hap1",
    "HG002v1.1.hap2": "HG002.hap2",
    "CN1v1.hap0": "CN1v1",
    "GRCh38.hap0": "GRCh38",
}

# Corrections preserved from the original analysis. These paths contain hybrid
# bundle configurations whose bundle-7 and/or paralog annotations would
# otherwise be double-counted when generating A/B/C/X/Y labels.
PATH_COUNT_CORRECTIONS = {
    "0+2+3+4+5+2+3+4+5+2+3+4+5+2+6+7+4+5+2+6+7+1": {"C": -1},
    "0+2+3+4+5+2+3+4+5+2+3+4+5+2+6+7+4+5+2+6+7+4+5+2+6+7+1": {"C": -2},
    "0+2+3+4+5+2+3+4+5+2+3+4+5+2+3+4+5+2+6+7+4+5+2+3+4+5+2+6+7+1": {"C": -1},
    "0+2+3+4+5+2+3+4+5+2+3+4+7+1": {"A": -1, "C": -1},
    "0+2+3+4+5+2+6+7+4+5+2+6+7+1": {"C": -1},
    "0+2+3+4+5+2+3+4+7+1": {"A": -1, "C": -1},
    "0+2+3+4+7+1": {"B": -1, "C": -1},
}

# Final label swap used in the current analysis.
LABEL_SWAP = {
    "A1B2C1.1": "A1B2C1.2",
    "A1B2C1.2": "A1B2C1.1",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Assign PGA structural-haplotype labels from curated principal bundles and paralog annotations."
    )
    parser.add_argument("principal_bundle_bed", help="Curated principal-bundle BED.")
    parser.add_argument(
        "paralog_table",
        help="TSV containing Sample, Hap, Source and PGAs columns.",
    )
    parser.add_argument("output_tsv", help="Output structural-haplotype table.")
    return parser.parse_args()


def normalize_sample_hap(sample, hap):
    sample_hap = f"{sample}.{hap}"
    sample_hap = SAMPLE_ALIASES.get(sample_hap, sample_hap)
    if sample_hap.startswith("apr"):
        sample_hap = sample_hap.replace("hap", "")
    return sample_hap


def base_label_counts(pgas, bundle_path):
    counts = {
        "A": pgas.count("PGA34A"),
        "B": pgas.count("PGA34B"),
        "C": bundle_path.count("7"),
        "X": bundle_path.count("2+3+4+7"),
        "Y": bundle_path.count("2+6+7+4"),
    }

    for key, delta in PATH_COUNT_CORRECTIONS.get(bundle_path, {}).items():
        counts[key] += delta

    if any(value < 0 for value in counts.values()):
        raise ValueError(
            f"Negative SH component count for PGAs={pgas}, bundle_path={bundle_path}: {counts}"
        )

    return counts


def label_prefix(counts):
    return "".join(f"{key}{counts[key]}" for key in ["A", "B", "C", "X", "Y"] if counts[key] > 0)


def main():
    args = parse_args()

    bundles = pd.read_csv(
        args.principal_bundle_bed,
        sep="\t",
        header=None,
        names=BUNDLE_COLUMNS,
        dtype=str,
    )

    bundles["sample_hap"] = bundles["sample_hap"].replace(SAMPLE_ALIASES)
    bundle_paths = (
        bundles[["sample_hap", "bundle_path"]]
        .drop_duplicates()
        .copy()
    )
    duplicated = bundle_paths["sample_hap"].duplicated(keep=False)
    if duplicated.any():
        bad = bundle_paths.loc[duplicated, "sample_hap"].unique().tolist()
        raise ValueError(f"Multiple bundle paths found for haplotypes: {bad[:10]}")

    paralogs = pd.read_csv(args.paralog_table, sep="\t")
    required = {"Sample", "Hap", "Source", "PGAs"}
    missing = required - set(paralogs.columns)
    if missing:
        raise ValueError(f"Missing columns in paralog table: {sorted(missing)}")

    paralogs = paralogs[["Sample", "Hap", "Source", "PGAs"]].copy()
    paralogs["sample_hap"] = [
        normalize_sample_hap(sample, hap)
        for sample, hap in zip(paralogs["Sample"], paralogs["Hap"])
    ]

    merged = paralogs.merge(bundle_paths, on="sample_hap", how="left", validate="one_to_one")
    if merged["bundle_path"].isna().any():
        missing_haps = merged.loc[merged["bundle_path"].isna(), "sample_hap"].tolist()
        raise ValueError(
            "No principal-bundle path found for: "
            + ", ".join(missing_haps[:10])
            + (" ..." if len(missing_haps) > 10 else "")
        )

    prefixes = []
    signatures = []
    for _, row in merged.iterrows():
        counts = base_label_counts(str(row["PGAs"]), str(row["bundle_path"]))
        prefixes.append(label_prefix(counts))
        signatures.append((str(row["PGAs"]), str(row["bundle_path"])))

    merged["sh_prefix"] = prefixes
    merged["_signature"] = signatures

    # If the same A/B/C/X/Y composition has multiple exact structural
    # configurations, append .1, .2, ... using the sorted exact signatures.
    signature_order = {}
    for prefix, group in merged.groupby("sh_prefix", sort=False):
        unique_signatures = sorted(set(group["_signature"]))
        signature_order[prefix] = {
            signature: i + 1 for i, signature in enumerate(unique_signatures)
        }

    final_labels = []
    for _, row in merged.iterrows():
        prefix = row["sh_prefix"]
        mapping = signature_order[prefix]
        if len(mapping) == 1:
            label = prefix
        else:
            label = f"{prefix}.{mapping[row['_signature']]}"
        final_labels.append(LABEL_SWAP.get(label, label))

    merged["final_hap_label"] = final_labels
    merged = merged.drop(columns=["sh_prefix", "_signature"])

    merged.to_csv(args.output_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
