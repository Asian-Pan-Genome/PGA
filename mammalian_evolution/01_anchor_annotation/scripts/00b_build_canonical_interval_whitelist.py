#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import pandas as pd


BED12_COLS = [
    "chrom", "start", "end", "name", "score", "strand",
    "thickStart", "thickEnd", "rgb", "blockCount",
    "blockSizes", "blockStarts"
]


def split_aliases(x):
    """
    Convert comma-separated aliases into lowercase set.
    Example:
      VPS37C,Vps37c,vps37c
    """
    return set(i.strip().lower() for i in str(x).split(",") if i.strip())


def parse_gene_label(name):
    parts = str(name).split("#")
    if len(parts) >= 2:
        return parts[1]
    return "NA"


def parse_name(name, ref_prefix):
    """
    Input examples:
      ENST00000312403.10#PGA5
      hg38.ENST00000312403.10#PGA5#8382
      HLbosTau10.ENSBTAT00000013786#PAG7
      mm10.ENSMUST00000025647#Pga5

    Output:
      source_transcript_id = hg38.ENST00000312403.10
      gene_label           = PGA5
      source_name          = hg38.ENST00000312403.10#PGA5
    """
    parts = str(name).split("#")

    transcript = parts[0]
    gene = parts[1] if len(parts) >= 2 else "NA"

    if transcript.startswith(ref_prefix + "."):
        source_transcript_id = transcript
    else:
        source_transcript_id = ref_prefix + "." + transcript

    source_name = source_transcript_id + "#" + gene

    return source_transcript_id, gene, source_name


def choose_common_chrom(a1, a2):
    """
    Choose chromosome/contig shared by both anchors.
    If multiple exist, choose the one with the largest total number of anchor rows.
    """
    common = sorted(set(a1["chrom"]) & set(a2["chrom"]))

    if not common:
        return None

    if len(common) == 1:
        return common[0]

    best_chrom = None
    best_n = -1

    for chrom in common:
        n = int((a1["chrom"] == chrom).sum()) + int((a2["chrom"] == chrom).sum())
        if n > best_n:
            best_n = n
            best_chrom = chrom

    return best_chrom


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build canonical VPS37C-VWCE source whitelist from a reference "
            "TOGA transcript BED. Anchor matching is case-insensitive."
        )
    )

    parser.add_argument(
        "--bed",
        required=True,
        help="Reference transcript BED12, e.g. hg38.toga.transcripts.bed"
    )

    parser.add_argument(
        "--ref-prefix",
        required=True,
        help="Reference prefix used in TOGA names, e.g. hg38, mm10, HLbosTau10, HLeleMaxInd3A"
    )

    parser.add_argument(
        "--anchor1",
        default="VPS37C,Vps37c,vps37c",
        help="Comma-separated aliases for anchor1. Default: VPS37C,Vps37c,vps37c"
    )

    parser.add_argument(
        "--anchor2",
        default="VWCE,Vwce,vwce",
        help="Comma-separated aliases for anchor2. Default: VWCE,Vwce,vwce"
    )

    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Output whitelist TSV"
    )

    args = parser.parse_args()

    anchor1_aliases = split_aliases(args.anchor1)
    anchor2_aliases = split_aliases(args.anchor2)

    df = pd.read_csv(
        args.bed,
        sep="\t",
        header=None,
        names=BED12_COLS,
        dtype={
            "chrom": str,
            "name": str,
            "score": str,
            "strand": str,
            "rgb": str,
            "blockSizes": str,
            "blockStarts": str,
        }
    )

    if df.empty:
        raise ValueError(f"Empty BED file: {args.bed}")

    for c in ["start", "end", "thickStart", "thickEnd", "blockCount"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

    df["gene_label"] = df["name"].apply(parse_gene_label)
    df["gene_label_lower"] = df["gene_label"].astype(str).str.lower()

    a1 = df[df["gene_label_lower"].isin(anchor1_aliases)].copy()
    a2 = df[df["gene_label_lower"].isin(anchor2_aliases)].copy()

    if len(a1) == 0:
        all_near = sorted(
            x for x in set(df["gene_label"])
            if "vps" in str(x).lower() or "37" in str(x).lower()
        )
        sys.stderr.write(
            "[ERROR] Cannot find anchor1 gene label.\n"
            f"[ERROR] anchor1 aliases: {','.join(sorted(anchor1_aliases))}\n"
            f"[ERROR] possible related labels: {','.join(all_near) if all_near else 'NA'}\n"
        )
        raise ValueError(f"Cannot find anchor1 gene_label in {args.bed}")

    if len(a2) == 0:
        all_near = sorted(
            x for x in set(df["gene_label"])
            if "vwce" in str(x).lower() or "vwc" in str(x).lower()
        )
        sys.stderr.write(
            "[ERROR] Cannot find anchor2 gene label.\n"
            f"[ERROR] anchor2 aliases: {','.join(sorted(anchor2_aliases))}\n"
            f"[ERROR] possible related labels: {','.join(all_near) if all_near else 'NA'}\n"
        )
        raise ValueError(f"Cannot find anchor2 gene_label in {args.bed}")

    chrom = choose_common_chrom(a1, a2)

    if chrom is None:
        raise ValueError(
            "Anchor1 and anchor2 are not on the same chromosome/contig. "
            "This reference may have a split or naming problem."
        )

    a1_sub = a1[a1["chrom"] == chrom]
    a2_sub = a2[a2["chrom"] == chrom]

    a1_start = int(a1_sub["start"].min())
    a1_end = int(a1_sub["end"].max())
    a2_start = int(a2_sub["start"].min())
    a2_end = int(a2_sub["end"].max())

    region_start = min(a1_start, a2_start)
    region_end = max(a1_end, a2_end)

    # Use overlap with canonical anchor interval.
    # This includes the two anchors and genes between them.
    interval = df[
        (df["chrom"] == chrom)
        & (df["end"] > region_start)
        & (df["start"] < region_end)
    ].copy()

    rows = []

    for _, r in interval.iterrows():
        source_transcript_id, gene_label, source_name = parse_name(
            r["name"],
            args.ref_prefix
        )

        gene_lower = str(gene_label).lower()

        if gene_lower in anchor1_aliases or gene_lower in anchor2_aliases:
            target_status = "anchor"
        else:
            target_status = "candidate_interval_gene"

        rows.append({
            "ref_prefix": args.ref_prefix,
            "chrom": r["chrom"],
            "start": int(r["start"]),
            "end": int(r["end"]),
            "strand": r["strand"],
            "blockCount": int(r["blockCount"]),
            "raw_name": r["name"],
            "source_transcript_id": source_transcript_id,
            "gene_label": gene_label,
            "source_name": source_name,
            "target_status": target_status,
            "canonical_interval": f"{chrom}:{region_start}-{region_end}",
        })

    out = pd.DataFrame(rows).drop_duplicates()

    out = out.sort_values([
        "chrom",
        "start",
        "end",
        "target_status",
        "gene_label",
        "source_transcript_id",
    ])

    out.to_csv(args.out, sep="\t", index=False)

    gene_labels = ",".join(sorted(set(out["gene_label"])))

    print(f"[INFO] ref_prefix: {args.ref_prefix}")
    print(f"[INFO] canonical interval: {chrom}:{region_start}-{region_end}")
    print(f"[INFO] anchor1 labels found: {','.join(sorted(set(a1_sub['gene_label'])))}")
    print(f"[INFO] anchor2 labels found: {','.join(sorted(set(a2_sub['gene_label'])))}")
    print(f"[INFO] total rows: {len(out)}")
    print(f"[INFO] gene labels: {gene_labels}")
    print(f"[INFO] output: {args.out}")


if __name__ == "__main__":
    main()
