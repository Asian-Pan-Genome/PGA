#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "Arial"
plt.rcParams["mathtext.it"] = "Arial:italic"
plt.rcParams["mathtext.bf"] = "Arial:bold"
plt.rcParams["axes.unicode_minus"] = False


PALETTE = {
    "AFR": "#E64B35",
    "AFR-Other": "#E64B35",
    "AFR-E&S": "#FFD5C2",
    "AFR-W": "#741B11",
    "AFR-NE": "#FF9900",

    "ARB": "#4DBBD5",
    "EUR": "#3C5488",
    "WAS": "#8491B4",
    "MID": "#8491B4",

    "CAS": "#D9CCF2",
    "CSA": "#8467BD",
    "SAS": "#311B92",

    "EAS": "#00A087",
    "SEA": "#91D1C2",
    "OCE": "#71C6B1",

    "AMR": "#7F7F7F",
}


SUPERPOP_ORDER = [
    "AFR-Other", "AFR-E&S", "AFR-W", "AFR-NE",
    "ARB", "EUR", "WAS",
    "CAS", "CSA", "SAS",
    "EAS", "SEA", "OCE",
    "AMR",
]


PANEL_STYLE = {
    "1KGP panel": {
        "marker": "x",
        "s": 16,
        "alpha": 0.4,
        "linewidths": 0.7,
        "zorder": 1,
    },
    "Assembly panel": {
        "marker": "o",
        "s": 46,
        "alpha": 0.8,
        "linewidths": 0.0,
        "zorder": 3,
    },
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot joint diploid PCA of assembly panel and 1KGP panel."
    )
    parser.add_argument("--eigenvec", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--pc-x", type=int, default=1)
    parser.add_argument("--pc-y", type=int, default=2)
    parser.add_argument("--flip-x", action="store_true")
    parser.add_argument("--flip-y", action="store_true")
    parser.add_argument("--fig-width", type=float, default=8.2)
    parser.add_argument("--fig-height", type=float, default=5.8)
    parser.add_argument(
        "--kg-marker",
        choices=["x", "s"],
        default="x",
        help="Marker for 1KGP panel: x or s. Default: x.",
    )
    return parser.parse_args()


def read_eigenvec(eigenvec):
    df = pd.read_csv(eigenvec, sep=r"\s+", header=None, dtype=str)

    if df.iloc[0, 0].lower() in {"fid", "#fid"}:
        df = pd.read_csv(eigenvec, sep=r"\s+", header=0)
        cols = list(df.columns)
        df = df.rename(columns={cols[0]: "FID", cols[1]: "IID"})
        pc_cols = [f"PC{i}" for i in range(1, df.shape[1] - 1)]
        df.columns = ["FID", "IID"] + pc_cols
    else:
        n_pc = df.shape[1] - 2
        df.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, n_pc + 1)]

    for col in df.columns:
        if col.startswith("PC"):
            df[col] = pd.to_numeric(df[col], errors="raise")

    return df


def read_eigenval(eigenvec):
    eigenval = eigenvec.replace(".eigenvec", ".eigenval")
    if not os.path.exists(eigenval):
        return {}

    vals = []
    with open(eigenval) as f:
        for line in f:
            if line.strip():
                vals.append(float(line.strip()))

    total = sum(vals)
    if total <= 0:
        return {}

    return {i + 1: vals[i] / total * 100 for i in range(len(vals))}


def make_superpop_legend_handles(df, panel, marker, marker_size):
    handles = []
    sub_panel = df[df["Panel"] == panel]

    for sp in SUPERPOP_ORDER:
        n = sub_panel[sub_panel["Superpopulation"] == sp].shape[0]
        if n == 0:
            continue

        handles.append(
            Line2D(
                [0],
                [0],
                marker=marker,
                linestyle="",
                color=PALETTE.get(sp, "#BDBDBD"),
                markerfacecolor=PALETTE.get(sp, "#BDBDBD"),
                markeredgecolor=PALETTE.get(sp, "#BDBDBD"),
                markersize=marker_size,
                label=f"{sp} ($n$ = {n})",
            )
        )

    return handles


def main():
    args = parse_args()

    pcx = f"PC{args.pc_x}"
    pcy = f"PC{args.pc_y}"

    pca = read_eigenvec(args.eigenvec)
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str)

    required = {"IID", "Panel", "Superpopulation"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"metadata missing required columns: {sorted(missing)}")

    meta["Superpopulation"] = meta["Superpopulation"].replace(
        {"NA": pd.NA, "na": pd.NA, "N/A": pd.NA, "nan": pd.NA, "": pd.NA}
    )

    df = pca.merge(meta, on="IID", how="left")

    n_missing = df["Superpopulation"].isna().sum()
    if n_missing > 0:
        print(
            f"Warning: {n_missing} PCA samples do not have metadata and will be dropped.",
            file=sys.stderr,
        )

    df = df.dropna(subset=["Panel", "Superpopulation"]).copy()

    if pcx not in df.columns or pcy not in df.columns:
        raise ValueError(f"{pcx} or {pcy} not found in eigenvec.")

    if args.flip_x:
        df[pcx] = -df[pcx]
    if args.flip_y:
        df[pcy] = -df[pcy]

    merged_tsv = f"{args.out}.PC{args.pc_x}_vs_PC{args.pc_y}.merged.tsv"
    df.to_csv(merged_tsv, sep="\t", index=False)

    print("Samples used in plot:")
    print(df.groupby(["Panel", "Superpopulation"]).size())
    print(f"Merged PCA table written to: {merged_tsv}")

    var_exp = read_eigenval(args.eigenvec)

    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))

    panel_style = PANEL_STYLE.copy()
    panel_style["1KGP panel"] = panel_style["1KGP panel"].copy()
    panel_style["1KGP panel"]["marker"] = args.kg_marker

    for panel in ["1KGP panel", "Assembly panel"]:
        if panel not in set(df["Panel"]):
            continue

        style = panel_style[panel]

        for sp in SUPERPOP_ORDER:
            sub = df[(df["Panel"] == panel) & (df["Superpopulation"] == sp)]
            if sub.empty:
                continue

            color = PALETTE.get(sp, "#BDBDBD")

            if panel == "1KGP panel" and style["marker"] == "x":
                ax.scatter(
                    sub[pcx],
                    sub[pcy],
                    c=color,
                    marker="x",
                    s=style["s"],
                    alpha=style["alpha"],
                    linewidths=style["linewidths"],
                    zorder=style["zorder"],
                )
            else:
                ax.scatter(
                    sub[pcx],
                    sub[pcy],
                    c=color,
                    marker=style["marker"],
                    s=style["s"],
                    alpha=style["alpha"],
                    edgecolors="none",
                    linewidths=0,
                    zorder=style["zorder"],
                )

    xlabel = f"PC{args.pc_x}"
    ylabel = f"PC{args.pc_y}"

    if args.pc_x in var_exp:
        xlabel += f" ({var_exp[args.pc_x]:.2f}%)"
    if args.pc_y in var_exp:
        ylabel += f" ({var_exp[args.pc_y]:.2f}%)"

    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.tick_params(axis="both", labelsize=12)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    kg_handles = make_superpop_legend_handles(
        df=df,
        panel="1KGP panel",
        marker=args.kg_marker,
        marker_size=6.0,
    )

    asm_handles = make_superpop_legend_handles(
        df=df,
        panel="Assembly panel",
        marker="o",
        marker_size=7.0,
    )

    leg1 = ax.legend(
        handles=kg_handles,
        title="1KGP panel",
        loc="center left",
        bbox_to_anchor=(1.02, 0.68),
        frameon=False,
        fontsize=9,
        title_fontsize=11,
        handletextpad=0.4,
        borderaxespad=0.0,
    )
    ax.add_artist(leg1)

    ax.legend(
        handles=asm_handles,
        title="Assembly panel",
        loc="center left",
        bbox_to_anchor=(1.02, 0.25),
        frameon=False,
        fontsize=9,
        title_fontsize=11,
        handletextpad=0.4,
        borderaxespad=0.0,
    )

    fig.tight_layout()

    out_pdf = f"{args.out}.PC{args.pc_x}_vs_PC{args.pc_y}.pdf"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
