#!/usr/bin/env python3
"""Temporal analysis of ancient PGA copy-number estimates.

The input table is the final ancient-sample CN table produced after PGA CN
genotyping. The script performs (i) a continuous-age OLS regression and (ii)
1-, 2-, and 3-kyr bin-size sensitivity analyses, with samples older than the
specified oldest-bin threshold pooled into a terminal bin.

The legacy input convention is retained: Date_Mean is negative for samples in
the past (for example, -9000 for 9 kya), so values increase toward the present.
"""

import argparse
from math import ceil
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import linregress


def configure_matplotlib():
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["mathtext.fontset"] = "custom"
    plt.rcParams["mathtext.rm"] = "Arial"
    plt.rcParams["mathtext.it"] = "Arial:italic"
    plt.rcParams["mathtext.bf"] = "Arial:bold"
    plt.rcParams["axes.unicode_minus"] = False


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze temporal trends in ancient PGA copy number."
    )
    parser.add_argument("--input", required=True, help="Input ancient PGA CN TSV")
    parser.add_argument(
        "--out-dir",
        default="ancient_copies_PGA_sensitivity",
        help="Output directory",
    )
    parser.add_argument(
        "--regions",
        nargs="+",
        default=["East Asia", "West Eurasia"],
        help="Region labels in the input table",
    )
    parser.add_argument(
        "--bin-sizes",
        nargs="+",
        type=int,
        default=[1000, 2000, 3000],
        help="Temporal bin sizes in years (default: 1000 2000 3000)",
    )
    parser.add_argument(
        "--oldest-bin-start",
        type=int,
        default=9000,
        help="Pool samples older than this age (years BP) into the terminal bin",
    )
    return parser.parse_args()


def load_data(input_path):
    df = pd.read_csv(input_path, sep="\t")
    required = ["Sample", "Region", "Date_Mean", "Pred_PGA_Total"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    df = df.copy()
    df["Date_Mean"] = pd.to_numeric(df["Date_Mean"], errors="coerce")
    df["Pred_PGA_Total"] = pd.to_numeric(df["Pred_PGA_Total"], errors="coerce")
    df = df.dropna(subset=["Region", "Date_Mean", "Pred_PGA_Total"])
    df["Date_kyr_toward_present"] = df["Date_Mean"] / 1000.0
    df["Age_BP"] = -df["Date_Mean"]
    return df


def get_common_ylim(df, regions):
    values = df.loc[df["Region"].isin(regions), "Pred_PGA_Total"].dropna()
    if values.empty:
        raise ValueError("No samples found for the requested regions")
    return 2, ceil(values.max() + 1)


def fit_linear(x, y):
    fit_df = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(fit_df) < 3 or fit_df["x"].nunique() < 2:
        return {
            "slope": np.nan,
            "intercept": np.nan,
            "p_value": np.nan,
            "r2": np.nan,
            "n": len(fit_df),
        }

    result = linregress(fit_df["x"], fit_df["y"])
    return {
        "slope": float(result.slope),
        "intercept": float(result.intercept),
        "p_value": float(result.pvalue),
        "r2": float(result.rvalue**2),
        "n": len(fit_df),
    }


def p_text(p_value):
    if pd.isna(p_value):
        return "$P$ = NA"
    return f"$P$ = {p_value:.2e}".replace("e-0", "e-").replace("e+0", "e+")


def add_fit_text(ax, fit):
    if pd.isna(fit["p_value"]):
        text = "Not enough data for linear fit"
    else:
        text = f"{p_text(fit['p_value'])}\n$R^2$ = {fit['r2']:.3f}"
    ax.text(0.04, 0.96, text, transform=ax.transAxes, fontsize=10, va="top", ha="left")


def make_bin_labels(oldest_bin_start, bin_size):
    if bin_size <= 0:
        raise ValueError("Bin size must be positive")
    if oldest_bin_start <= 0:
        raise ValueError("--oldest-bin-start must be positive")

    # Keep the terminal boundary fixed at the requested age. If the bin size
    # does not divide the boundary exactly (e.g. 2-kyr bins with a 9-kyr
    # boundary), the final closed interval is correspondingly shorter.
    upper_edges = list(range(bin_size, oldest_bin_start, bin_size))
    upper_edges.append(oldest_bin_start)
    upper_edges = np.asarray(upper_edges, dtype=float)
    bins = np.concatenate([[-np.inf], upper_edges, [np.inf]])

    labels = []
    previous = 0
    for upper in upper_edges:
        if previous == 0:
            labels.append(f"< {upper / 1000:g} Kya")
        else:
            labels.append(f"{previous / 1000:g}-{upper / 1000:g} Kya")
        previous = upper
    labels.append(f"> {oldest_bin_start / 1000:g} Kya")
    return bins, labels


def assign_time_bins(df_region, bin_size, oldest_bin_start):
    out = df_region.copy()
    bins, labels = make_bin_labels(oldest_bin_start, bin_size)
    plot_labels = list(reversed(labels))

    out["Date_label_sensitivity"] = pd.cut(
        out["Age_BP"], bins=bins, labels=labels, include_lowest=True
    )
    out = out.dropna(subset=["Date_label_sensitivity"]).copy()

    label_map = {label: index for index, label in enumerate(plot_labels)}
    out["Date_label_sensitivity"] = pd.Categorical(
        out["Date_label_sensitivity"].astype(str),
        categories=plot_labels,
        ordered=True,
    )
    out["Plot_Code"] = out["Date_label_sensitivity"].map(label_map)
    return out, plot_labels


def plot_unbinned(ax, df_region, region, y_limits):
    sns.scatterplot(
        ax=ax,
        data=df_region,
        x="Date_kyr_toward_present",
        y="Pred_PGA_Total",
        color="lightgray",
        edgecolor="none",
        alpha=0.8,
        s=34,
    )

    fit = fit_linear(df_region["Date_kyr_toward_present"], df_region["Pred_PGA_Total"])
    if not pd.isna(fit["p_value"]):
        sns.regplot(
            ax=ax,
            data=df_region,
            x="Date_kyr_toward_present",
            y="Pred_PGA_Total",
            scatter=False,
            color="black",
            ci=95,
            truncate=False,
            line_kws={"linestyle": "--", "linewidth": 1.5, "alpha": 0.6},
        )

    ax.set_ylim(y_limits)
    ax.set_title(region, fontsize=14)
    ax.set_xlabel("Time (kya)", fontsize=12)
    ax.set_ylabel("PGA copy number", fontsize=12)
    ax.tick_params(axis="both", labelsize=10)
    add_fit_text(ax, fit)
    sns.despine(ax=ax)
    return fit


def plot_binned(ax, df_region, bin_size, oldest_bin_start, y_limits):
    binned_df, labels = assign_time_bins(df_region, bin_size, oldest_bin_start)
    fit = fit_linear(binned_df["Plot_Code"], binned_df["Pred_PGA_Total"])

    sns.boxplot(
        ax=ax,
        data=binned_df,
        x="Date_label_sensitivity",
        y="Pred_PGA_Total",
        order=labels,
        showfliers=False,
        showcaps=False,
        fill=False,
    )
    sns.swarmplot(
        ax=ax,
        data=binned_df,
        x="Date_label_sensitivity",
        y="Pred_PGA_Total",
        order=labels,
        color="lightgray",
        alpha=0.8,
        size=4,
    )

    if not pd.isna(fit["p_value"]):
        sns.regplot(
            ax=ax,
            data=binned_df,
            x="Plot_Code",
            y="Pred_PGA_Total",
            scatter=False,
            color="black",
            ci=95,
            truncate=False,
            line_kws={"linestyle": "--", "linewidth": 1.5, "alpha": 0.6},
        )

    ticklabels = []
    for label in labels:
        count = (binned_df["Date_label_sensitivity"] == label).sum()
        ticklabels.append(f"{label}\n($n$ = {count})")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(ticklabels, rotation=45, ha="center", fontsize=10)

    ax.set_title(f"{bin_size / 1000:g}-kyr bins", fontsize=14)
    ax.set_xlabel("")
    ax.set_ylabel("PGA copy number", fontsize=12)
    ax.set_ylim(y_limits)
    add_fit_text(ax, fit)
    sns.despine(ax=ax)
    return fit, len(labels)


def plot_region(df, region, bin_sizes, oldest_bin_start, output_pdf, y_limits):
    df_region = df.loc[df["Region"] == region].copy()
    if df_region.empty:
        raise ValueError(f"No samples found for region: {region}")

    fit_rows = []
    with PdfPages(output_pdf) as pdf:
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        fit = plot_unbinned(ax, df_region, region, y_limits)
        fit_rows.append({
            "Region": region,
            "Analysis": "continuous",
            "Bin_size_years": 0,
            "N_samples": fit["n"],
            "N_time_bins": np.nan,
            "Slope": fit["slope"],
            "Intercept": fit["intercept"],
            "P_value": fit["p_value"],
            "R2": fit["r2"],
        })
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        ncols = len(bin_sizes)
        fig, axes = plt.subplots(1, ncols, figsize=(5.0 * ncols, 4.8), squeeze=False)
        for i, bin_size in enumerate(bin_sizes):
            fit, n_time_bins = plot_binned(
                axes[0, i], df_region, bin_size, oldest_bin_start, y_limits
            )
            fit_rows.append({
                "Region": region,
                "Analysis": "binned",
                "Bin_size_years": bin_size,
                "N_samples": fit["n"],
                "N_time_bins": n_time_bins,
                "Slope": fit["slope"],
                "Intercept": fit["intercept"],
                "P_value": fit["p_value"],
                "R2": fit["r2"],
            })
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    return fit_rows


def main():
    configure_matplotlib()
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load_data(args.input)
    y_limits = get_common_ylim(df, args.regions)

    all_fit_rows = []
    for region in args.regions:
        output_pdf = out_dir / f"ancient_copies.PGA_sensitivity.{region.replace(' ', '_')}.pdf"
        all_fit_rows.extend(
            plot_region(
                df,
                region,
                args.bin_sizes,
                args.oldest_bin_start,
                output_pdf,
                y_limits,
            )
        )
        print(f"Wrote {output_pdf}")

    summary_path = out_dir / "ancient_copies.PGA_sensitivity.fit_summary.tsv"
    pd.DataFrame(all_fit_rows).to_csv(summary_path, sep="\t", index=False)
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
