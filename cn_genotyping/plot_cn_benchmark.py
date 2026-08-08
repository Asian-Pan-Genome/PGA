#!/usr/bin/env python3

import argparse
import math

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


TARGETS = ("PGA", "PGA34", "PGA34A", "PGA34B", "PGA5")


def parse_args():
    parser = argparse.ArgumentParser(description="Plot true-versus-predicted PGA copy-number heatmaps.")
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--out-prefix", required=True)
    return parser.parse_args()


def load_truth(path):
    truth = pd.read_csv(path, sep="\t")
    required = {"Sample", "PGA34A", "PGA34B", "PGA5"}
    missing = required - set(truth.columns)
    if missing:
        raise ValueError(f"Truth file is missing columns: {', '.join(sorted(missing))}")

    truth["Sample"] = truth["Sample"].astype(str)
    counts = truth["Sample"].value_counts()
    if counts.max() > 1:
        valid = counts[counts == 2].index
        truth = (
            truth[truth["Sample"].isin(valid)]
            .groupby("Sample", as_index=False)[["PGA34A", "PGA34B", "PGA5"]]
            .sum()
        )
    else:
        truth = truth[["Sample", "PGA34A", "PGA34B", "PGA5"]].copy()

    truth["PGA34"] = truth["PGA34A"] + truth["PGA34B"]
    truth["PGA"] = truth["PGA34"] + truth["PGA5"]
    return truth


def plot_panel(ax, merged, target):
    pred_column = f"Pred_{target}"
    values = merged[[target, pred_column]].dropna().astype(int)
    max_cn = int(max(values[target].max(), values[pred_column].max()))
    states = range(max_cn + 1)
    table = pd.crosstab(values[target], values[pred_column], dropna=False)
    table = table.reindex(index=states, columns=states, fill_value=0)

    sns.heatmap(table, annot=True, fmt="d", cbar=True, square=True, ax=ax)
    ax.set_title(target)
    ax.set_xlabel("Predicted CN")
    ax.set_ylabel("True CN")


def main():
    args = parse_args()
    predictions = pd.read_csv(args.predictions, sep="\t")
    truth = load_truth(args.truth)
    merged = pd.merge(predictions, truth, on="Sample", how="inner")
    if merged.empty:
        raise ValueError("No overlapping samples between predictions and truth.")

    ncols = 3
    nrows = math.ceil(len(TARGETS) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 5.5 * nrows))
    axes = axes.flatten()

    for ax, target in zip(axes, TARGETS):
        plot_panel(ax, merged, target)
    for ax in axes[len(TARGETS):]:
        ax.axis("off")

    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf", bbox_inches="tight")
    fig.savefig(f"{args.out_prefix}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
