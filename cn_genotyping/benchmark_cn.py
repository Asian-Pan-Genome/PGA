#!/usr/bin/env python3

import argparse

import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, mean_squared_error


TARGETS = ("PGA", "PGA34", "PGA34A", "PGA34B", "PGA5")


def parse_args():
    parser = argparse.ArgumentParser(description="Benchmark predicted PGA copy numbers against assembly-resolved truth labels.")
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--train-label", default=None)
    parser.add_argument("--predict-label", default=None)
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


def main():
    args = parse_args()
    predictions = pd.read_csv(args.predictions, sep="\t")
    truth = load_truth(args.truth)
    merged = pd.merge(predictions, truth, on="Sample", how="inner")
    if merged.empty:
        raise ValueError("No overlapping samples between predictions and truth.")

    rows = []
    for target in TARGETS:
        pred_column = f"Pred_{target}"
        if pred_column not in merged.columns:
            raise ValueError(f"Missing prediction column: {pred_column}")

        y_true = merged[target].to_numpy(dtype=float)
        y_pred = merged[pred_column].to_numpy(dtype=float)
        row = {
            "Gene": target,
            "Accuracy": accuracy_score(y_true, y_pred),
            "RMSE": float(np.sqrt(mean_squared_error(y_true, y_pred))),
            "N": len(y_true),
        }
        if args.train_label is not None:
            row["Train"] = args.train_label
        if args.predict_label is not None:
            row["Predict"] = args.predict_label
        rows.append(row)

    columns = []
    if args.train_label is not None:
        columns.append("Train")
    if args.predict_label is not None:
        columns.append("Predict")
    columns.extend(["Gene", "Accuracy", "RMSE", "N"])

    output = pd.DataFrame(rows)[columns]
    output.to_csv(args.out, sep="\t", index=False)
    print(output.to_string(index=False))
    print(f"Output: {args.out}")


if __name__ == "__main__":
    main()
