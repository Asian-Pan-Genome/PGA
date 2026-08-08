#!/usr/bin/env python3

import argparse
from functools import partial
from multiprocessing import Pool

import joblib
import numpy as np
import pandas as pd
import pysam
from sklearn.linear_model import LinearRegression


PARALOGS = ("PGA34A", "PGA34B", "PGA5")
MODEL_SPECS = {
    "PGA": ("Raw_PGA", "PGA"),
    "PGA34": ("Raw_PGA34", "PGA34"),
    "PGA5": ("Raw_PGA5", "PGA5"),
    "PGA34B": ("Raw_PGA34B", "PGA34B"),
}


def parse_manifest(path):
    table = pd.read_csv(path, sep="\t")
    required = {"CHROM", "POS", "REF", "ALT", "SPECIFIC", "GF_SPECIFICITY", "PGA34B"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Manifest is missing columns: {', '.join(sorted(missing))}")

    sites = {}
    for _, row in table.iterrows():
        specific = [x for x in str(row["SPECIFIC"]).split(",") if x and x != "nan"]
        if "PGA34B" not in specific:
            continue

        gt_text = str(row["PGA34B"]).split(":", 1)[0]
        if gt_text not in {"0", "1"}:
            continue

        sites[(str(row["CHROM"]), int(row["POS"]))] = {
            "ref": str(row["REF"]),
            "alt": str(row["ALT"]),
            "gt": int(gt_text),
            "weight": float(row["GF_SPECIFICITY"]),
            "shared_count": len(specific),
        }
    return sites


def read_depth(path):
    table = pd.read_csv(path, sep="\t")
    required = {"REF_Copy", "Baseline_depth", "PGA_depth", "PGA34_depth", "PGA5_depth"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Depth file {path} is missing columns: {', '.join(sorted(missing))}")
    if table.empty:
        raise ValueError(f"Depth file is empty: {path}")
    return table.iloc[0]


def extract_one_sample(sample_info, reference_fasta, sites):
    sample, alignment_path, depth_path = sample_info

    try:
        depth = read_depth(depth_path)
        ref_copy = int(depth["REF_Copy"])
        per_copy_depth = float(depth["Baseline_depth"]) / 2.0
        if per_copy_depth <= 0:
            raise ValueError("Baseline_depth must be > 0")

        raw_pga = float(depth["PGA_depth"]) / per_copy_depth * ref_copy
        raw_pga34 = float(depth["PGA34_depth"]) / per_copy_depth * max(ref_copy - 1, 0)
        raw_pga5 = float(depth["PGA5_depth"]) / per_copy_depth
    except Exception as exc:
        print(f"[ERROR] {sample}: could not read depth features: {exc}")
        return None

    ratios = []
    weights = []
    mode = "rc" if alignment_path.lower().endswith(".cram") else "rb"
    kwargs = {"reference_filename": reference_fasta} if mode == "rc" else {}

    try:
        with pysam.AlignmentFile(alignment_path, mode, **kwargs) as alignment, \
             pysam.FastaFile(reference_fasta) as reference:
            for (chrom, pos), info in sites.items():
                counts = {"A": 0, "C": 0, "G": 0, "T": 0}
                total = 0

                columns = alignment.pileup(
                    chrom,
                    pos - 1,
                    pos,
                    truncate=True,
                    stepper="samtools",
                    fastafile=reference,
                    min_base_quality=0,
                    ignore_overlaps=True,
                    ignore_orphans=False,
                )

                for column in columns:
                    for pileup_read in column.pileups:
                        if pileup_read.is_del or pileup_read.query_position is None:
                            continue
                        base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                        if base in counts:
                            counts[base] += 1
                        total += 1
                    break

                if total == 0:
                    continue

                target_allele = info["ref"] if info["gt"] == 0 else info["alt"]
                target_count = counts.get(target_allele, 0)
                ratios.append(target_count / (total * info["shared_count"]))
                weights.append(info["weight"])
    except Exception as exc:
        print(f"[ERROR] {sample}: could not extract informative-site support: {exc}")
        return None

    if ratios:
        pga34b_ratio = float(np.average(ratios, weights=weights))
    else:
        pga34b_ratio = 0.0

    return {
        "Sample": sample,
        "Raw_PGA": raw_pga,
        "Raw_PGA34": raw_pga34,
        "Raw_PGA5": raw_pga5,
        "Raw_PGA34B_Ratio": pga34b_ratio,
        "Raw_PGA34B": raw_pga * pga34b_ratio,
    }


def extract_features(sample_list, reference_fasta, manifest, output, threads):
    sites = parse_manifest(manifest)
    if not sites:
        raise ValueError("No PGA34B-informative sites were found in the manifest.")

    tasks = []
    with open(sample_list) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[0] == "Sample":
                continue
            if len(fields) < 3:
                raise ValueError("Sample list requires three tab-delimited columns: Sample, Alignment, Depth.")
            tasks.append(tuple(fields[:3]))

    worker = partial(extract_one_sample, reference_fasta=reference_fasta, sites=sites)
    results = []
    with Pool(processes=threads) as pool:
        for result in pool.imap_unordered(worker, tasks):
            if result is not None:
                results.append(result)

    pd.DataFrame(results).sort_values("Sample").to_csv(output, sep="\t", index=False)
    print(f"Features: {output} ({len(results)} samples)")


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


class ConstrainedGMM:
    def __init__(self, max_iter=100):
        self.max_iter = max_iter
        self.mu = None
        self.pi = None
        self.sigma = None
        self.offset = None
        self.scaler = None

    def fit(self, x_labeled, y_labeled, x_unlabeled=None):
        x_labeled = np.asarray(x_labeled, dtype=float)
        y_labeled = np.asarray(y_labeled, dtype=float)
        if len(x_labeled) == 0:
            raise ValueError("At least one labeled sample is required.")

        self.scaler = LinearRegression()
        self.scaler.fit(x_labeled.reshape(-1, 1), y_labeled)

        labeled_corrected = self.scaler.predict(x_labeled.reshape(-1, 1))
        if x_unlabeled is not None and len(x_unlabeled) > 0:
            x_pool = self.scaler.predict(np.asarray(x_unlabeled, dtype=float).reshape(-1, 1))
        else:
            x_pool = labeled_corrected

        min_cn = max(0, int(np.floor(x_pool.min())))
        max_cn = int(np.ceil(x_pool.max())) + 1
        self.mu = np.arange(min_cn, max_cn + 1, dtype=float)
        self.pi = np.ones(len(self.mu), dtype=float) / len(self.mu)
        self.offset = 0.0

        residuals = labeled_corrected - y_labeled
        self.sigma = max(float(np.std(residuals)), 0.1)

        x_column = x_pool.reshape(-1, 1)
        for _ in range(self.max_iter):
            means = self.mu + self.offset
            log_probs = (
                -0.5 * ((x_column - means) / self.sigma) ** 2
                - np.log(self.sigma)
                - 0.5 * np.log(2.0 * np.pi)
            )
            log_weighted = log_probs + np.log(np.clip(self.pi, 1e-300, None))
            max_log = np.max(log_weighted, axis=1, keepdims=True)
            log_sum = max_log + np.log(np.sum(np.exp(log_weighted - max_log), axis=1, keepdims=True))
            gamma = np.exp(log_weighted - log_sum)

            self.pi = np.sum(gamma, axis=0) / len(x_pool)
            self.offset = float(np.sum(gamma * (x_column - self.mu)) / len(x_pool))

            means = self.mu + self.offset
            variance = float(np.sum(gamma * (x_column - means) ** 2) / len(x_pool))
            self.sigma = max(np.sqrt(variance), 0.05)

    def predict(self, x):
        x = np.asarray(x, dtype=float)
        corrected = self.scaler.predict(x.reshape(-1, 1))
        x_column = corrected.reshape(-1, 1)
        means = self.mu + self.offset
        log_probs = -0.5 * ((x_column - means) / self.sigma) ** 2
        log_weighted = log_probs + np.log(np.clip(self.pi, 1e-300, None))
        return self.mu[np.argmax(log_weighted, axis=1)].astype(int)


def fit_models(features_path, truth_path, model_out, unlabeled_features_path=None):
    features = pd.read_csv(features_path, sep="\t")
    truth = load_truth(truth_path)
    labeled = pd.merge(features, truth, on="Sample", how="inner")
    if labeled.empty:
        raise ValueError("No samples overlap between the feature and truth files.")

    unlabeled = None
    if unlabeled_features_path:
        unlabeled = pd.read_csv(unlabeled_features_path, sep="\t")

    models = {}
    for name, (feature_column, truth_column) in MODEL_SPECS.items():
        model = ConstrainedGMM()
        x_unlabeled = None if unlabeled is None else unlabeled[feature_column].to_numpy(dtype=float)
        model.fit(
            labeled[feature_column].to_numpy(dtype=float),
            labeled[truth_column].to_numpy(dtype=float),
            x_unlabeled=x_unlabeled,
        )
        models[name] = model
        print(f"{name}: n_labeled={len(labeled)}, sigma={model.sigma:.4f}, offset={model.offset:.4f}")

    joblib.dump(models, model_out)
    print(f"Models: {model_out}")


def predict(features_path, model_path, output):
    features = pd.read_csv(features_path, sep="\t")
    models = joblib.load(model_path)

    result = features[["Sample"]].copy()
    for name, (feature_column, _) in MODEL_SPECS.items():
        result[f"Pred_{name}"] = models[name].predict(features[feature_column].to_numpy(dtype=float))

    result["Pred_PGA34A"] = (result["Pred_PGA34"] - result["Pred_PGA34B"]).clip(lower=0)
    result = result[["Sample", "Pred_PGA", "Pred_PGA34", "Pred_PGA34A", "Pred_PGA34B", "Pred_PGA5"]]
    result.to_csv(output, sep="\t", index=False)
    print(f"Predictions: {output}")


def build_parser():
    parser = argparse.ArgumentParser(description="PGA copy-number feature extraction and constrained-GMM genotyping.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_extract = subparsers.add_parser("extract", help="Extract continuous CN features.")
    p_extract.add_argument("--sample-list", required=True)
    p_extract.add_argument("--reference", required=True, help="Pseudo-PGA FASTA.")
    p_extract.add_argument("--manifest", required=True)
    p_extract.add_argument("--out", required=True)
    p_extract.add_argument("--threads", type=int, default=4)

    p_fit = subparsers.add_parser("fit", help="Fit constrained GMMs using assembly-resolved truth labels.")
    p_fit.add_argument("--features", required=True, help="Feature table for labeled samples.")
    p_fit.add_argument("--truth", required=True)
    p_fit.add_argument("--unlabeled-features", default=None, help="Optional target-cohort feature table used for unsupervised mixture fitting.")
    p_fit.add_argument("--model-out", required=True)

    p_predict = subparsers.add_parser("predict", help="Predict integer copy numbers with fitted models.")
    p_predict.add_argument("--features", required=True)
    p_predict.add_argument("--model", required=True)
    p_predict.add_argument("--out", required=True)

    return parser


def main():
    args = build_parser().parse_args()
    if args.command == "extract":
        extract_features(args.sample_list, args.reference, args.manifest, args.out, args.threads)
    elif args.command == "fit":
        fit_models(args.features, args.truth, args.model_out, args.unlabeled_features)
    elif args.command == "predict":
        predict(args.features, args.model, args.out)


if __name__ == "__main__":
    main()
