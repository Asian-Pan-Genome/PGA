#!/usr/bin/env python3

import argparse
import itertools
import multiprocessing as mp
import re
from functools import partial

import numpy as np
import pandas as pd


POPULATIONS = [
    "AFR",
    "AFR-E&S",
    "AFR-W",
    "AFR-NE",
    "ARB",
    "EUR",
    "WAS",
    "CAS",
    "CSA",
    "SAS",
    "EAS",
    "SEA",
    "AMR",
]

AFR_SUBGROUPS = ["AFR-E&S", "AFR-NE", "AFR-W"]

DEFAULT_EXCLUDE_PATTERN = (
    r"^(?:LINC|USP17L|NBPF|ZNF|FAM90A|IG[HKL]|TR[ABDG]|PCDH|H[234][ABC]|"
    r"CT4[75]|GOLGA|KRT|GAGE|OR|TBC1D3|NPIP)"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Rank pairwise population differentiation in gene-family copy number "
            "using Vst x absolute mean CN difference."
        )
    )
    parser.add_argument("--cn-matrix", required=True, help="Genome-wide family CN matrix")
    parser.add_argument(
        "--sample-table",
        required=True,
        help="TSV defining the PGA analysis cohort; requires Sample and Hap columns",
    )
    parser.add_argument(
        "--population-table",
        required=True,
        help="TSV with Sample and New_superpop columns",
    )
    parser.add_argument(
        "--target-gene",
        default="PGA3",
        help="Row identifier used for the target family in the CN matrix (default: PGA3)",
    )
    parser.add_argument(
        "--target-label",
        default="PGA",
        help="Human-readable target label written to outputs (default: PGA)",
    )
    parser.add_argument(
        "--target-cn-table",
        default=None,
        help=(
            "Optional paralog-resolved CN table used to replace --target-gene. "
            "Requires Sample and Hap columns."
        ),
    )
    parser.add_argument(
        "--target-cn-columns",
        nargs="+",
        default=None,
        help=(
            "Columns from --target-cn-table to sum for the replacement target CN, "
            "for example PGA34A1 PGA34A2."
        ),
    )
    parser.add_argument(
        "--exclude-pattern",
        default=DEFAULT_EXCLUDE_PATTERN,
        help="Regex applied to representative gene names before ranking",
    )
    parser.add_argument(
        "--full-pair",
        nargs=2,
        metavar=("POP1", "POP2"),
        default=["EAS", "EUR"],
        help="Population pair for which the full genome-wide ranking table is written",
    )
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--out-prefix", required=True)
    return parser.parse_args()


def normalize_sample_names(df):
    df = df.copy()
    df["Sample"] = df["Sample"].astype(str).str.replace("v1.1", "", regex=False)
    return df


def load_population_samples(sample_table_path, population_table_path):
    samples = pd.read_csv(sample_table_path, sep="\t", header=0)
    required = {"Sample", "Hap"}
    missing = required - set(samples.columns)
    if missing:
        raise ValueError(f"{sample_table_path} missing columns: {sorted(missing)}")

    if "Superpop" in samples.columns:
        samples = samples.dropna(subset=["Superpop"]).copy()

    population_table = pd.read_csv(population_table_path, sep="\t", header=0)
    required = {"Sample", "New_superpop"}
    missing = required - set(population_table.columns)
    if missing:
        raise ValueError(f"{population_table_path} missing columns: {sorted(missing)}")

    samples = pd.merge(samples, population_table, on="Sample", how="left")
    samples = samples.dropna(subset=["New_superpop"]).copy()
    samples = normalize_sample_names(samples)
    samples["sample_hap"] = samples["Sample"] + "." + samples["Hap"].astype(str)

    aggregate = samples[samples["New_superpop"].isin(AFR_SUBGROUPS)].copy()
    aggregate["New_superpop"] = "AFR"
    samples = pd.concat([samples, aggregate], ignore_index=True)

    pop_samples = {}
    for pop in POPULATIONS:
        pop_samples[pop] = (
            samples.loc[samples["New_superpop"] == pop, "sample_hap"]
            .drop_duplicates()
            .tolist()
        )

    empty = [pop for pop, ids in pop_samples.items() if not ids]
    if empty:
        raise ValueError(f"No haplotypes found for populations: {', '.join(empty)}")

    return pop_samples


def replace_target_cn(cnv_table, args):
    if args.target_cn_table is None:
        if args.target_cn_columns is not None:
            raise ValueError("--target-cn-columns requires --target-cn-table")
        return cnv_table

    if not args.target_cn_columns:
        raise ValueError("--target-cn-table requires --target-cn-columns")

    target = pd.read_csv(args.target_cn_table, sep="\t", header=0)
    required = {"Sample", "Hap", *args.target_cn_columns}
    missing = required - set(target.columns)
    if missing:
        raise ValueError(f"{args.target_cn_table} missing columns: {sorted(missing)}")

    target = normalize_sample_names(target)
    target["sample_hap"] = target["Sample"] + "." + target["Hap"].astype(str)
    target[args.target_label] = target[args.target_cn_columns].sum(axis=1)

    if target["sample_hap"].duplicated().any():
        raise ValueError("Duplicate sample_hap values in target CN table")

    replacement = target.set_index("sample_hap")[args.target_label]

    cnv_table = cnv_table[cnv_table["gene"] != args.target_gene].copy()
    row = {"gene": args.target_gene, "gene_family": args.target_label}
    for column in cnv_table.columns:
        if column in {"gene", "gene_family"}:
            continue
        row[column] = replacement.get(column, np.nan)

    return pd.concat([cnv_table, pd.DataFrame([row])], ignore_index=True)


def compute_vst(values_g1, values_g2):
    if len(values_g1) < 2 or len(values_g2) < 2:
        return np.nan

    combined = np.concatenate([values_g1, values_g2])
    if np.max(combined) - np.min(combined) <= 0:
        return np.nan

    v1 = np.var(values_g1, ddof=1)
    v2 = np.var(values_g2, ddof=1)
    vtotal = np.var(combined, ddof=1)
    if vtotal == 0:
        return np.nan

    vwithin = (
        v1 * (len(values_g1) - 1) + v2 * (len(values_g2) - 1)
    ) / ((len(values_g1) - 1) + (len(values_g2) - 1))

    vst = (vtotal - vwithin) / vtotal
    return np.nan if vst < 0 else float(vst)


def rank_pair(pop1, pop2, cnv_table, pop_samples, target_gene, full_pair):
    columns1 = pop_samples[pop1]
    columns2 = pop_samples[pop2]

    results = []
    for _, row in cnv_table.iterrows():
        values1 = row[columns1].to_numpy(dtype=float)
        values2 = row[columns2].to_numpy(dtype=float)

        vst = compute_vst(values1, values2)
        mean_diff = abs(float(np.mean(values1)) - float(np.mean(values2)))
        results.append((row["gene"], vst, mean_diff))

    ranked = pd.DataFrame(results, columns=["Gene", "Vst", "Average_Diff"])
    ranked = ranked.dropna(subset=["Vst", "Average_Diff"])
    ranked = ranked[ranked["Average_Diff"] != 0].copy()
    ranked["Combined_Metric"] = ranked["Vst"] * ranked["Average_Diff"]
    ranked = ranked.sort_values("Combined_Metric", ascending=False).reset_index(drop=True)
    ranked["Rank"] = (ranked.index + 1) / len(ranked)
    ranked["Rank_Percentile"] = ranked["Rank"] * 100

    target = ranked[ranked["Gene"] == target_gene]
    if target.empty:
        summary = {
            "Population_1": pop1,
            "Population_2": pop2,
            "Status": "target_not_ranked",
            "Vst": np.nan,
            "Average_Diff": np.nan,
            "Combined_Metric": np.nan,
            "Rank": np.nan,
            "Rank_Percentile": np.nan,
            "Gene_Index": np.nan,
            "Total_Genes": len(ranked),
        }
    else:
        index = int(target.index[0])
        row = target.iloc[0]
        summary = {
            "Population_1": pop1,
            "Population_2": pop2,
            "Status": "assessed",
            "Vst": float(row["Vst"]),
            "Average_Diff": float(row["Average_Diff"]),
            "Combined_Metric": float(row["Combined_Metric"]),
            "Rank": float(row["Rank"]),
            "Rank_Percentile": float(row["Rank_Percentile"]),
            "Gene_Index": index + 1,
            "Total_Genes": len(ranked),
        }

    requested_pair = frozenset(full_pair)
    pair_table = ranked if frozenset((pop1, pop2)) == requested_pair else None
    return summary, pair_table


def is_overlapping_afr_comparison(pop1, pop2):
    return "AFR" in {pop1, pop2} and bool({pop1, pop2} & set(AFR_SUBGROUPS))


def main():
    args = parse_args()

    cnv_table = pd.read_csv(args.cn_matrix, sep="\t", header=0)
    if "gene" not in cnv_table.columns:
        raise ValueError(f"{args.cn_matrix} must contain a gene column")
    if "gene_family" not in cnv_table.columns:
        cnv_table.insert(1, "gene_family", cnv_table["gene"])

    cnv_table = cnv_table[
        ~cnv_table["gene"].astype(str).str.contains(args.exclude_pattern, regex=True)
    ].copy()
    cnv_table = replace_target_cn(cnv_table, args)

    if cnv_table["gene"].duplicated().any():
        duplicated = cnv_table.loc[cnv_table["gene"].duplicated(), "gene"].tolist()
        raise ValueError(f"Duplicated gene-family identifiers: {duplicated[:10]}")

    pop_samples = load_population_samples(args.sample_table, args.population_table)
    selected_samples = sorted({sample for ids in pop_samples.values() for sample in ids})

    missing_columns = sorted(set(selected_samples) - set(cnv_table.columns))
    if missing_columns:
        raise ValueError(
            "CN matrix is missing selected haplotype columns, for example: "
            + ", ".join(missing_columns[:10])
        )

    numeric = cnv_table[selected_samples].apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        n_missing = int(numeric.isna().sum().sum())
        raise ValueError(
            f"Found {n_missing} missing/non-numeric CN values in the selected cohort. "
            "The original ranking assumes a complete genome-wide CN matrix; resolve "
            "missing values before ranking rather than changing the denominator silently."
        )
    cnv_table[selected_samples] = numeric

    if args.target_gene not in set(cnv_table["gene"]):
        raise ValueError(f"Target row not found after filtering: {args.target_gene}")

    full_pair = tuple(args.full_pair)
    if any(pop not in POPULATIONS for pop in full_pair):
        raise ValueError(f"Unknown population in --full-pair: {full_pair}")

    assessed_tasks = []
    summaries = []
    for pop1, pop2 in itertools.combinations(POPULATIONS, 2):
        if is_overlapping_afr_comparison(pop1, pop2):
            summaries.append(
                {
                    "Population_1": pop1,
                    "Population_2": pop2,
                    "Status": "not_assessed_overlapping_AFR",
                    "Vst": np.nan,
                    "Average_Diff": np.nan,
                    "Combined_Metric": np.nan,
                    "Rank": np.nan,
                    "Rank_Percentile": np.nan,
                    "Gene_Index": np.nan,
                    "Total_Genes": np.nan,
                }
            )
        else:
            assessed_tasks.append((pop1, pop2))

    worker = partial(
        rank_pair,
        cnv_table=cnv_table,
        pop_samples=pop_samples,
        target_gene=args.target_gene,
        full_pair=full_pair,
    )

    if args.processes == 1:
        results = [worker(*task) for task in assessed_tasks]
    else:
        with mp.Pool(processes=args.processes) as pool:
            results = pool.starmap(worker, assessed_tasks)

    full_pair_table = None
    for summary, pair_table in results:
        summaries.append(summary)
        if pair_table is not None:
            full_pair_table = pair_table

    summary_df = pd.DataFrame(summaries)
    summary_df.insert(0, "Target", args.target_label)
    summary_df.to_csv(f"{args.out_prefix}.pairwise_rank.tsv", sep="\t", index=False)

    if full_pair_table is None:
        raise RuntimeError(f"Requested full pair was not assessed: {full_pair}")

    full_pair_table.insert(0, "Population_1", full_pair[0])
    full_pair_table.insert(1, "Population_2", full_pair[1])
    full_pair_table.to_csv(
        f"{args.out_prefix}.{full_pair[0]}_vs_{full_pair[1]}.genomewide.tsv",
        sep="\t",
        index=False,
    )

    n_assessed = int((summary_df["Status"] == "assessed").sum())
    n_not_assessed = int(summary_df["Status"].str.startswith("not_assessed").sum())
    print(f"Filtered background families: {len(cnv_table)}")
    print(f"Assessed target comparisons: {n_assessed}")
    print(f"Not assessed overlapping AFR comparisons: {n_not_assessed}")
    print(f"Wrote {args.out_prefix}.pairwise_rank.tsv")
    print(
        f"Wrote {args.out_prefix}.{full_pair[0]}_vs_{full_pair[1]}.genomewide.tsv"
    )


if __name__ == "__main__":
    main()
