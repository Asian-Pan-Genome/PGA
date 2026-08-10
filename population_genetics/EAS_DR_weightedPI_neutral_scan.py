#!/usr/bin/env python3
"""Neutral scan for the EAS DR weighted mean windowed nucleotide diversity."""

import argparse
import multiprocessing as mp
import time
from dataclasses import dataclass

import matplotlib.pyplot as plt
import msprime
import numpy as np
import pandas as pd


DEFAULT_REGION = "chr11:61090000-61506000"
DEFAULT_SEARCH_LENGTH_BP = 906_000
DEFAULT_WINDOW_BP = 10_000
DEFAULT_STEP_BP = 1_000
DEFAULT_N_DIPLOID = 223
DEFAULT_MUTATION_RATE = 1.25e-8
DEFAULT_RECOMBINATION_RATE = 1e-8
DEFAULT_MAF = 0.05
DEFAULT_REPLICATES = 10_000


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


def parse_region(region_str: str) -> Region:
    chrom, coords = region_str.split(":")
    start_s, end_s = coords.replace(",", "").split("-")
    region = Region(chrom=chrom, start=int(start_s), end=int(end_s))
    if region.end <= region.start:
        raise ValueError(f"Invalid region: {region_str}")
    return region


def load_relate_demography(csv_path: str, population: str) -> msprime.Demography:
    """Load a Relate-inferred population-size history from long or wide CSV."""
    df = pd.read_csv(csv_path)
    df.columns = [c.strip() for c in df.columns]

    if population in df.columns and "gens_ago" in df.columns:
        history = df[["gens_ago", population]].rename(
            columns={population: "population_size"}
        )
    elif {"region", "gens_ago", "population_size"}.issubset(df.columns):
        history = df.loc[
            df["region"] == population, ["gens_ago", "population_size"]
        ].copy()
    elif "region" in df.columns and "gens_ago" in df.columns:
        size_cols = [c for c in df.columns if c not in {"region", "gens_ago"}]
        if len(size_cols) != 1:
            raise ValueError(
                "Cannot infer the population-size column from the demography CSV."
            )
        history = df.loc[
            df["region"] == population, ["gens_ago", size_cols[0]]
        ].rename(columns={size_cols[0]: "population_size"})
    else:
        raise ValueError(
            "Expected either wide format (gens_ago plus population columns) or "
            "long format (region, gens_ago, population_size)."
        )

    if history.empty:
        raise ValueError(f"Population {population!r} not found in {csv_path}")

    history["gens_ago"] = pd.to_numeric(history["gens_ago"], errors="raise")
    history["population_size"] = pd.to_numeric(
        history["population_size"], errors="raise"
    )
    history = (
        history.groupby("gens_ago", as_index=False)["population_size"]
        .mean()
        .sort_values("gens_ago")
    )

    demography = msprime.Demography()
    demography.add_population(
        name=population,
        initial_size=float(history.iloc[0]["population_size"]),
    )
    for _, row in history.iloc[1:].iterrows():
        demography.add_population_parameters_change(
            time=float(row["gens_ago"]),
            initial_size=float(row["population_size"]),
            population=population,
        )
    demography.sort_events()
    return demography


def read_windowed_pi(pi_path: str) -> pd.DataFrame:
    df = pd.read_csv(pi_path, sep="\t")
    required = {"CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "PI"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{pi_path} missing columns: {sorted(missing)}")

    for col in ["BIN_START", "BIN_END", "N_VARIANTS", "PI"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(
        subset=["BIN_START", "BIN_END", "N_VARIANTS", "PI"]
    ).copy()
    df["BIN_START"] = df["BIN_START"].astype(int)
    df["BIN_END"] = df["BIN_END"].astype(int)

    # VCFtools .windowed.pi coordinates are 1-based and inclusive.
    df["WINDOW_LEN"] = df["BIN_END"] - df["BIN_START"] + 1
    df["CENTER"] = (df["BIN_START"] + df["BIN_END"]) / 2.0
    return df


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    keep = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(keep):
        return np.nan
    return float(np.average(values[keep], weights=weights[keep]))


def observed_weighted_pi(
    pi_path: str,
    region: Region,
    window_bp: int,
    out_prefix: str,
) -> dict:
    df = read_windowed_pi(pi_path)
    sub = df.loc[
        (df["CHROM"].astype(str) == region.chrom)
        & (df["WINDOW_LEN"] == window_bp)
        & (df["CENTER"] >= region.start)
        & (df["CENTER"] < region.end)
    ].copy()

    if sub.empty:
        raise ValueError(
            "No centered windows found in the observed PI file. "
            f"Check region={region.chrom}:{region.start}-{region.end} and {pi_path}."
        )

    sub.to_csv(
        f"{out_prefix}.observed_centered_windows.tsv", sep="\t", index=False
    )

    pi = sub["PI"].to_numpy(dtype=float)
    n_variants = sub["N_VARIANTS"].to_numpy(dtype=float)
    observed = weighted_mean(pi, n_variants)
    if not np.isfinite(observed):
        raise ValueError("Observed weighted mean PI is not finite.")

    return {
        "observed_weighted_mean_pi": observed,
        "observed_region": f"{region.chrom}:{region.start}-{region.end}",
        "observed_region_length_bp": region.length,
        "observed_n_windows": int(len(sub)),
        "observed_sum_n_variants": float(np.nansum(n_variants)),
        "observed_mean_n_variants": float(np.nanmean(n_variants)),
        "observed_min_n_variants": float(np.nanmin(n_variants)),
        "observed_max_n_variants": float(np.nanmax(n_variants)),
    }


def make_padded_search_design(
    search_length_bp: int,
    window_bp: int,
    step_bp: int,
):
    """Create sliding windows around a padded neutral search interval."""
    sequence_length = search_length_bp + window_bp
    search_start = window_bp / 2.0
    starts = np.arange(0, sequence_length - window_bp + 1, step_bp, dtype=float)
    ends = starts + window_bp
    centers = (starts + ends) / 2.0
    return sequence_length, search_start, starts, ends, centers


def extract_maf_filtered_biallelic_sites(mts, maf: float):
    n_haploid = mts.num_samples
    positions = []
    pi_contrib = []

    for var in mts.variants():
        if len(var.alleles) != 2:
            continue

        genotypes = var.genotypes
        if np.any((genotypes != 0) & (genotypes != 1)):
            continue

        alt_count = int(np.sum(genotypes))
        if alt_count <= 0 or alt_count >= n_haploid:
            continue

        site_maf = min(
            alt_count / n_haploid,
            1.0 - alt_count / n_haploid,
        )
        if site_maf < maf:
            continue

        pi_site = (
            2.0
            * alt_count
            * (n_haploid - alt_count)
            / (n_haploid * (n_haploid - 1.0))
        )
        positions.append(float(var.site.position))
        pi_contrib.append(pi_site)

    if not positions:
        return np.array([], dtype=float), np.array([], dtype=float)

    positions = np.asarray(positions, dtype=float)
    pi_contrib = np.asarray(pi_contrib, dtype=float)
    order = np.argsort(positions)
    return positions[order], pi_contrib[order]


def windowed_pi_from_sites(
    positions: np.ndarray,
    pi_contrib: np.ndarray,
    starts: np.ndarray,
    ends: np.ndarray,
    window_bp: int,
):
    n_windows = len(starts)
    if len(positions) == 0:
        return np.zeros(n_windows, dtype=float), np.zeros(n_windows, dtype=int)

    left = np.searchsorted(positions, starts, side="left")
    right = np.searchsorted(positions, ends, side="left")

    prefix_pi = np.concatenate(([0.0], np.cumsum(pi_contrib)))
    pi_sum = prefix_pi[right] - prefix_pi[left]
    n_variants = right - left
    pi_window = pi_sum / float(window_bp)
    return pi_window.astype(float), n_variants.astype(int)


def scan_min_weighted_pi(
    pi_window: np.ndarray,
    n_variants: np.ndarray,
    centers: np.ndarray,
    search_start: float,
    search_length_bp: int,
    target_length_bp: int,
    step_bp: int,
):
    max_offset = search_length_bp - target_length_bp
    if max_offset < 0:
        raise ValueError("search_length_bp must be >= target_length_bp")

    offsets = np.arange(0, max_offset + 1, step_bp, dtype=float)
    numerator = np.where(
        np.isfinite(pi_window), pi_window * n_variants, 0.0
    )
    denominator = np.where(np.isfinite(pi_window), n_variants, 0).astype(float)

    prefix_num = np.concatenate(([0.0], np.cumsum(numerator)))
    prefix_den = np.concatenate(([0.0], np.cumsum(denominator)))
    prefix_nvar = np.concatenate(([0.0], np.cumsum(n_variants.astype(float))))

    best_stat = np.nan
    best_offset = np.nan
    best_n_windows = 0
    best_sum_n_variants = np.nan

    for offset in offsets:
        interval_start = search_start + offset
        interval_end = interval_start + target_length_bp
        left = np.searchsorted(centers, interval_start, side="left")
        right = np.searchsorted(centers, interval_end, side="left")

        weight_sum = prefix_den[right] - prefix_den[left]
        if weight_sum <= 0:
            continue

        statistic = (prefix_num[right] - prefix_num[left]) / weight_sum
        if not np.isfinite(best_stat) or statistic < best_stat:
            best_stat = float(statistic)
            best_offset = float(offset)
            best_n_windows = int(right - left)
            best_sum_n_variants = float(prefix_nvar[right] - prefix_nvar[left])

    return {
        "sim_min_weighted_mean_pi": best_stat,
        "sim_best_interval_start_offset": best_offset,
        "sim_best_interval_end_offset": best_offset + target_length_bp,
        "sim_best_interval_n_windows": best_n_windows,
        "sim_best_interval_sum_n_variants": best_sum_n_variants,
        "sim_n_candidate_intervals": int(len(offsets)),
    }


def run_batch(
    seed: int,
    n_reps: int,
    demography_csv: str,
    population: str,
    sequence_length: int,
    search_start: float,
    search_length_bp: int,
    target_length_bp: int,
    starts: np.ndarray,
    ends: np.ndarray,
    centers: np.ndarray,
    n_diploid: int,
    mutation_rate: float,
    recombination_rate: float,
    maf: float,
    window_bp: int,
    step_bp: int,
):
    rng = np.random.default_rng(seed)
    demography = load_relate_demography(demography_csv, population)
    rows = []

    replicates = msprime.sim_ancestry(
        samples=[msprime.SampleSet(n_diploid, population=population, ploidy=2)],
        demography=demography,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        num_replicates=n_reps,
        random_seed=int(seed),
    )

    for local_rep_idx, ts in enumerate(replicates):
        mutation_seed = int(rng.integers(1, 2**32 - 1))
        mts = msprime.sim_mutations(
            ts,
            rate=mutation_rate,
            random_seed=mutation_seed,
        )

        positions, pi_contrib = extract_maf_filtered_biallelic_sites(mts, maf)
        pi_window, n_variants = windowed_pi_from_sites(
            positions,
            pi_contrib,
            starts,
            ends,
            window_bp,
        )
        scan = scan_min_weighted_pi(
            pi_window,
            n_variants,
            centers,
            search_start,
            search_length_bp,
            target_length_bp,
            step_bp,
        )

        rows.append(
            {
                "seed": seed,
                "local_rep_idx": local_rep_idx,
                "n_maf_filtered_biallelic_sites": int(len(positions)),
                **scan,
            }
        )

    return rows


def empirical_p_leq(sim_values: np.ndarray, observed: float) -> dict:
    valid = np.isfinite(sim_values)
    values = sim_values[valid]
    if len(values) == 0:
        raise ValueError("No valid simulated values.")

    count = int(np.sum(values <= observed))
    p_value = count / len(values)
    p_label = f"< {1.0 / len(values):.6g}" if count == 0 else f"{p_value:.6g}"

    return {
        "p_raw": p_value,
        "p_label": p_label,
        "n_valid_simulations": int(len(values)),
        "n_invalid_simulations": int((~valid).sum()),
        "n_sim_le_observed": count,
        "p_upper_bound_if_zero": float(1.0 / len(values)),
    }


def plot_histogram(
    sim_values: np.ndarray,
    observed: float,
    p_label: str,
    out_pdf: str,
):
    values = sim_values[np.isfinite(sim_values)]

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.hist(values, bins=50)
    ax.axvline(observed, linestyle="--", linewidth=2)
    ax.set_xlabel("Simulated minimum weighted mean windowed π")
    ax.set_ylabel("Replicates")
    ax.text(
        0.04,
        0.96,
        (
            f"Observed = {observed:.6g}\n"
            f"Empirical lower-tail P = {p_label}\n"
            f"Valid simulations = {len(values)}"
        ),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=10,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_pdf, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Neutral scan for weighted mean windowed PI in the EAS DR."
    )
    parser.add_argument("--pi-file", required=True)
    parser.add_argument("--demography-csv", default="./pop_size.Relate.csv")
    parser.add_argument("--target-region", default="EAS")
    parser.add_argument("--region", default=DEFAULT_REGION)
    parser.add_argument("--search-length-bp", type=int, default=DEFAULT_SEARCH_LENGTH_BP)
    parser.add_argument("--window-bp", type=int, default=DEFAULT_WINDOW_BP)
    parser.add_argument("--step-bp", type=int, default=DEFAULT_STEP_BP)
    parser.add_argument("--n-diploid", type=int, default=DEFAULT_N_DIPLOID)
    parser.add_argument("--mutation-rate", type=float, default=DEFAULT_MUTATION_RATE)
    parser.add_argument(
        "--recombination-rate", type=float, default=DEFAULT_RECOMBINATION_RATE
    )
    parser.add_argument("--maf", type=float, default=DEFAULT_MAF)
    parser.add_argument("--replicates", type=int, default=DEFAULT_REPLICATES)
    parser.add_argument("--processes", type=int, default=max(1, mp.cpu_count() - 2))
    parser.add_argument("--out-prefix", default="EAS.DR.weightedMeanPI.maf5.scanSVR")
    args = parser.parse_args()

    region = parse_region(args.region)
    if args.search_length_bp < region.length:
        parser.error("--search-length-bp must be >= the target region length")
    if args.window_bp <= 0 or args.step_bp <= 0:
        parser.error("--window-bp and --step-bp must be positive")
    if args.n_diploid <= 0 or args.replicates <= 0 or args.processes <= 0:
        parser.error("sample size, replicates, and processes must be positive")
    if not 0 <= args.maf <= 0.5:
        parser.error("--maf must be between 0 and 0.5")

    sequence_length, search_start, starts, ends, centers = make_padded_search_design(
        args.search_length_bp,
        args.window_bp,
        args.step_bp,
    )

    print("=" * 80)
    print("EAS DR neutral scan for weighted mean windowed PI")
    print("=" * 80)
    print(f"Observed PI file:       {args.pi_file}")
    print(f"Demography CSV:         {args.demography_csv}")
    print(f"Population:             {args.target_region}")
    print(f"Observed region:        {region.chrom}:{region.start}-{region.end}")
    print(f"Target length:          {region.length} bp")
    print(f"Search length:          {args.search_length_bp} bp")
    print(f"Padded sequence length: {sequence_length} bp")
    print(f"Window / step:          {args.window_bp} / {args.step_bp} bp")
    print(f"MAF threshold:          {args.maf}")
    print(f"Sample size:            {args.n_diploid} diploid")
    print(f"Replicates:             {args.replicates}")
    print("=" * 80)

    observed_info = observed_weighted_pi(
        args.pi_file,
        region,
        args.window_bp,
        args.out_prefix,
    )
    observed = observed_info["observed_weighted_mean_pi"]

    print("\nObserved statistic")
    for key, value in observed_info.items():
        print(f"  {key}: {value}")

    start_time = time.time()
    n_processes = min(args.processes, args.replicates)
    base_reps = args.replicates // n_processes

    tasks = []
    for i in range(n_processes):
        n_reps = base_reps + (1 if i < args.replicates % n_processes else 0)
        seed = 42 + i * 10_000
        tasks.append(
            (
                seed,
                n_reps,
                args.demography_csv,
                args.target_region,
                sequence_length,
                search_start,
                args.search_length_bp,
                region.length,
                starts,
                ends,
                centers,
                args.n_diploid,
                args.mutation_rate,
                args.recombination_rate,
                args.maf,
                args.window_bp,
                args.step_bp,
            )
        )

    if len(tasks) == 1:
        batches = [run_batch(*tasks[0])]
    else:
        with mp.Pool(processes=len(tasks)) as pool:
            batches = pool.starmap(run_batch, tasks)

    rows = [row for batch in batches for row in batch]
    sim_df = pd.DataFrame(rows)
    sim_df.insert(0, "replicate", np.arange(len(sim_df), dtype=int))

    simulated_tsv = f"{args.out_prefix}.simulated.tsv"
    sim_df.to_csv(simulated_tsv, sep="\t", index=False)

    sim_values = sim_df["sim_min_weighted_mean_pi"].to_numpy(dtype=float)
    p_info = empirical_p_leq(sim_values, observed)
    elapsed = time.time() - start_time

    summary = {
        **observed_info,
        **p_info,
        "region": f"{region.chrom}:{region.start}-{region.end}",
        "target_length_bp": region.length,
        "search_length_bp": args.search_length_bp,
        "padded_sequence_length_bp": sequence_length,
        "window_bp": args.window_bp,
        "step_bp": args.step_bp,
        "n_diploid": args.n_diploid,
        "mutation_rate": args.mutation_rate,
        "recombination_rate": args.recombination_rate,
        "maf": args.maf,
        "requested_replicates": args.replicates,
        "elapsed_seconds": elapsed,
        "sim_mean_min_weighted_mean_pi": float(np.nanmean(sim_values)),
        "sim_sd_min_weighted_mean_pi": float(np.nanstd(sim_values, ddof=1)),
        "sim_min_min_weighted_mean_pi": float(np.nanmin(sim_values)),
        "sim_max_min_weighted_mean_pi": float(np.nanmax(sim_values)),
    }

    summary_tsv = f"{args.out_prefix}.summary.tsv"
    pd.DataFrame([summary]).to_csv(summary_tsv, sep="\t", index=False)

    histogram_pdf = f"{args.out_prefix}.histogram.pdf"
    plot_histogram(sim_values, observed, p_info["p_label"], histogram_pdf)

    print("\nSimulation summary")
    for key, value in summary.items():
        print(f"  {key}: {value}")

    print("\nOutput files")
    print(f"  Observed centered windows: {args.out_prefix}.observed_centered_windows.tsv")
    print(f"  Simulated values:          {simulated_tsv}")
    print(f"  Summary:                   {summary_tsv}")
    print(f"  Histogram:                 {histogram_pdf}")

    if p_info["n_sim_le_observed"] == 0:
        print(
            "\nReporting note: empirical count is 0. "
            f"Report as P < {p_info['p_upper_bound_if_zero']:.6g}, not P = 0."
        )


if __name__ == "__main__":
    main()
