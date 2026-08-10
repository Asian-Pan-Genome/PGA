#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AADR v66 2M ancient SNP-panel windowed diversity analysis for the PGA desert region.

Purpose
-------
This script computes time-resolved, windowed SNP-panel diversity statistics for:

  East Eurasia (80E core) across the PGA DR
  West Eurasia across the left portion of the PGA DR (DR-L)

using chr11-wide AADR EIGENSTRAT genotypes.

Major choices implemented from the project discussion
-----------------------------------------------------
1. East Eurasia = longitude >= 80E and latitude >= -10.
2. West Eurasia = Europe + West Asia up to 80E + Greenland.
3. East target = DR span chr11:60,857,472-61,273,472, but observed DR
   statistics exclude windows centered inside the PGA/SVR gene cluster
   chr11:60,970,987-61,018,916.
4. West target = DR-L chr11:60,857,472-60,970,987.
5. chr11 empirical background uses same-length candidate intervals across the
   entire chromosome 11, without excluding PGA/SVR.
6. Windowed statistics use 10 kb windows with 1 kb step and center-based
   inclusion.
7. No MAF filtering is applied. Monomorphic callable SNPs are retained.
   Windows with no callable SNPs are treated as NA, not zero.
8. Metrics:
     - weighted_windowed_pi, weighted by n_callable_sites per window
     - mean_windowed_tajimaD_like
9. Bootstrap CIs are individual-level bootstrap CIs on the observed target only.
   Background candidate intervals are not bootstrapped.
10. Rolling analysis is the unbinned-equivalent trend:
      rolling width = 1000 years, step = 100 years by default.
11. Binned sensitivity analysis uses 1-, 2-, and 3-kyr bins by default;
    sparsely sampled oldest intervals are pooled into a terminal bin.

Coordinate convention
---------------------
All regions are interpreted as half-open intervals [start, end). Window centers
are included when start <= center < end.

Input expectations
------------------
--geno/--snp/--ind should be chr11-wide EIGENSTRAT output for valid samples.
--metadata should contain at least:
    Genetic ID, Individual ID, Latitude, Longitude, Date_BP
--anno is used only to flag/remove AADR Ignore_/outlier samples.

Interpretation warning
----------------------
AADR ancient genotypes are mostly pseudo-haploid and the 2M panel is an
ascertained SNP panel. Tajima's D-like values should be interpreted as
panel-based allele-frequency-spectrum summaries, not classical neutrality tests.
"""

import argparse
import gzip
import hashlib
import math
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import linregress
from statsmodels.nonparametric.smoothers_lowess import lowess as sm_lowess


# -----------------------------
# Constants and palettes
# -----------------------------

CHR11_HG19_END = 135_086_622  # inclusive hg19 chr11 length

DR_START = 60_857_472
DR_END = 61_273_472
DRL_START = 60_857_472
DRL_END = 60_970_987
SVR_START = 60_970_987
SVR_END = 61_018_916

PALETTE = {
    'AFR':     '#E64B35',
    'AFR-E&S': '#FFD5C2',
    'AFR-W':   '#741B11',
    'AFR-NE':  '#FF9900',
    'ARB':     '#4DBBD5',
    'EUR':     '#3C5488',
    'WAS':     '#8491B4',
    'CAS':     '#D9CCF2',
    'CSA':     '#8467BD',
    'SAS':     '#311B92',
    'EAS':     '#00A087',
    'SEA':     '#91D1C2',
    'AMR':     '#7F7F7F',
}

GROUP_STYLE = {
    "East Eurasia (DR)": {
        "color": PALETTE["EAS"],
        "direction": "lower",
        "region_name": "DR",
    },
    "West Eurasia (DR-L)": {
        "color": PALETTE["EUR"],
        "direction": "upper",
        "region_name": "DR-L",
    },
}

METRIC_LABELS = {
    "weighted_windowed_pi": "Weighted windowed $\\pi$",
    "mean_windowed_tajimaD_like": "Mean windowed Tajima's $D$-like",
}

# Labels used in plots; region names are omitted from legends/text to keep figures clean.
GROUP_DISPLAY_LABELS = {
    "East Eurasia (DR)": "East Eurasia",
    "West Eurasia (DR-L)": "West Eurasia",
}

METRIC_SIMPLE_TITLES = {
    "weighted_windowed_pi": "Weighted windowed $\\pi$",
    "mean_windowed_tajimaD_like": "Mean windowed Tajima's $D$-like",
}


@dataclass(frozen=True)
class Region:
    name: str
    start: int
    end: int
    background_length: int
    exclude_svr_for_observed: bool = False

    @property
    def length(self) -> int:
        return self.end - self.start


REGIONS = {
    "DR": Region(
        name="DR",
        start=DR_START,
        end=DR_END,
        background_length=DR_END - DR_START,
        exclude_svr_for_observed=True,
    ),
    "DR-L": Region(
        name="DR-L",
        start=DRL_START,
        end=DRL_END,
        background_length=DRL_END - DRL_START,
        exclude_svr_for_observed=False,
    ),
}


# -----------------------------
# Utility functions
# -----------------------------

def configure_matplotlib():
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["mathtext.fontset"] = "custom"
    plt.rcParams["mathtext.rm"] = "Arial"
    plt.rcParams["mathtext.it"] = "Arial:italic"
    plt.rcParams["mathtext.bf"] = "Arial:bold"
    plt.rcParams["axes.unicode_minus"] = False


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def norm_col(x):
    return " ".join(str(x).replace("\xa0", " ").split()).lower()


def find_col(cols, *, startswith=None, exact=None, required=True):
    hits = []
    for c in cols:
        nc = norm_col(c)
        ok = True
        if startswith is not None:
            ok &= nc.startswith(startswith.lower())
        if exact is not None:
            ok &= nc == exact.lower()
        if ok:
            hits.append(c)
    if len(hits) == 1:
        return hits[0]
    if required:
        raise KeyError(f"Expected one matching column, got {len(hits)}: {hits}")
    return None


def safe_mean(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan
    return float(np.mean(x))


def safe_weighted_mean(values, weights):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    ok = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if ok.sum() == 0:
        return np.nan
    return float(np.average(values[ok], weights=weights[ok]))


def empirical_percentiles(values, observed):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0 or not np.isfinite(observed):
        return {
            "n_candidate_regions": int(len(values)),
            "percentile_lower_is_low": np.nan,
            "empirical_p_lower_tail": np.nan,
            "empirical_p_upper_tail": np.nan,
        }
    lower_p = float(np.mean(values <= observed))
    upper_p = float(np.mean(values >= observed))
    return {
        "n_candidate_regions": int(len(values)),
        "percentile_lower_is_low": float(100.0 * lower_p),
        "empirical_p_lower_tail": lower_p,
        "empirical_p_upper_tail": upper_p,
    }


def star_for_empirical(row, direction):
    if direction == "lower":
        p = row.get("empirical_p_lower_tail", np.nan)
    else:
        p = row.get("empirical_p_upper_tail", np.nan)
    if not np.isfinite(p):
        return ""
    if p <= 0.01:
        return "***"
    if p <= 0.05:
        return "**"
    if p <= 0.10:
        return "*"
    return ""




def fit_linear(x, y):
    df = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(df) < 3 or df["x"].nunique() < 2:
        return {
            "slope": np.nan,
            "intercept": np.nan,
            "p_value": np.nan,
            "r2": np.nan,
            "n_points": int(len(df)),
        }
    res = linregress(df["x"].to_numpy(float), df["y"].to_numpy(float))
    return {
        "slope": float(res.slope),
        "intercept": float(res.intercept),
        "p_value": float(res.pvalue),
        "r2": float(res.rvalue ** 2),
        "n_points": int(len(df)),
    }


def compact_sci(x, digits=2):
    """Format scientific notation as 1.00e-4 rather than 1.00e-04."""
    if not np.isfinite(x):
        return "NA"
    txt = f"{float(x):.{digits}e}"
    return re.sub(r"e([+-])0*(\d+)$", r"e\1\2", txt)


def p_text(p):
    if not np.isfinite(p):
        return r"$P$ = NA"
    return rf"$P$ = {compact_sci(p, 2)}"


def stable_seed(base_seed, *parts):
    """Generate a reproducible 31-bit seed independent of PYTHONHASHSEED."""
    payload = "|".join([str(base_seed), *(str(part) for part in parts)]).encode()
    digest = hashlib.blake2b(payload, digest_size=8).digest()
    return int.from_bytes(digest, "little") % (2**31 - 1) or 1


# -----------------------------
# Input readers
# -----------------------------

def read_ind(ind_path):
    samples = []
    with open_text(ind_path) as f:
        for line in f:
            if line.strip():
                samples.append(line.split()[0])
    return samples


def read_snp(snp_path):
    df = pd.read_csv(
        snp_path,
        sep=r"\s+",
        header=None,
        names=["snp_id", "chrom", "genetic_pos", "pos", "ref", "alt"],
        dtype={"snp_id": str, "chrom": str, "genetic_pos": str, "pos": int, "ref": str, "alt": str},
    )
    df = df[df["chrom"].astype(str).isin(["11", "chr11"])].copy()
    df = df.reset_index(drop=False).rename(columns={"index": "original_index"})
    return df


def load_genotypes(geno_path, n_samples, keep_original_indices):
    keep = set(int(i) for i in keep_original_indices)
    rows = []
    with open_text(geno_path) as f:
        for i, line in enumerate(f):
            if i not in keep:
                continue
            g = line.strip()
            if len(g) != n_samples:
                raise RuntimeError(f"EIGENSTRAT line {i+1} length {len(g)} != n_samples {n_samples}")
            arr = np.empty(n_samples, dtype=np.float32)
            # EIGENSTRAT code is number of REF alleles. Convert to ALT dosage.
            for j, ch in enumerate(g):
                if ch == "2":
                    arr[j] = 0.0
                elif ch == "1":
                    arr[j] = 0.5
                elif ch == "0":
                    arr[j] = 1.0
                elif ch == "9":
                    arr[j] = np.nan
                else:
                    raise RuntimeError(f"Unexpected EIGENSTRAT code {ch!r} at line {i+1}")
            rows.append(arr)
    if len(rows) == 0:
        return np.empty((0, n_samples), dtype=np.float32)
    return np.vstack(rows)


def detect_ignore_outlier(anno_path):
    anno = pd.read_csv(anno_path, sep="\t", dtype=str, low_memory=False)
    gid_col = find_col(anno.columns, startswith="genetic id")
    pattern = re.compile(r"(Ignore_|outlier)", flags=re.IGNORECASE)
    text_cols = [c for c in anno.columns if anno[c].dtype == object]
    flags = []
    for _, row in anno[text_cols].iterrows():
        hit = False
        for v in row.values:
            if isinstance(v, str) and pattern.search(v):
                hit = True
                break
        flags.append(hit)
    return pd.DataFrame({"Genetic ID": anno[gid_col].astype(str), "flag_ignore_or_outlier": flags})


def assign_group(lat, lon):
    if pd.isna(lat) or pd.isna(lon):
        return "Unassigned"
    # East has priority at the 80E boundary.
    if lon >= 80 and lat >= -10:
        return "East Eurasia (DR)"

    mask_europe = (lon >= -25) and (lon <= 45) and (lat >= 35)
    mask_west_asia = (lon >= 25) and (lon <= 80) and (lat >= 12)
    mask_greenland = (lon >= -75) and (lon <= -10) and (lat >= 58)
    if mask_europe or mask_west_asia or mask_greenland:
        return "West Eurasia (DR-L)"
    return "Unassigned"


def prepare_metadata(metadata_path, anno_path, samples):
    meta = pd.read_csv(metadata_path, sep="\t", dtype=str)
    required = ["Genetic ID", "Individual ID", "Latitude", "Longitude", "Date_BP"]
    missing = [c for c in required if c not in meta.columns]
    if missing:
        raise ValueError(f"Metadata missing required columns: {missing}")

    meta["Latitude"] = pd.to_numeric(meta["Latitude"], errors="coerce")
    meta["Longitude"] = pd.to_numeric(meta["Longitude"], errors="coerce")
    meta["Date_BP"] = pd.to_numeric(meta["Date_BP"], errors="coerce")
    meta["Age_BP"] = meta["Date_BP"]

    flags = detect_ignore_outlier(anno_path)
    meta = meta.merge(flags, on="Genetic ID", how="left")
    meta["flag_ignore_or_outlier"] = meta["flag_ignore_or_outlier"].fillna(False).astype(bool)

    sample_to_idx = {s: i for i, s in enumerate(samples)}
    meta = meta[meta["Genetic ID"].isin(sample_to_idx)].copy()
    meta["sample_index"] = meta["Genetic ID"].map(sample_to_idx).astype(int)

    meta["group"] = meta.apply(lambda x: assign_group(x["Latitude"], x["Longitude"]), axis=1)
    meta["qc_missing_date"] = meta["Age_BP"].isna()
    meta["qc_missing_location"] = meta["Latitude"].isna() | meta["Longitude"].isna()
    meta["qc_unassigned_group"] = meta["group"].eq("Unassigned")
    meta["qc_flagged"] = meta["flag_ignore_or_outlier"]
    meta["pass_qc"] = ~(
        meta["qc_missing_date"]
        | meta["qc_missing_location"]
        | meta["qc_unassigned_group"]
        | meta["qc_flagged"]
    )
    return meta


# -----------------------------
# Window and statistic functions
# -----------------------------

def make_windows(chrom_start, chrom_end_inclusive, window_bp, step_bp):
    max_half_open = chrom_end_inclusive + 1
    starts = np.arange(chrom_start, max_half_open - window_bp + 1, step_bp, dtype=np.int64)
    ends = starts + window_bp
    centers = (starts + ends) / 2.0
    return starts, ends, centers


def site_stats_for_samples(G, sample_idx):
    n_sites = G.shape[0]
    if len(sample_idx) == 0 or n_sites == 0:
        return {
            "n_called": np.zeros(n_sites, dtype=np.float32),
            "site_pi": np.full(n_sites, np.nan, dtype=np.float32),
            "seg": np.zeros(n_sites, dtype=np.float32),
            "valid": np.zeros(n_sites, dtype=np.float32),
        }

    X = G[:, sample_idx]
    called = np.isfinite(X)
    n_called = called.sum(axis=1).astype(np.float32)
    alt_sum = np.nansum(X, axis=1).astype(np.float32)

    valid = n_called >= 2
    p = np.full(n_sites, np.nan, dtype=np.float32)
    p[valid] = alt_sum[valid] / n_called[valid]

    site_pi = np.full(n_sites, np.nan, dtype=np.float32)
    site_pi[valid] = 2.0 * p[valid] * (1.0 - p[valid]) * (n_called[valid] / (n_called[valid] - 1.0))
    seg = (valid & (p > 0.0) & (p < 1.0)).astype(np.float32)

    return {
        "n_called": n_called,
        "site_pi": site_pi,
        "seg": seg,
        "valid": valid.astype(np.float32),
    }


def cumsum0(x):
    x = np.asarray(x, dtype=np.float64)
    return np.concatenate([[0.0], np.cumsum(x)])


def tajima_d_from_arrays(pi_sum, S, n_eff):
    pi_sum = np.asarray(pi_sum, dtype=float)
    S = np.asarray(S, dtype=float)
    n_eff = np.asarray(n_eff, dtype=float)
    out = np.full(pi_sum.shape, np.nan, dtype=float)

    ok = np.isfinite(pi_sum) & np.isfinite(S) & np.isfinite(n_eff) & (S >= 2) & (n_eff >= 4)
    if ok.sum() == 0:
        return out

    n_round = np.full(n_eff.shape, -1, dtype=int)
    n_round[ok] = np.round(n_eff[ok]).astype(int)
    unique_n = np.unique(n_round[ok])
    for n in unique_n:
        idx = ok & (n_round == n)
        if idx.sum() == 0 or n < 4:
            continue
        a1 = sum(1.0 / i for i in range(1, n))
        a2 = sum(1.0 / (i * i) for i in range(1, n))
        b1 = (n + 1.0) / (3.0 * (n - 1.0))
        b2 = 2.0 * (n * n + n + 3.0) / (9.0 * n * (n - 1.0))
        c1 = b1 - 1.0 / a1
        c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1)
        e1 = c1 / a1
        e2 = c2 / (a1 * a1 + a2)
        denom = np.sqrt(e1 * S[idx] + e2 * S[idx] * (S[idx] - 1.0))
        good = denom > 0
        vals = np.full(idx.sum(), np.nan, dtype=float)
        vals[good] = (pi_sum[idx][good] - S[idx][good] / a1) / denom[good]
        out[idx] = vals
    return out


def window_stats_from_site_stats(pos, site_stats, starts, ends):
    pos = np.asarray(pos, dtype=np.int64)
    start_idx = np.searchsorted(pos, starts, side="left")
    end_idx = np.searchsorted(pos, ends, side="left")

    valid = site_stats["valid"]
    site_pi = np.where(np.isfinite(site_stats["site_pi"]), site_stats["site_pi"], 0.0)
    seg = site_stats["seg"]
    n_called = np.where(site_stats["valid"] > 0, site_stats["n_called"], 0.0)

    cum_valid = cumsum0(valid)
    cum_pi = cumsum0(site_pi)
    cum_seg = cumsum0(seg)
    cum_n = cumsum0(n_called)

    n_callable_sites = cum_valid[end_idx] - cum_valid[start_idx]
    pi_sum = cum_pi[end_idx] - cum_pi[start_idx]
    S = cum_seg[end_idx] - cum_seg[start_idx]
    n_sum = cum_n[end_idx] - cum_n[start_idx]

    window_pi = np.full(len(starts), np.nan, dtype=float)
    ok = n_callable_sites > 0
    window_pi[ok] = pi_sum[ok] / n_callable_sites[ok]

    mean_n = np.full(len(starts), np.nan, dtype=float)
    mean_n[ok] = n_sum[ok] / n_callable_sites[ok]
    n_eff = np.round(mean_n)
    window_D = tajima_d_from_arrays(pi_sum=pi_sum, S=S, n_eff=n_eff)

    return pd.DataFrame({
        "window_start": starts,
        "window_end": ends,
        "center": (starts + ends) / 2.0,
        "n_callable_sites": n_callable_sites.astype(int),
        "pi_sum": pi_sum,
        "mean_n_called": mean_n,
        "n_segregating_sites": S.astype(int),
        "window_pi": window_pi,
        "window_tajimaD_like": window_D,
    })


def target_window_mask(window_df, region: Region):
    c = window_df["center"].to_numpy(float)
    mask = (c >= region.start) & (c < region.end)
    if region.exclude_svr_for_observed:
        mask &= ~((c >= SVR_START) & (c < SVR_END))
    return mask


def target_site_mask(pos, region: Region, window_bp):
    # Include sites that can contribute to windows centered in the target region.
    # A broad span is used for n_samples_with_any_call diagnostics and bootstrap speed.
    pad = window_bp
    mask = (pos >= region.start - pad) & (pos < region.end + pad)
    if region.exclude_svr_for_observed:
        mask &= ~((pos >= SVR_START) & (pos < SVR_END))
    return mask


def summarize_region_from_window_df(window_df, region: Region):
    mask = target_window_mask(window_df, region)
    d = window_df.loc[mask].copy()
    if d.empty:
        return {
            "n_windows_total": 0,
            "n_windows_callable_pi": 0,
            "n_windows_callable_tajima": 0,
            "weighted_windowed_pi": np.nan,
            "mean_windowed_tajimaD_like": np.nan,
            "sum_callable_sites_in_windows": 0,
            "mean_callable_sites_per_window": np.nan,
            "mean_called_samples_per_callable_site": np.nan,
        }

    pi_vals = d["window_pi"].to_numpy(float)
    weights = d["n_callable_sites"].to_numpy(float)
    d_vals = d["window_tajimaD_like"].to_numpy(float)

    return {
        "n_windows_total": int(len(d)),
        "n_windows_callable_pi": int(np.isfinite(pi_vals).sum()),
        "n_windows_callable_tajima": int(np.isfinite(d_vals).sum()),
        "weighted_windowed_pi": safe_weighted_mean(pi_vals, weights),
        "mean_windowed_tajimaD_like": safe_mean(d_vals),
        "sum_callable_sites_in_windows": int(np.nansum(weights)),
        "mean_callable_sites_per_window": safe_mean(weights[weights > 0]),
        "mean_called_samples_per_callable_site": safe_mean(d["mean_n_called"].to_numpy(float)),
    }


def samples_with_any_call(G, pos, sample_idx, region: Region, window_bp):
    if len(sample_idx) == 0:
        return 0
    sm = target_site_mask(pos, region, window_bp)
    if sm.sum() == 0:
        return 0
    X = G[sm, :][:, sample_idx]
    return int((np.isfinite(X).sum(axis=0) > 0).sum())


def localize_windows_and_sites_for_region(G, pos, starts, ends, region: Region):
    """Return local G/pos/windows sufficient for observed target-region summaries.

    The original implementation recomputed all chr11 windows for every bootstrap
    replicate. For target-region bootstrap CIs this is unnecessary and extremely
    slow. We keep only windows whose centers fall in the target span and SNPs
    that can contribute to those windows. The region-specific SVR exclusion is
    still applied later by target_window_mask().
    """
    centers = (starts + ends) / 2.0
    wmask = (centers >= region.start) & (centers < region.end)
    if not np.any(wmask):
        return G[:0, :], pos[:0], starts[:0], ends[:0]
    ls = starts[wmask]
    le = ends[wmask]
    min_start = int(np.min(ls))
    max_end = int(np.max(le))
    smask = (pos >= min_start) & (pos < max_end)
    return G[smask, :], pos[smask], ls, le


def compute_target_summary(G, pos, sample_idx, starts, ends, region: Region, window_bp):
    G_loc, pos_loc, starts_loc, ends_loc = localize_windows_and_sites_for_region(G, pos, starts, ends, region)
    stats = site_stats_for_samples(G_loc, sample_idx)
    wdf = window_stats_from_site_stats(pos_loc, stats, starts_loc, ends_loc)
    out = summarize_region_from_window_df(wdf, region)
    out["n_samples"] = int(len(sample_idx))
    out["n_samples_with_any_call"] = samples_with_any_call(G, pos, sample_idx, region, window_bp)
    return out


def bootstrap_target_ci(G, pos, sample_idx, starts, ends, region, window_bp, n_boot, threads, seed):
    if n_boot <= 0 or len(sample_idx) == 0:
        return {
            "weighted_windowed_pi_boot_low": np.nan,
            "weighted_windowed_pi_boot_high": np.nan,
            "mean_windowed_tajimaD_like_boot_low": np.nan,
            "mean_windowed_tajimaD_like_boot_high": np.nan,
        }

    rng = np.random.default_rng(seed)
    sample_idx = np.array(sample_idx, dtype=int)
    boot_indices = [rng.choice(sample_idx, size=len(sample_idx), replace=True) for _ in range(n_boot)]

    # Localize once for the target region. This is the main speed-up: bootstrap
    # no longer recomputes all chr11 windows.
    G_loc, pos_loc, starts_loc, ends_loc = localize_windows_and_sites_for_region(G, pos, starts, ends, region)

    def one(bidx):
        stats = site_stats_for_samples(G_loc, bidx)
        wdf = window_stats_from_site_stats(pos_loc, stats, starts_loc, ends_loc)
        out = summarize_region_from_window_df(wdf, region)
        return out

    rows = []
    if threads <= 1:
        for b in boot_indices:
            rows.append(one(b))
    else:
        with ThreadPoolExecutor(max_workers=threads) as ex:
            futs = [ex.submit(one, b) for b in boot_indices]
            for fut in as_completed(futs):
                rows.append(fut.result())

    df = pd.DataFrame(rows)
    out = {}
    for metric in ["weighted_windowed_pi", "mean_windowed_tajimaD_like"]:
        vals = df[metric].to_numpy(float) if metric in df else np.array([])
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            out[f"{metric}_boot_low"] = np.nan
            out[f"{metric}_boot_high"] = np.nan
        else:
            out[f"{metric}_boot_low"] = float(np.percentile(vals, 2.5))
            out[f"{metric}_boot_high"] = float(np.percentile(vals, 97.5))
    return out


# -----------------------------
# Candidate-region background
# -----------------------------

def candidate_region_summaries(window_df, region_length, chrom_start, chrom_end_inclusive, candidate_step_bp):
    max_half_open = chrom_end_inclusive + 1
    cand_starts = np.arange(chrom_start, max_half_open - region_length + 1, candidate_step_bp, dtype=np.int64)
    cand_ends = cand_starts + region_length
    centers = window_df["center"].to_numpy(float)

    left = np.searchsorted(centers, cand_starts, side="left")
    right = np.searchsorted(centers, cand_ends, side="left")

    def interval_mean(values):
        values = np.asarray(values, dtype=float)
        valid = np.isfinite(values).astype(float)
        vals = np.where(np.isfinite(values), values, 0.0)
        cv = cumsum0(vals)
        cc = cumsum0(valid)
        s = cv[right] - cv[left]
        c = cc[right] - cc[left]
        out = np.full(len(cand_starts), np.nan, dtype=float)
        ok = c > 0
        out[ok] = s[ok] / c[ok]
        return out, c

    pi_vals = window_df["window_pi"].to_numpy(float)
    d_vals = window_df["window_tajimaD_like"].to_numpy(float)
    ncall = window_df["n_callable_sites"].to_numpy(float)

    mean_d, n_d_windows = interval_mean(d_vals)

    valid_w = np.isfinite(pi_vals) & np.isfinite(ncall) & (ncall > 0)
    wv = np.where(valid_w, pi_vals * ncall, 0.0)
    ww = np.where(valid_w, ncall, 0.0)
    cwv = cumsum0(wv)
    cww = cumsum0(ww)
    sum_wv = cwv[right] - cwv[left]
    sum_w = cww[right] - cww[left]
    weighted_pi = np.full(len(cand_starts), np.nan, dtype=float)
    ok = sum_w > 0
    weighted_pi[ok] = sum_wv[ok] / sum_w[ok]

    return pd.DataFrame({
        "candidate_start": cand_starts,
        "candidate_end": cand_ends,
        "n_windows_for_tajima": n_d_windows.astype(int),
        "weighted_windowed_pi": weighted_pi,
        "mean_windowed_tajimaD_like": mean_d,
    })


def add_empirical_percentiles(summary, candidate_df):
    out = dict(summary)
    for metric in ["weighted_windowed_pi", "mean_windowed_tajimaD_like"]:
        emp = empirical_percentiles(candidate_df[metric].to_numpy(float), summary.get(metric, np.nan))
        for k, v in emp.items():
            out[f"{metric}_{k}"] = v
    return out



# -----------------------------
# Time binning and rolling
# -----------------------------



def regular_bin_edges(max_age_bp, bin_size_years):
    end_age = int(math.ceil(float(max_age_bp) / bin_size_years) * bin_size_years)
    return np.arange(0, end_age + bin_size_years, bin_size_years, dtype=int)


def format_bin_label(start_bp, end_bp, is_tail=False):
    if is_tail:
        return f"> {start_bp / 1000:g} Kya"
    if start_bp == 0:
        return f"< {end_bp / 1000:g} Kya"
    return f"{start_bp / 1000:g}-{end_bp / 1000:g} Kya"


def add_target_any_call_flags(meta, G, pos, window_bp):
    """Add has_any_call_target using the actual target region for each plotted group.

    This is used only for adaptive deep-time bin merging and n diagnostics. The
    metric itself is still computed from all samples in the bin; samples with no
    target calls simply contribute no allele observations.
    """
    out = meta.copy()
    out["has_any_call_target"] = False
    for group_name in GROUP_STYLE:
        region = REGIONS[GROUP_STYLE[group_name]["region_name"]]
        mask = out["group"].eq(group_name)
        if not mask.any():
            continue
        sample_idx = out.loc[mask, "sample_index"].to_numpy(dtype=int)
        sm = target_site_mask(pos, region, window_bp)
        if sm.sum() == 0 or len(sample_idx) == 0:
            continue
        X = G[sm, :][:, sample_idx]
        flags = (np.isfinite(X).sum(axis=0) > 0)
        out.loc[mask, "has_any_call_target"] = flags
    return out


def assign_adaptive_bins_for_group(gmeta, max_age_bp, bin_size_years, low_n_threshold=10):
    """Assign bins for one group with deep-time tail merging.

    Starting from recent to ancient regular bins, the first bin with fewer than
    `low_n_threshold` samples that have any target-region call defines the tail
    boundary. That bin and all older samples are merged into a single "> X Kya"
    bin. This implements the requested "cutting-end" rule and prevents long
    sparse Palaeolithic bins from dominating the visual trajectory.
    """
    tmp = gmeta.copy()
    if tmp.empty:
        return tmp
    edges = regular_bin_edges(max_age_bp, bin_size_years)
    # Count samples that can actually contribute to the target metric.
    tail_start = None
    for start_bp, end_bp in zip(edges[:-1], edges[1:]):
        m = (tmp["Age_BP"] >= start_bp) & (tmp["Age_BP"] < end_bp)
        n_any = int(tmp.loc[m, "has_any_call_target"].sum()) if m.any() else 0
        if n_any < low_n_threshold:
            tail_start = int(start_bp)
            break
    # If all regular bins are sufficiently supported, do not merge unless there
    # are samples older than the last edge.
    if tail_start is None:
        tail_start = int(edges[-1])

    labels = []
    bin_start = []
    bin_end = []
    is_tail = []
    for age in tmp["Age_BP"].to_numpy(float):
        if age >= tail_start:
            labels.append(format_bin_label(tail_start, np.inf, is_tail=True))
            bin_start.append(float(tail_start))
            bin_end.append(np.inf)
            is_tail.append(True)
        else:
            k = int(math.floor(age / bin_size_years))
            st = int(k * bin_size_years)
            en = int((k + 1) * bin_size_years)
            labels.append(format_bin_label(st, en, is_tail=False))
            bin_start.append(float(st))
            bin_end.append(float(en))
            is_tail.append(False)
    tmp["time_label"] = labels
    tmp["bin_start_years"] = bin_start
    tmp["bin_end_years"] = bin_end
    tmp["is_tail_bin"] = is_tail
    # Continuous plotting coordinate; older-left axis is handled by inverting x.
    # The final point is placed by mean_age_bp in the summary, not by category.
    # Keep a simple plot code for bookkeeping only.
    order = (tmp[["time_label", "bin_start_years"]]
             .drop_duplicates()
             .sort_values("bin_start_years", ascending=False)
             .reset_index(drop=True))
    code_map = {r.time_label: i for i, r in order.iterrows()}
    tmp["plot_code"] = tmp["time_label"].map(code_map).astype(int)
    return tmp


def rolling_centers(max_age_bp, width_bp, step_bp):
    max_age = int(math.ceil(max_age_bp / step_bp) * step_bp)
    return np.arange(0, max_age + step_bp, step_bp, dtype=int)


# -----------------------------
# Core analysis functions
# -----------------------------

def analyze_one_time_unit(
    G,
    pos,
    starts,
    ends,
    samples_meta,
    group_name,
    region: Region,
    time_label,
    plot_code,
    unit_type,
    bin_size_years,
    n_boot,
    threads,
    seed,
    do_background,
    candidate_step_bp,
    chrom_start,
    chrom_end,
    window_bp,
):
    sample_idx = samples_meta["sample_index"].to_numpy(dtype=int)
    summary = compute_target_summary(G, pos, sample_idx, starts, ends, region, window_bp)
    summary.update({
        "unit_type": unit_type,
        "group": group_name,
        "region": region.name,
        "time_label": time_label,
        "plot_code": plot_code,
        "bin_size_years": bin_size_years,
        "mean_age_bp": float(samples_meta["Age_BP"].mean()) if len(samples_meta) else np.nan,
        "median_age_bp": float(samples_meta["Age_BP"].median()) if len(samples_meta) else np.nan,
        "min_age_bp": float(samples_meta["Age_BP"].min()) if len(samples_meta) else np.nan,
        "max_age_bp": float(samples_meta["Age_BP"].max()) if len(samples_meta) else np.nan,
        "target_start": region.start,
        "target_end": region.end,
        "background_length_bp": region.background_length,
        "target_excludes_svr_centered_windows": int(region.exclude_svr_for_observed),
    })


    # Bootstrap observed target only.
    ci = bootstrap_target_ci(
        G=G,
        pos=pos,
        sample_idx=sample_idx,
        starts=starts,
        ends=ends,
        region=region,
        window_bp=window_bp,
        n_boot=n_boot,
        threads=threads,
        seed=seed,
    )
    summary.update(ci)

    candidate_df = None
    if do_background and len(sample_idx) > 0:
        site_stats = site_stats_for_samples(G, sample_idx)
        wdf = window_stats_from_site_stats(pos, site_stats, starts, ends)
        candidate_df = candidate_region_summaries(
            window_df=wdf,
            region_length=region.background_length,
            chrom_start=chrom_start,
            chrom_end_inclusive=chrom_end,
            candidate_step_bp=candidate_step_bp,
        )
        summary = add_empirical_percentiles(summary, candidate_df)
    return summary, candidate_df


# -----------------------------
# Plotting
# -----------------------------

def get_metric_ci_cols(metric):
    return f"{metric}_boot_low", f"{metric}_boot_high"


def marker_for_summary(row, marker_mode):
    if marker_mode != "empirical":
        return ""
    direction = GROUP_STYLE[row["group"]]["direction"]
    metric = row.get("_plot_metric")
    if not metric:
        return ""
    empirical = {
        "empirical_p_lower_tail": row.get(f"{metric}_empirical_p_lower_tail", np.nan),
        "empirical_p_upper_tail": row.get(f"{metric}_empirical_p_upper_tail", np.nan),
    }
    return star_for_empirical(empirical, direction)



def smooth_fit_values(x, y, method="lowess", frac=0.45):
    """Return fitted values at the observed x positions for a smooth trend."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x0 = x[ok]
    y0 = y[ok]
    yhat = np.full_like(y, np.nan, dtype=float)
    if len(x0) < 3 or len(np.unique(x0)) < 2:
        return yhat
    if method == "linear" or len(np.unique(x0)) < 3:
        fit = fit_linear(-x0, y0)
        if np.isfinite(fit.get("slope", np.nan)):
            yhat[ok] = fit["intercept"] + fit["slope"] * (-x0)
        return yhat
    if sm_lowess is not None and len(x0) >= 4 and len(np.unique(x0)) >= 4:
        order = np.argsort(x0)
        xs = x0[order]
        ys = y0[order]
        this_frac = 1.0 if len(xs) <= 6 else min(max(frac, 0.2), 1.0)
        sm = sm_lowess(ys, xs, frac=this_frac, it=1, return_sorted=True)
        yhat[ok] = np.interp(x0, sm[:, 0], sm[:, 1])
        return yhat
    deg = 2 if len(x0) >= 4 and len(np.unique(x0)) >= 3 else 1
    coef = np.polyfit(x0, y0, deg=deg)
    yhat[ok] = np.polyval(coef, x0)
    return yhat


def pseudo_r2_from_yhat(y, yhat):
    y = np.asarray(y, dtype=float)
    yhat = np.asarray(yhat, dtype=float)
    ok = np.isfinite(y) & np.isfinite(yhat)
    if ok.sum() < 3:
        return np.nan
    yy = y[ok]
    yh = yhat[ok]
    sst = float(np.sum((yy - np.mean(yy)) ** 2))
    if sst <= 0:
        return np.nan
    sse = float(np.sum((yy - yh) ** 2))
    return max(0.0, min(1.0, 1.0 - sse / sst))


def lowess_permutation_test(x, y, n_perm=999, seed=1, frac=0.45):
    """Permutation P for any age-associated smooth pattern.

    Statistic = LOWESS pseudo-R2 against an intercept-only null. The test is
    non-directional and does not assume a linear trajectory. Bootstrap replicates
    are intentionally not used as independent observations.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    if len(x) < 4 or len(np.unique(x)) < 3:
        return {"slope": np.nan, "intercept": np.nan, "r_value": np.nan,
                "p_value": np.nan, "r2": np.nan, "n_points": int(len(x)),
                "fit_test": "lowess_perm", "test_stat": np.nan,
                "n_permutations": int(n_perm)}
    yhat = smooth_fit_values(x, y, method="lowess", frac=frac)
    stat_obs = pseudo_r2_from_yhat(y, yhat)
    if not np.isfinite(stat_obs):
        return {"slope": np.nan, "intercept": np.nan, "r_value": np.nan,
                "p_value": np.nan, "r2": np.nan, "n_points": int(len(x)),
                "fit_test": "lowess_perm", "test_stat": np.nan,
                "n_permutations": int(n_perm)}
    rng = np.random.default_rng(seed)
    ge = 0
    n_done = 0
    for _ in range(int(n_perm)):
        yp = rng.permutation(y)
        yhp = smooth_fit_values(x, yp, method="lowess", frac=frac)
        st = pseudo_r2_from_yhat(yp, yhp)
        if np.isfinite(st):
            n_done += 1
            if st >= stat_obs - 1e-15:
                ge += 1
    pval = (ge + 1.0) / (n_done + 1.0) if n_done > 0 else np.nan
    lf = fit_linear(-x, y)  # descriptive linear direction only
    return {"slope": lf.get("slope", np.nan), "intercept": lf.get("intercept", np.nan),
            "r_value": lf.get("r_value", np.nan), "p_value": float(pval),
            "r2": float(stat_obs), "n_points": int(len(x)),
            "fit_test": "lowess_perm", "test_stat": float(stat_obs),
            "n_permutations": int(n_done)}


def fit_points(df, metric, x_col, fit_min_n=0, trend_test="lowess_perm",
               n_perm=999, seed=1):
    d = df.copy()
    d = d[np.isfinite(pd.to_numeric(d[metric], errors="coerce"))]
    d = d[np.isfinite(pd.to_numeric(d[x_col], errors="coerce"))]
    if fit_min_n > 0 and "n_samples_with_any_call" in d.columns:
        d = d[pd.to_numeric(d["n_samples_with_any_call"], errors="coerce") >= fit_min_n]
    x = pd.to_numeric(d[x_col], errors="coerce").to_numpy(float)
    y = pd.to_numeric(d[metric], errors="coerce").to_numpy(float)
    if trend_test == "none":
        return {"slope": np.nan, "intercept": np.nan, "r_value": np.nan,
                "p_value": np.nan, "r2": np.nan, "n_points": int(np.isfinite(y).sum()),
                "fit_test": "none", "test_stat": np.nan, "n_permutations": 0}
    if trend_test == "ols":
        out = fit_linear(-x, y)
        out["fit_test"] = "ols"
        out["test_stat"] = out.get("r2", np.nan)
        out["n_permutations"] = 0
        return out
    return lowess_permutation_test(x, y, n_perm=n_perm, seed=seed)


def add_fit_text(ax, fit_rows, font_size=12, x=0.02, y=0.98,
                 ha="left", va="top", line_step=0.08):
    for row in fit_rows:
        group = row["group"]
        label = GROUP_DISPLAY_LABELS.get(group, group)
        color = GROUP_STYLE[group]["color"]
        p = row.get("p_value", np.nan)
        if np.isfinite(row.get("r2", np.nan)):
            txt = rf"{label}: {p_text(p)}, $R^2$ = {row['r2']:.3f}"
        else:
            txt = f"{label}: fit NA"
        ax.text(x, y, txt, transform=ax.transAxes, va=va, ha=ha, fontsize=font_size, color=color)
        y -= line_step


def lowess_curve(x, y, frac=0.45):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    if len(x) < 3 or len(np.unique(x)) < 3:
        return None, None
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    if sm_lowess is not None and len(x) >= 4:
        this_frac = 1.0 if len(x) <= 6 else min(max(frac, 0.2), 1.0)
        sm = sm_lowess(y, x, frac=this_frac, it=1, return_sorted=True)
        return sm[:, 0], sm[:, 1]
    deg = 2 if len(x) >= 4 else 1
    coef = np.polyfit(x, y, deg=deg)
    xs = np.linspace(np.nanmin(x), np.nanmax(x), 200)
    ys = np.polyval(coef, xs)
    return xs, ys


def plot_smooth(ax, x, y, color, method="lowess", alpha=0.85):
    if method == "none":
        return
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    if len(x) < 3 or len(np.unique(x)) < 2:
        return
    if method == "linear":
        fit = fit_linear(-x, y)
        if np.isfinite(fit["slope"]):
            xs = np.linspace(np.nanmin(x), np.nanmax(x), 200)
            ys = fit["intercept"] + fit["slope"] * (-xs)
            ax.plot(xs, ys, color=color, linestyle="--", linewidth=1.4, alpha=alpha)
        return
    xs, ys = lowess_curve(x, y)
    if xs is not None:
        ax.plot(xs, ys, color=color, linestyle="--", linewidth=1.6, alpha=alpha)


def recompute_plot_age(df, tail_bin_position="align-oldest"):
    """Set plot_age_kya for binned and rolling rows.

    tail_bin_position:
      mean         : tail bin at its sample mean age.
      boundary     : tail bin at its own lower boundary (>X plotted at X).
      align-oldest : all group-specific tail bins for a bin size share the
                     oldest tail boundary in that panel, e.g. East >9 and West
                     >13 are both plotted at 13 kya.
    """
    d = df.copy()
    d["plot_age_kya"] = np.nan
    is_roll = d["unit_type"].eq("rolling")
    if "age_center_years" in d.columns:
        d.loc[is_roll, "plot_age_kya"] = pd.to_numeric(d.loc[is_roll, "age_center_years"], errors="coerce") / 1000.0
    is_bin = d["unit_type"].eq("binned")
    if is_bin.any():
        d.loc[is_bin, "plot_age_kya"] = pd.to_numeric(d.loc[is_bin, "mean_age_bp"], errors="coerce") / 1000.0
        if "is_tail_bin" in d.columns and "bin_start_years" in d.columns:
            tail = is_bin & (pd.to_numeric(d["is_tail_bin"], errors="coerce").fillna(0).astype(int) == 1)
            if tail.any():
                if tail_bin_position == "mean":
                    pass
                elif tail_bin_position == "boundary":
                    d.loc[tail, "plot_age_kya"] = pd.to_numeric(d.loc[tail, "bin_start_years"], errors="coerce") / 1000.0
                elif tail_bin_position == "align-oldest":
                    for bs, subidx in d.loc[tail].groupby("bin_size_years").groups.items():
                        idx = list(subidx)
                        starts = pd.to_numeric(d.loc[idx, "bin_start_years"], errors="coerce")
                        if starts.notna().any():
                            d.loc[idx, "plot_age_kya"] = float(starts.max()) / 1000.0
                else:
                    raise ValueError(f"Unknown tail_bin_position: {tail_bin_position}")
    d["x_age_kya"] = d["plot_age_kya"]
    return d


def prepare_plot_df(df, max_age_bp=None, rolling_max_age_bp=None, tail_bin_position="align-oldest"):
    d = recompute_plot_age(df, tail_bin_position=tail_bin_position)
    if max_age_bp is not None:
        max_kya = float(max_age_bp) / 1000.0
        d = d[pd.to_numeric(d["plot_age_kya"], errors="coerce") <= max_kya]
    if rolling_max_age_bp is not None:
        is_roll = d["unit_type"].eq("rolling")
        keep_roll = pd.to_numeric(d.get("age_center_years", np.nan), errors="coerce") <= float(rolling_max_age_bp)
        d = d[(~is_roll) | (is_roll & keep_roll)]
    return d.copy()


def compute_metric_ylim(df, metric, include_ci=True, force_min=None):
    vals = []
    if metric in df:
        vals.extend(pd.to_numeric(df[metric], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna().tolist())
    if include_ci:
        low_col, high_col = get_metric_ci_cols(metric)
        for col in (low_col, high_col):
            if col in df:
                vals.extend(pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna().tolist())
    if not vals:
        return None
    ymin = float(np.nanmin(vals))
    ymax = float(np.nanmax(vals))
    if force_min is not None:
        ymin = min(ymin, float(force_min))
    if not np.isfinite(ymin) or not np.isfinite(ymax):
        return None
    if ymin == ymax:
        pad = max(abs(ymin) * 0.1, 0.05)
    else:
        pad = (ymax - ymin) * 0.10
    return ymin - pad, ymax + pad


def plot_binned_panel(ax, df, metric, marker_mode="none", annotate_n=False,
                      annotate_stars=True, low_n_threshold=10, trend_method="lowess",
                      fit_min_n=0, trend_test="lowess_perm", trend_permutations=999,
                      trend_seed=1, y_limits=None):
    fit_rows = []
    metric_label = METRIC_LABELS[metric]
    for group, gd in df.groupby("group", sort=False):
        color = GROUP_STYLE[group]["color"]
        gd = gd.copy()
        gd = gd[np.isfinite(pd.to_numeric(gd[metric], errors="coerce"))]
        if gd.empty:
            continue
        if "x_age_kya" not in gd.columns:
            gd["x_age_kya"] = pd.to_numeric(gd["mean_age_bp"], errors="coerce") / 1000.0
        gd = gd[np.isfinite(gd["x_age_kya"])].sort_values("x_age_kya", ascending=False)
        if gd.empty:
            continue
        x = gd["x_age_kya"].to_numpy(float)
        y = gd[metric].to_numpy(float)
        low_col, high_col = get_metric_ci_cols(metric)
        ylow = gd[low_col].to_numpy(float) if low_col in gd else np.full(len(gd), np.nan)
        yhigh = gd[high_col].to_numpy(float) if high_col in gd else np.full(len(gd), np.nan)
        yerr_low = np.nan_to_num(np.maximum(0.0, y - ylow), nan=0.0, posinf=0.0, neginf=0.0)
        yerr_high = np.nan_to_num(np.maximum(0.0, yhigh - y), nan=0.0, posinf=0.0, neginf=0.0)
        ax.errorbar(
            x, y, yerr=[yerr_low, yerr_high], fmt="o", linestyle="none",
            color=color, ecolor=color, markersize=7.0, markeredgewidth=0,
            elinewidth=1.2, capsize=0, alpha=0.90,
            label=GROUP_DISPLAY_LABELS.get(group, group), zorder=3,
        )
        plot_smooth(ax, x, y, color, method=trend_method, alpha=0.90)
        if annotate_n:
            for _, r in gd.iterrows():
                n = int(r.get("n_samples_with_any_call", 0))
                suffix = "†" if n < low_n_threshold else ""
                ax.text(r["x_age_kya"], r[metric], f"n={n}{suffix}", fontsize=5.8,
                        color="#666666", ha="center", va="bottom")
        if annotate_stars and marker_mode == "empirical":
            y_range = np.nanmax(y) - np.nanmin(y) if len(y) else 0
            offset = 0.04 * y_range if np.isfinite(y_range) and y_range > 0 else 0.02
            for _, r in gd.iterrows():
                rr = r.to_dict(); rr["_plot_metric"] = metric
                star = marker_for_summary(rr, marker_mode)
                if star:
                    ax.text(r["x_age_kya"], r[metric] + offset, star, fontsize=8.5,
                            color="black", ha="center", va="bottom", fontweight="bold")
        seed_i = stable_seed(trend_seed, group, metric, gd.get("bin_size_years", pd.Series([np.nan])).iloc[0])
        fit = fit_points(gd, metric, "x_age_kya", fit_min_n=fit_min_n,
                         trend_test=trend_test, n_perm=trend_permutations, seed=seed_i)
        fit.update({"group": group, "metric": metric})
        fit_rows.append(fit)
    ax.set_ylabel(metric_label, fontsize=14)
    ax.set_xlabel("Age (kya BP)", fontsize=14)
    if y_limits is not None:
        ax.set_ylim(*y_limits)
    ax.invert_xaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=12)
    return fit_rows


def plot_rolling_page(pdf, roll_df, metric, title_suffix, marker_mode="none",
                      low_n_threshold=10, trend_method="lowess", fit_min_n=0,
                      rolling_plot_min_n=1, show_rolling_ci=True,
                      trend_test="lowess_perm", trend_permutations=999,
                      trend_seed=1):
    d = roll_df.copy()
    d = d[np.isfinite(pd.to_numeric(d[metric], errors="coerce"))]
    if d.empty:
        return []
    if "x_age_kya" not in d.columns:
        d["x_age_kya"] = pd.to_numeric(d["age_center_years"], errors="coerce") / 1000.0
    d = d[np.isfinite(d["x_age_kya"])]
    dplot = d[pd.to_numeric(d["n_samples_with_any_call"], errors="coerce") >= rolling_plot_min_n].copy()
    if dplot.empty:
        dplot = d.copy()

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    fit_rows = []
    for group, gd in dplot.groupby("group", sort=False):
        color = GROUP_STYLE[group]["color"]
        gd = gd.sort_values("x_age_kya")
        x = gd["x_age_kya"].to_numpy(float)
        y = gd[metric].to_numpy(float)
        low_col, high_col = get_metric_ci_cols(metric)
        if show_rolling_ci and low_col in gd and high_col in gd:
            ylow = gd[low_col].to_numpy(float)
            yhigh = gd[high_col].to_numpy(float)
            yerr_low = np.nan_to_num(np.maximum(0.0, y - ylow), nan=0.0, posinf=0.0, neginf=0.0)
            yerr_high = np.nan_to_num(np.maximum(0.0, yhigh - y), nan=0.0, posinf=0.0, neginf=0.0)
            ax.errorbar(
                x, y, yerr=[yerr_low, yerr_high], fmt="o", linestyle="none",
                color=color, ecolor=color, markersize=6.0, markeredgewidth=0,
                elinewidth=0.8, capsize=0, alpha=0.80,
                label=GROUP_DISPLAY_LABELS.get(group, group), zorder=2,
            )
        else:
            ax.scatter(x, y, color=color, s=36, alpha=0.85, linewidths=0,
                       label=GROUP_DISPLAY_LABELS.get(group, group), zorder=2)
        plot_smooth(ax, x, y, color, method=trend_method, alpha=0.95)
        seed_i = stable_seed(trend_seed, group, metric, "rolling")
        fit = fit_points(gd, metric, "x_age_kya", fit_min_n=fit_min_n,
                         trend_test=trend_test, n_perm=trend_permutations, seed=seed_i)
        fit.update({"unit_type": "rolling", "bin_size_years": np.nan, "group": group, "metric": metric})
        fit_rows.append(fit)
    add_fit_text(ax, fit_rows, font_size=12, x=0.98, y=0.56, ha="right", va="center", line_step=0.10)
    ax.set_title(title_suffix, fontsize=16)
    ax.set_xlabel("Age (kya BP)", fontsize=14)
    ax.set_ylabel(METRIC_LABELS[metric], fontsize=14)
    if metric == "mean_windowed_tajimaD_like":
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(bottom=-3.0, top=ymax if ymax > -2.5 else 2.0)
    ax.invert_xaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=12)
    ax.legend(frameon=False, fontsize=12, loc="upper right")
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)
    return fit_rows


def plot_binned_pages(pdf, binned_df, metric, marker_mode, title_suffix,
                      panels_per_page=6, annotate_n=False, annotate_stars=True,
                      low_n_threshold=10, trend_method="lowess", fit_min_n=0,
                      trend_test="lowess_perm", trend_permutations=999,
                      trend_seed=1):
    bin_sizes = sorted(binned_df["bin_size_years"].dropna().unique())
    all_fit_rows = []
    for page_start in range(0, len(bin_sizes), panels_per_page):
        subset_sizes = bin_sizes[page_start:page_start + panels_per_page]
        n_panels = len(subset_sizes)
        ncols = 3
        nrows = int(math.ceil(n_panels / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5.0 * ncols, 4.6 * nrows), squeeze=False)
        axes_flat = axes.ravel()
        page_df = binned_df[binned_df["bin_size_years"].isin(subset_sizes)].copy()
        y_limits = compute_metric_ylim(page_df, metric, include_ci=True)
        for i, bs in enumerate(subset_sizes):
            ax = axes_flat[i]
            d = binned_df[binned_df["bin_size_years"] == bs].copy()
            if d.empty:
                ax.set_axis_off(); continue
            rows = plot_binned_panel(ax, d, metric, marker_mode=marker_mode,
                                     annotate_n=annotate_n, annotate_stars=annotate_stars,
                                     low_n_threshold=low_n_threshold, trend_method=trend_method,
                                     fit_min_n=fit_min_n, trend_test=trend_test,
                                     trend_permutations=trend_permutations,
                                     trend_seed=trend_seed, y_limits=y_limits)
            if i % ncols != 0:
                ax.set_ylabel("")
                ax.tick_params(axis="y", labelleft=False)
            for r in rows:
                r.update({"unit_type": "binned", "bin_size_years": bs, "marker_mode": marker_mode})
            all_fit_rows.extend(rows)
            ax.set_title(f"{int(bs // 1000)}ky bins", fontsize=16)
            add_fit_text(ax, rows, font_size=12, x=0.98, y=0.56, ha="right", va="center", line_step=0.10)
        for ax in axes_flat[n_panels:]:
            ax.set_axis_off()
        handles, labels = axes_flat[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.995, 0.995), frameon=False, fontsize=12)
        fig.suptitle(title_suffix, fontsize=16)
        fig.tight_layout(rect=[0, 0.01, 1, 0.97])
        pdf.savefig(fig)
        plt.close(fig)
    return all_fit_rows


def make_metric_pdf(out_pdf, roll_df, binned_df, metric, marker_mode, title_suffix,
                    annotate_n=False, low_n_threshold=10, trend_method="lowess",
                    fit_min_n=0, rolling_plot_min_n=1, include_rolling=True,
                    show_rolling_ci=False, trend_test="lowess_perm",
                    trend_permutations=999, trend_seed=1):
    fit_rows = []
    with PdfPages(out_pdf) as pdf:
        if include_rolling:
            rolling_suffix = title_suffix
            fit_rows.extend(plot_rolling_page(pdf, roll_df, metric, rolling_suffix,
                                              marker_mode="none", low_n_threshold=low_n_threshold,
                                              trend_method=trend_method, fit_min_n=fit_min_n,
                                              rolling_plot_min_n=rolling_plot_min_n,
                                              show_rolling_ci=show_rolling_ci,
                                              trend_test=trend_test,
                                              trend_permutations=trend_permutations,
                                              trend_seed=trend_seed))
        fit_rows.extend(plot_binned_pages(pdf, binned_df, metric, marker_mode, title_suffix,
                                          annotate_n=annotate_n, annotate_stars=True,
                                          low_n_threshold=low_n_threshold, trend_method=trend_method,
                                          fit_min_n=fit_min_n,
                                          trend_test=trend_test,
                                          trend_permutations=trend_permutations,
                                          trend_seed=trend_seed))
    return fit_rows


def generate_plots_from_results(res, outdir, args):
    outdir = Path(outdir)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)

    plot_res = prepare_plot_df(
        res,
        max_age_bp=args.plot_max_age_years,
        rolling_max_age_bp=None,
        tail_bin_position="align-oldest",
    )
    plot_res.to_csv(
        outdir / "05.plot_ready.windowed_region_summary.tsv",
        sep="\t",
        index=False,
    )

    roll_df = plot_res[plot_res["unit_type"] == "rolling"].copy()
    bin_df = plot_res[plot_res["unit_type"] == "binned"].copy()

    plot_jobs = [
        (
            "weighted_windowed_pi",
            METRIC_SIMPLE_TITLES["weighted_windowed_pi"],
            "AADR.PGA.weighted_windowed_pi.pdf",
        ),
        (
            "mean_windowed_tajimaD_like",
            METRIC_SIMPLE_TITLES["mean_windowed_tajimaD_like"],
            "AADR.PGA.mean_windowed_tajimaD_like.pdf",
        ),
    ]

    all_fit_rows = []
    for metric, title_suffix, filename in plot_jobs:
        rows = make_metric_pdf(
            out_pdf=outdir / "plots" / filename,
            roll_df=roll_df,
            binned_df=bin_df,
            metric=metric,
            marker_mode="empirical",
            title_suffix=title_suffix,
            annotate_n=False,
            low_n_threshold=args.low_n_threshold,
            trend_method="lowess",
            fit_min_n=args.fit_min_n,
            rolling_plot_min_n=args.rolling_plot_min_n,
            include_rolling=True,
            show_rolling_ci=True,
            trend_test="lowess_perm",
            trend_permutations=args.trend_permutations,
            trend_seed=args.trend_seed,
        )
        for row in rows:
            row.update({
                "pdf": filename,
                "trend_method": "lowess",
                "trend_test": "lowess_perm",
                "trend_permutations": args.trend_permutations,
                "plot_max_age_years": args.plot_max_age_years,
                "fit_min_n": args.fit_min_n,
                "tail_bin_position": "align-oldest",
            })
        all_fit_rows.extend(rows)

    fit_df = pd.DataFrame(all_fit_rows)
    fit_df.to_csv(outdir / "06.plot_fit_summary.tsv", sep="\t", index=False)

    pd.DataFrame([{
        "plot_max_age_years": args.plot_max_age_years,
        "fit_min_n": args.fit_min_n,
        "trend_method": "lowess",
        "trend_test": "lowess_perm",
        "trend_permutations": args.trend_permutations,
        "trend_seed": args.trend_seed,
        "tail_bin_position": "align-oldest",
        "rolling_plot_min_n": args.rolling_plot_min_n,
        "low_n_threshold": args.low_n_threshold,
        "note": (
            "LOWESS curves are descriptive. P values are permutation-based tests "
            "of the smooth temporal association using observed time-unit summaries; "
            "bootstrap replicates are not treated as independent observations."
        ),
    }]).to_csv(outdir / "07.plot_config.tsv", sep="\t", index=False)
    return fit_df

# -----------------------------
# Main
# -----------------------------


def parse_args():
    parser = argparse.ArgumentParser(
        description="AADR v66 2M temporal diversity analysis at the PGA desert region."
    )
    parser.add_argument("--geno", required=True, help="Chromosome 11 EIGENSTRAT .geno")
    parser.add_argument("--snp", required=True, help="Chromosome 11 EIGENSTRAT .snp")
    parser.add_argument("--ind", required=True, help="Chromosome 11 EIGENSTRAT .ind")
    parser.add_argument("--metadata", required=True, help="Metadata from aadr_prepare_subset.py")
    parser.add_argument("--anno", required=True, help="AADR v66 annotation table")
    parser.add_argument("--outdir", required=True, help="Output directory")

    # Genomic windows and chromosome-wide matched-region background.
    parser.add_argument("--window-bp", type=int, default=10_000)
    parser.add_argument("--window-step-bp", type=int, default=1_000)
    parser.add_argument("--candidate-step-bp", type=int, default=1_000)

    # Temporal summaries. Values are years before present, not genomic base pairs.
    parser.add_argument("--rolling-width-years", type=int, default=1_000)
    parser.add_argument("--rolling-step-years", type=int, default=100)
    parser.add_argument("--rolling-max-age-years", type=float, default=None)
    parser.add_argument(
        "--bin-sizes",
        nargs="+",
        type=int,
        default=[1_000, 2_000, 3_000],
        help="Temporal bin sizes in years (default: 1000 2000 3000)",
    )

    parser.add_argument("--n-boot", type=int, default=1_000)
    parser.add_argument(
        "--threads",
        type=int,
        default=max(1, min(8, os.cpu_count() or 1)),
    )
    parser.add_argument("--seed", type=int, default=20260511)
    parser.add_argument(
        "--low-n-threshold",
        type=int,
        default=10,
        help="Merge this and all older bins when target-callable n falls below the threshold",
    )

    # Plot-level temporal trend test.
    parser.add_argument("--plot-max-age-years", type=float, default=None)
    parser.add_argument("--fit-min-n", type=int, default=0)
    parser.add_argument("--rolling-plot-min-n", type=int, default=10)
    parser.add_argument("--trend-permutations", type=int, default=9_999)
    parser.add_argument("--trend-seed", type=int, default=20260511)
    return parser.parse_args()

def main():
    configure_matplotlib()
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "plots").mkdir(exist_ok=True)

    print("Reading EIGENSTRAT sample and SNP files...")
    samples = read_ind(args.ind)
    snp = read_snp(args.snp)
    n_samples = len(samples)

    print(f"Loading chr11 genotypes: {len(snp)} SNPs x {n_samples} samples")
    G = load_genotypes(args.geno, n_samples, snp["original_index"].to_numpy())

    order = np.argsort(snp["pos"].to_numpy())
    snp = snp.iloc[order].reset_index(drop=True)
    G = G[order, :]
    pos = snp["pos"].to_numpy(dtype=np.int64)

    print("Preparing metadata and sample QC...")
    meta = prepare_metadata(args.metadata, args.anno, samples)
    meta.to_csv(outdir / "01.samples_qc.tsv", sep="\t", index=False)

    meta_pass = meta[meta["pass_qc"]].copy()
    meta_pass = meta_pass[meta_pass["group"].isin(GROUP_STYLE)].copy()
    if meta_pass.empty:
        raise RuntimeError("No East/West Eurasian samples remain after QC.")

    print("Calculating target-region sample callability...")
    meta_pass = add_target_any_call_flags(meta_pass, G, pos, args.window_bp)
    meta_pass.to_csv(outdir / "01.samples_qc.with_target_callability.tsv", sep="\t", index=False)

    max_age_for_bins = float(meta_pass["Age_BP"].max())
    max_age_for_rolling = (
        float(args.rolling_max_age_years)
        if args.rolling_max_age_years is not None
        else max_age_for_bins
    )
    bin_sizes = list(args.bin_sizes)
    if not bin_sizes or any(size <= 0 for size in bin_sizes):
        raise ValueError("--bin-sizes must contain positive values")
    print(f"Temporal bin sizes (years): {bin_sizes}")

    starts, ends, _ = make_windows(
        1,
        CHR11_HG19_END,
        args.window_bp,
        args.window_step_bp,
    )
    print(
        f"Generated {len(starts)} chromosome 11 windows "
        f"({args.window_bp} bp; step {args.window_step_bp} bp)"
    )

    all_rows = []
    rng = np.random.default_rng(args.seed)

    print("Running 1-3-kyr binned analyses with adaptive oldest-bin pooling...")
    for bin_size in bin_sizes:
        for group_name, raw_group_meta in meta_pass.groupby("group", sort=False):
            region = REGIONS[GROUP_STYLE[group_name]["region_name"]]
            group_meta = assign_adaptive_bins_for_group(
                raw_group_meta,
                max_age_bp=max_age_for_bins,
                bin_size_years=bin_size,
                low_n_threshold=args.low_n_threshold,
            )
            bin_order = (
                group_meta[["time_label", "bin_start_years", "bin_end_years", "is_tail_bin"]]
                .drop_duplicates()
                .sort_values("bin_start_years", ascending=False)
            )

            for _, bin_row in bin_order.iterrows():
                label = str(bin_row["time_label"])
                subset = group_meta[group_meta["time_label"].astype(str) == label].copy()
                if subset.empty:
                    continue

                summary, _ = analyze_one_time_unit(
                    G=G,
                    pos=pos,
                    starts=starts,
                    ends=ends,
                    samples_meta=subset,
                    group_name=group_name,
                    region=region,
                    time_label=label,
                    plot_code=int(subset["plot_code"].iloc[0]),
                    unit_type="binned",
                    bin_size_years=bin_size,
                    n_boot=args.n_boot,
                    threads=args.threads,
                    seed=int(rng.integers(1, 2**31 - 1)),
                    do_background=True,
                    candidate_step_bp=args.candidate_step_bp,
                    chrom_start=1,
                    chrom_end=CHR11_HG19_END,
                    window_bp=args.window_bp,
                )
                summary["bin_start_years"] = float(subset["bin_start_years"].iloc[0])
                summary["bin_end_years"] = (
                    float(subset["bin_end_years"].iloc[0])
                    if np.isfinite(subset["bin_end_years"].iloc[0])
                    else np.inf
                )
                summary["is_tail_bin"] = int(bool(subset["is_tail_bin"].iloc[0]))
                summary["adaptive_tail_n_threshold"] = args.low_n_threshold
                all_rows.append(summary)

    print("Running 1-kyr rolling analyses at 100-year steps...")
    centers_years = rolling_centers(
        max_age_for_rolling,
        args.rolling_width_years,
        args.rolling_step_years,
    )
    half_width = args.rolling_width_years / 2.0

    for group_name, group_meta in meta_pass.groupby("group", sort=False):
        region = REGIONS[GROUP_STYLE[group_name]["region_name"]]
        for center_years in centers_years:
            low = center_years - half_width
            high = center_years + half_width
            subset = group_meta[
                (group_meta["Age_BP"] >= low) & (group_meta["Age_BP"] < high)
            ].copy()
            if subset.empty:
                continue

            summary, _ = analyze_one_time_unit(
                G=G,
                pos=pos,
                starts=starts,
                ends=ends,
                samples_meta=subset,
                group_name=group_name,
                region=region,
                time_label=f"rolling_{center_years / 1000:.2f}kya",
                plot_code=-float(center_years),
                unit_type="rolling",
                bin_size_years=args.rolling_width_years,
                n_boot=args.n_boot,
                threads=args.threads,
                seed=int(rng.integers(1, 2**31 - 1)),
                do_background=False,
                candidate_step_bp=args.candidate_step_bp,
                chrom_start=1,
                chrom_end=CHR11_HG19_END,
                window_bp=args.window_bp,
            )
            summary["age_center_years"] = float(center_years)
            summary["rolling_width_years"] = args.rolling_width_years
            summary["rolling_step_years"] = args.rolling_step_years
            all_rows.append(summary)

    results = pd.DataFrame(all_rows)
    results.to_csv(outdir / "02.windowed_region_summary.tsv", sep="\t", index=False)

    print("Generating final temporal plots and LOWESS permutation summaries...")
    generate_plots_from_results(results, outdir, args)

    print("Done.")
    print(f"Output directory: {outdir}")
    print("Main outputs:")
    print(f"  {outdir / '01.samples_qc.tsv'}")
    print(f"  {outdir / '01.samples_qc.with_target_callability.tsv'}")
    print(f"  {outdir / '02.windowed_region_summary.tsv'}")
    print(f"  {outdir / '05.plot_ready.windowed_region_summary.tsv'}")
    print(f"  {outdir / '06.plot_fit_summary.tsv'}")
    print(f"  {outdir / '07.plot_config.tsv'}")
    print(f"  {outdir / 'plots'}")


if __name__ == "__main__":
    main()
