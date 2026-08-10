#!/usr/bin/env python3

"""Association between regional mean diploid PGA copy number and plant-derived protein fraction."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import permutation_test, rankdata, spearmanr, t


PC_COLUMNS = ["PC1", "PC2", "PC3"]

# 1KGP population-to-FAOSTAT mapping used in the final analysis.
# ACB and ASW have mappings but are excluded because New_superpop is undefined.
# PUR is intentionally left unmapped.
POP_TO_REGION = {
    "ACB": "Barbados",
    "ASW": "USA",
    "BEB": "Bangladesh",
    "CEU": "USA",
    "CHB": "China, mainland",
    "CHS": "China, mainland",
    "CDX": "China, mainland",
    "CLM": "Colombia",
    "ESN": "Nigeria",
    "FIN": "Finland",
    "GBR": "UK",
    "GIH": "India",
    "GWD": "Gambia",
    "IBS": "Spain",
    "ITU": "India",
    "JPT": "Japan",
    "KHV": "Viet Nam",
    "LWK": "Kenya",
    "MSL": "Sierra Leone",
    "MXL": "Mexico",
    "PEL": "Peru",
    "PJL": "Pakistan",
    "STU": "Sri Lanka",
    "TSI": "Italy",
    "YRI": "Nigeria",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Test the regional association between mean diploid PGA CN and "
            "plant-derived protein fraction in the 1KGP cohort."
        )
    )
    parser.add_argument("--cnv", type=Path, required=True, help="1KGP/HGDP PGA CN table")
    parser.add_argument("--metadata", type=Path, required=True, help="1KGP metadata table")
    parser.add_argument("--pca", type=Path, required=True, help="Joint diploid PCA table")
    parser.add_argument("--fao", type=Path, required=True, help="Prepared FAOSTAT baseline CSV")
    parser.add_argument("--out-prefix", type=Path, required=True, help="Output prefix")
    parser.add_argument("--n-permutations", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def require_columns(df: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = sorted(set(columns) - set(df.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def require_unique(df: pd.DataFrame, column: str, label: str) -> None:
    duplicated = df[column].duplicated(keep=False)
    if duplicated.any():
        examples = df.loc[duplicated, column].astype(str).head(10).tolist()
        raise ValueError(f"{label} has duplicate {column} values: {examples}")


def load_cnv(path: Path) -> pd.DataFrame:
    cn = pd.read_csv(path, sep="\t", dtype={"Sample": str})
    require_columns(cn, ["Sample", "PGA34", "PGA5", "Population"], "CN table")
    cn = cn.loc[~cn["Sample"].str.startswith("HGDP", na=False)].copy()
    require_unique(cn, "Sample", "1KGP CN table")

    cn["PGA34"] = pd.to_numeric(cn["PGA34"], errors="coerce")
    cn["PGA5"] = pd.to_numeric(cn["PGA5"], errors="coerce")
    cn["Diploid_PGA_CN"] = cn["PGA34"] + cn["PGA5"]
    return cn[["Sample", "Population", "Diploid_PGA_CN"]].copy()


def load_metadata(path: Path) -> pd.DataFrame:
    metadata = pd.read_csv(path, sep="\t", dtype={"Sample": str})
    require_columns(metadata, ["Sample", "Population", "New_superpop"], "1KGP metadata")
    require_unique(metadata, "Sample", "1KGP metadata")
    return metadata[["Sample", "Population", "New_superpop"]].copy()


def load_pca(path: Path) -> pd.DataFrame:
    pca = pd.read_csv(path, sep="\t", dtype={"Basename": str})
    require_columns(pca, ["Basename", "Panel", *PC_COLUMNS], "PCA table")
    pca = pca.loc[pca["Panel"].eq("1KGP panel")].copy()
    if pca.empty:
        raise ValueError("No rows with Panel == '1KGP panel' were found in the PCA table.")
    require_unique(pca, "Basename", "1KGP PCA table")

    for column in PC_COLUMNS:
        pca[column] = pd.to_numeric(pca[column], errors="coerce")

    return pca[["Basename", *PC_COLUMNS]].rename(columns={"Basename": "Sample"})


def load_fao(path: Path) -> pd.DataFrame:
    fao = pd.read_csv(path)
    require_columns(fao, ["region", "plant_ratio"], "FAOSTAT baseline table")
    fao["region"] = fao["region"].astype(str).str.strip()
    require_unique(fao, "region", "FAOSTAT baseline table")
    fao["plant_ratio"] = pd.to_numeric(fao["plant_ratio"], errors="coerce")
    return fao[["region", "plant_ratio"]].copy()


def exclusion_reason(row: pd.Series) -> str:
    if pd.isna(row.get("Diploid_PGA_CN")):
        return "missing_cn"
    if pd.isna(row.get("New_superpop")) or not str(row.get("New_superpop")).strip():
        return "missing_new_superpop"
    if pd.isna(row.get("Population")) or not str(row.get("Population")).strip():
        return "missing_population"
    if pd.isna(row.get("Mapped_Region")) or not str(row.get("Mapped_Region")).strip():
        return "unmapped_population"
    if pd.isna(row.get("Plant_Ratio")):
        return "missing_diet"
    if any(pd.isna(row.get(pc)) for pc in PC_COLUMNS):
        return "missing_pc"
    return "included"


def prepare_samples(
    cn: pd.DataFrame,
    metadata: pd.DataFrame,
    pca: pd.DataFrame,
    fao: pd.DataFrame,
) -> pd.DataFrame:
    # Start from the complete 1KGP PCA panel and attach CN/metadata by sample ID.
    samples = pca.merge(metadata, on="Sample", how="left", validate="one_to_one")
    samples = samples.merge(
        cn,
        on="Sample",
        how="left",
        suffixes=("_meta", "_cn"),
        validate="one_to_one",
    )

    disagreement = (
        samples["Population_meta"].notna()
        & samples["Population_cn"].notna()
        & samples["Population_meta"].ne(samples["Population_cn"])
    )
    if disagreement.any():
        examples = samples.loc[
            disagreement, ["Sample", "Population_meta", "Population_cn"]
        ].head(10)
        raise ValueError(
            "1KGP CN and metadata population labels disagree: "
            f"{examples.to_dict('records')}"
        )

    samples["Population"] = samples["Population_meta"].fillna(samples["Population_cn"])
    samples = samples.drop(columns=["Population_meta", "Population_cn"])

    samples["Mapped_Region"] = samples["Population"].map(POP_TO_REGION)
    diet_lookup = dict(zip(fao["region"], fao["plant_ratio"]))
    samples["Plant_Ratio"] = samples["Mapped_Region"].map(diet_lookup)
    samples["Status"] = samples.apply(exclusion_reason, axis=1)
    return samples


def join_unique(values: pd.Series) -> str:
    return ";".join(sorted(values.dropna().astype(str).unique()))


def aggregate_regions(samples: pd.DataFrame) -> pd.DataFrame:
    included = samples.loc[samples["Status"].eq("included")].copy()
    if included.empty:
        raise ValueError("No samples remain after filtering.")

    regions = (
        included.groupby("Mapped_Region", as_index=False, sort=True)
        .agg(
            Populations=("Population", join_unique),
            Sample_Count=("Sample", "nunique"),
            Plant_Ratio=("Plant_Ratio", "first"),
            Mean_Diploid_PGA_CN=("Diploid_PGA_CN", "mean"),
            Mean_PC1=("PC1", "mean"),
            Mean_PC2=("PC2", "mean"),
            Mean_PC3=("PC3", "mean"),
        )
        .sort_values("Mapped_Region")
        .reset_index(drop=True)
    )
    return regions


def pearson_correlation(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x - x.mean()
    y = y - y.mean()
    denominator = np.linalg.norm(x) * np.linalg.norm(y)
    if denominator <= np.finfo(float).eps:
        raise ValueError("Correlation is undefined because one variable is constant.")
    return float(np.dot(x, y) / denominator)


def permutation_pvalue(
    x: np.ndarray,
    y: np.ndarray,
    observed: float,
    n_permutations: int,
    seed: int,
    residual_maker: np.ndarray | None = None,
    chunk_size: int = 5_000,
) -> float:
    rng = np.random.default_rng(seed)
    x_centered = np.asarray(x, dtype=float) - np.mean(x)
    y = np.asarray(y, dtype=float)
    x_norm = np.linalg.norm(x_centered)
    if x_norm <= np.finfo(float).eps:
        raise ValueError("Permutation predictor is constant.")

    extreme = 0
    completed = 0
    tolerance = 1e-12
    while completed < n_permutations:
        current = min(chunk_size, n_permutations - completed)
        permuted = np.stack([rng.permutation(y) for _ in range(current)])
        if residual_maker is not None:
            permuted = permuted @ residual_maker.T
        permuted = permuted - permuted.mean(axis=1, keepdims=True)
        norms = np.linalg.norm(permuted, axis=1)
        correlations = np.divide(
            permuted @ x_centered,
            norms * x_norm,
            out=np.full(current, np.nan),
            where=norms > np.finfo(float).eps,
        )
        extreme += int(np.sum(np.abs(correlations) >= abs(observed) - tolerance))
        completed += current

    return float((extreme + 1) / (n_permutations + 1))


def raw_spearman(
    plant_ratio: np.ndarray,
    mean_cn: np.ndarray,
    n_permutations: int,
    seed: int,
) -> tuple[float, float]:
    x = np.asarray(plant_ratio, dtype=float)
    y = np.asarray(mean_cn, dtype=float)
    rho = float(spearmanr(x, y).statistic)

    def statistic(x_permuted: np.ndarray) -> float:
        return float(spearmanr(x_permuted, y).statistic)

    result = permutation_test(
        (x,),
        statistic,
        permutation_type="pairings",
        alternative="two-sided",
        n_resamples=n_permutations,
        random_state=seed,
    )
    return rho, float(result.pvalue)


def partial_spearman(
    plant_ratio: np.ndarray,
    mean_cn: np.ndarray,
    pcs: np.ndarray,
    n_permutations: int,
    seed: int,
) -> tuple[float, float, int, int]:
    x_rank = rankdata(np.asarray(plant_ratio, dtype=float), method="average")
    y_rank = rankdata(np.asarray(mean_cn, dtype=float), method="average")
    pcs = np.asarray(pcs, dtype=float)

    n = len(x_rank)
    nuisance = np.column_stack([np.ones(n), pcs])
    nuisance_rank = int(np.linalg.matrix_rank(nuisance))
    if nuisance_rank != nuisance.shape[1]:
        raise ValueError("Regional PC covariate matrix is rank deficient.")

    full_model = np.column_stack([nuisance, x_rank])
    model_rank = int(np.linalg.matrix_rank(full_model))
    residual_df = n - model_rank
    if model_rank != full_model.shape[1] or residual_df <= 0:
        raise ValueError("Partial Spearman model is rank deficient.")

    residual_maker = np.eye(n) - nuisance @ np.linalg.pinv(nuisance)
    x_residual = residual_maker @ x_rank
    y_residual = residual_maker @ y_rank
    rho = pearson_correlation(x_residual, y_residual)

    # Freedman-Lane residual permutation under the reduced model.
    pvalue = permutation_pvalue(
        x_residual,
        y_residual,
        rho,
        n_permutations,
        seed,
        residual_maker=residual_maker,
    )
    return rho, pvalue, model_rank, residual_df


def add_linear_trend(ax: plt.Axes, x: np.ndarray, y: np.ndarray) -> None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3 or np.unique(x).size < 2:
        return

    slope, intercept = np.polyfit(x, y, 1)
    x_grid = np.linspace(x.min(), x.max(), 200)
    y_fit = intercept + slope * x_grid

    fitted = intercept + slope * x
    residuals = y - fitted
    df = len(x) - 2
    if df > 0:
        s2 = np.sum(residuals**2) / df
        ssx = np.sum((x - x.mean()) ** 2)
        if ssx > 0:
            se = np.sqrt(s2 * (1 / len(x) + (x_grid - x.mean()) ** 2 / ssx))
            critical = t.ppf(0.975, df)
            ax.fill_between(x_grid, y_fit - critical * se, y_fit + critical * se, alpha=0.18)

    ax.plot(x_grid, y_fit, linestyle="--", linewidth=1.4)


def format_pvalue(value: float) -> str:
    if value < 0.001:
        return f"{value:.2e}".replace("e-0", "e-").replace("e+0", "e+")
    return f"{value:.4f}"


def plot_association(
    regions: pd.DataFrame,
    rho: float,
    pvalue: float,
    label: str,
    output: Path,
) -> None:
    x = regions["Plant_Ratio"].to_numpy(float)
    y = regions["Mean_Diploid_PGA_CN"].to_numpy(float)

    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    ax.scatter(x, y, s=45, zorder=3)
    add_linear_trend(ax, x, y)

    for _, row in regions.iterrows():
        ax.annotate(
            row["Mapped_Region"],
            (row["Plant_Ratio"], row["Mean_Diploid_PGA_CN"]),
            xytext=(4, 3),
            textcoords="offset points",
            fontsize=8,
        )

    ax.text(
        0.04,
        0.96,
        f"{label} $\\rho$ = {rho:.3f}\n$P_{{perm}}$ = {format_pvalue(pvalue)}",
        transform=ax.transAxes,
        va="top",
        fontsize=10,
    )
    ax.set_xlabel("Plant-derived protein fraction")
    ax.set_ylabel(r"Mean diploid $\mathit{PGA}$ CN")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if args.n_permutations < 1:
        raise ValueError("--n-permutations must be at least 1")

    for path in [args.cnv, args.metadata, args.pca, args.fao]:
        if not path.is_file():
            raise FileNotFoundError(path)

    cn = load_cnv(args.cnv)
    metadata = load_metadata(args.metadata)
    pca = load_pca(args.pca)
    fao = load_fao(args.fao)

    samples = prepare_samples(cn, metadata, pca, fao)
    regions = aggregate_regions(samples)

    plant_ratio = regions["Plant_Ratio"].to_numpy(float)
    mean_cn = regions["Mean_Diploid_PGA_CN"].to_numpy(float)
    pcs = regions[["Mean_PC1", "Mean_PC2", "Mean_PC3"]].to_numpy(float)

    raw_rho, raw_p = raw_spearman(
        plant_ratio,
        mean_cn,
        args.n_permutations,
        args.seed,
    )
    adjusted_rho, adjusted_p, model_rank, residual_df = partial_spearman(
        plant_ratio,
        mean_cn,
        pcs,
        args.n_permutations,
        args.seed,
    )

    prefix = args.out_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)

    sample_path = Path(f"{prefix}.sample_match.tsv")
    region_path = Path(f"{prefix}.region_data.tsv")
    stats_path = Path(f"{prefix}.statistics.tsv")
    raw_figure = Path(f"{prefix}.raw.pdf")
    adjusted_figure = Path(f"{prefix}.pc_adjusted.pdf")

    samples.to_csv(sample_path, sep="\t", index=False, na_rep="NA")
    regions.to_csv(region_path, sep="\t", index=False, na_rep="NA")

    statistics = pd.DataFrame(
        [
            {
                "n_individuals": int(samples["Status"].eq("included").sum()),
                "n_regions": int(len(regions)),
                "n_pcs": 3,
                "model_rank": model_rank,
                "residual_df": residual_df,
                "n_permutations": args.n_permutations,
                "seed": args.seed,
                "raw_spearman_rho": raw_rho,
                "raw_permutation_p": raw_p,
                "pc_adjusted_partial_spearman_rho": adjusted_rho,
                "pc_adjusted_freedman_lane_p": adjusted_p,
            }
        ]
    )
    statistics.to_csv(stats_path, sep="\t", index=False)

    plot_association(regions, raw_rho, raw_p, "Spearman", raw_figure)
    plot_association(
        regions,
        adjusted_rho,
        adjusted_p,
        "PC1-PC3 adjusted",
        adjusted_figure,
    )

    print("Sample inclusion:")
    print(samples["Status"].value_counts(dropna=False).to_string())
    print("\nRegions:", len(regions))
    print("\nStatistics:")
    print(statistics.to_string(index=False))
    print("\nOutputs:")
    for path in [sample_path, region_path, stats_path, raw_figure, adjusted_figure]:
        print(path)


if __name__ == "__main__":
    main()
