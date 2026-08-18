#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
方案 A: 精确 k-mer 集合 + Mash 距离 + UPGMA 聚类
============================================================
对一个 fasta 内的所有序列做 alignment-free 聚类，用于分子家族划分。

流程 (每个 k 独立):
  1. 逐条序列切 canonical k-mer (presence/absence, 跳过含 N 的窗口)
  2. 构建稀疏二值矩阵 (序列 × k-mer)，可选过滤 ubiquitous (>阈值) 的 k-mer
  3. 精确计算两两 Jaccard 相似度
  4. Mash 距离变换:  D = -(1/k) * ln(2J / (1+J))   (分歧度校正, 不饱和)
  5. Average Linkage (UPGMA) 层次聚类 + Silhouette 选最优簇数
  6. PCoA (classical MDS)
  7. 输出: 距离矩阵(.npz) / linkage / Newick(可上 iTOL) / 元数据 /
           聚类表 / PCoA 坐标+图 / 树状图 / iTOL 颜色条 / summary.json

Header 解析 (以 '#' 分隔):
  >ENST00000325558.11#PGA3#hg38               -> gene=PGA3,  species=hg38
  >mm10.ENS...#Pga5#72#Aonyx_cinereus__...    -> gene=Pga5,  species=Aonyx_...
  约定: gene = 第 2 个字段;  species = 最后一个字段。

用法:
  python 01_run_kmer_classification.py \
      --fasta 487.all.toga.PGA_like.local.v3.assign_candidate_ids.nucleotide.fa \
      --outdir output \
      --k_values 13,15,17,19
"""

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram, to_tree
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score, adjusted_rand_score

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- 2-bit 编码表 ---
_CODE = np.full(256, -1, dtype=np.int8)
_CODE[ord("A")] = 0
_CODE[ord("C")] = 1
_CODE[ord("G")] = 2
_CODE[ord("T")] = 3
_CODE[ord("a")] = 0
_CODE[ord("c")] = 1
_CODE[ord("g")] = 2
_CODE[ord("t")] = 3
_COMP = {0: 3, 1: 2, 2: 1, 3: 0}


# ----------------------------------------------------------------------
# fasta 读取 + header 解析
# ----------------------------------------------------------------------
def read_fasta(path):
    """返回 list[(header, seq)]，header 不含 '>'。"""
    records = []
    header = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == ">":
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def parse_header(header):
    """gene = 第二字段, species = 最后字段。"""
    parts = header.split("#")
    gene = parts[1] if len(parts) >= 2 else "NA"
    species = parts[-1] if len(parts) >= 2 else "NA"
    return gene, species


def sanitize(label):
    """Newick / iTOL 安全的 ID。"""
    out = []
    for ch in label:
        if ch in " ()[]:;,'\t#":
            out.append("_")
        else:
            out.append(ch)
    return "".join(out)


# ----------------------------------------------------------------------
# canonical k-mer 提取 (presence/absence, 跳过 N)
# ----------------------------------------------------------------------
def seq_kmer_codes(seq, k):
    """返回该序列出现过的 canonical k-mer 编码集合 (set[int])。"""
    mask = (1 << (2 * k)) - 1
    shift = 2 * (k - 1)
    code = 0
    rc = 0
    length = 0
    kmers = set()
    b = seq.encode("ascii", "replace")
    arr = np.frombuffer(b, dtype=np.uint8)
    vals = _CODE[arr].tolist()  # -1 表示非 ACGT (含 N); list 迭代更快
    for v in vals:
        if v < 0:
            length = 0
            code = 0
            rc = 0
            continue
        code = ((code << 2) | v) & mask
        rc = (rc >> 2) | (_COMP[v] << shift)
        length += 1
        if length >= k:
            kmers.add(code if code < rc else rc)
    return kmers


def build_matrix(records, k, ubiquitous_threshold):
    """构建 序列×k-mer 稀疏二值矩阵 (CSR)。返回 (X, n_total_kmers, n_kept)。"""
    n = len(records)
    rows_list = []
    codes_list = []
    for i, (_, seq) in enumerate(records):
        kset = seq_kmer_codes(seq, k)
        if not kset:
            continue
        codes_arr = np.fromiter(kset, dtype=np.uint64, count=len(kset))
        codes_list.append(codes_arr)
        rows_list.append(np.full(len(kset), i, dtype=np.int32))
        if (i + 1) % 500 == 0:
            print(f"    [k={k}] kmer 提取 {i + 1}/{n}", flush=True)

    all_codes = np.concatenate(codes_list)
    all_rows = np.concatenate(rows_list)
    # 把 k-mer 编码 factorize 成列索引
    uniq, col_idx = np.unique(all_codes, return_inverse=True)
    n_total = uniq.size
    data = np.ones(all_rows.size, dtype=np.float32)
    X = csr_matrix((data, (all_rows, col_idx.astype(np.int64))), shape=(n, n_total))
    X.data[:] = 1.0  # 二值化 (去重已在 set 完成, 保险)

    # 过滤 ubiquitous: 出现在 > 阈值 比例序列中的 k-mer
    col_counts = np.asarray((X > 0).sum(axis=0)).ravel()
    keep = col_counts <= ubiquitous_threshold * n
    n_kept = int(keep.sum())
    if n_kept < n_total:
        X = X[:, keep]
    return X, n_total, n_kept


# ----------------------------------------------------------------------
# Jaccard + Mash 距离
# ----------------------------------------------------------------------
def jaccard_mash_distance(X, k):
    """精确 Jaccard -> Mash 距离。返回 n×n 距离矩阵 (float64)。"""
    Xb = (X > 0).astype(np.float32)
    inter = (Xb @ Xb.T).toarray().astype(np.float64)  # 交集大小
    sizes = np.asarray(Xb.sum(axis=1)).ravel()
    union = sizes[:, None] + sizes[None, :] - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        J = np.where(union > 0, inter / union, 0.0)
    np.fill_diagonal(J, 1.0)

    # Mash 距离:  D = -(1/k) * ln(2J/(1+J))
    with np.errstate(divide="ignore", invalid="ignore"):
        D = -(1.0 / k) * np.log(2.0 * J / (1.0 + J))
    # J=0 -> D=inf, 用有限最大值封顶, 保证可聚类
    finite = D[np.isfinite(D)]
    cap = finite.max() if finite.size else 1.0
    D[~np.isfinite(D)] = cap
    D[D < 0] = 0.0
    np.fill_diagonal(D, 0.0)
    D = (D + D.T) / 2.0  # 强制对称
    return D


# ----------------------------------------------------------------------
# PCoA (classical MDS)
# ----------------------------------------------------------------------
def pcoa(D, n_axes=10):
    n = D.shape[0]
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J.dot(D2).dot(J)
    B = (B + B.T) / 2.0
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    pos = vals > 1e-9
    vals_p = vals[pos]
    vecs_p = vecs[:, pos]
    coords = vecs_p * np.sqrt(vals_p)
    n_axes = min(n_axes, coords.shape[1])
    coords = coords[:, :n_axes]
    explained = vals_p[:n_axes] / vals_p.sum()
    return coords, explained


# ----------------------------------------------------------------------
# Newick (迭代式, 避免递归深度问题)
# ----------------------------------------------------------------------
def linkage_to_newick(Z, labels):
    tree, nodelist = to_tree(Z, rd=True)  # nodelist[i] = id 为 i 的节点
    parts = [None] * len(nodelist)
    for node in nodelist:  # id 升序: 子节点必先于父节点处理
        if node.is_leaf():
            parts[node.id] = sanitize(labels[node.id])
        else:
            l, r = node.left, node.right
            ld = max(node.dist - l.dist, 0.0)
            rd = max(node.dist - r.dist, 0.0)
            parts[node.id] = f"({parts[l.id]}:{ld:.6f},{parts[r.id]}:{rd:.6f})"
    return parts[tree.id] + ";"


# ----------------------------------------------------------------------
# iTOL 颜色条 (按 gene family)
# ----------------------------------------------------------------------
_PALETTE = [
    "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
    "#008080", "#e6beff", "#9a6324", "#800000", "#808000",
    "#000075", "#808080", "#ffd8b1", "#aaffc3", "#000000",
]


def gene_color_map(genes, top_n):
    counts = pd.Series(genes).value_counts()
    top = list(counts.index[:top_n])
    cmap = {g: _PALETTE[i % len(_PALETTE)] for i, g in enumerate(top)}
    return cmap, top


def write_itol_colorstrip(path, ids, genes, cmap, dataset_label="GeneFamily"):
    with open(path, "w") as fh:
        fh.write("DATASET_COLORSTRIP\n")
        fh.write("SEPARATOR TAB\n")
        fh.write(f"DATASET_LABEL\t{dataset_label}\n")
        fh.write("COLOR\t#000000\n")
        fh.write("STRIP_WIDTH\t30\n")
        fh.write("DATA\n")
        for sid, g in zip(ids, genes):
            col = cmap.get(g, "#cccccc")
            fh.write(f"{sid}\t{col}\t{g}\n")


# ----------------------------------------------------------------------
# 单个 k 的完整流程
# ----------------------------------------------------------------------
def run_one_k(records, k, outdir, args, meta):
    t0 = time.time()
    kdir = os.path.join(outdir, f"k{k}")
    os.makedirs(kdir, exist_ok=True)
    print(f"\n=== k = {k} ===", flush=True)

    X, n_total, n_kept = build_matrix(records, k, args.ubiquitous_threshold)
    print(f"    k-mer: 总 {n_total:,} -> 保留 {n_kept:,} "
          f"(过滤 >{args.ubiquitous_threshold:.0%} ubiquitous: {n_total - n_kept:,})", flush=True)

    D = jaccard_mash_distance(X, k)
    np.savez_compressed(os.path.join(kdir, "mash_dist.npz"), D=D)

    # --- 层次聚类 (UPGMA) ---
    condensed = squareform(D, checks=False)
    Z = linkage(condensed, method="average")
    np.save(os.path.join(kdir, "linkage.npy"), Z)

    # Newick (上 iTOL 看)
    newick = linkage_to_newick(Z, meta["id"].tolist())
    with open(os.path.join(kdir, "tree.nwk"), "w") as fh:
        fh.write(newick + "\n")

    # --- Silhouette 选最优簇数 ---
    n = D.shape[0]
    max_c = min(args.max_clusters, n - 1)
    sil_scores = {}
    for c in range(2, max_c + 1):
        labels_c = fcluster(Z, t=c, criterion="maxclust")
        if len(set(labels_c)) < 2:
            continue
        sil_scores[c] = float(silhouette_score(D, labels_c, metric="precomputed"))
    best_c = max(sil_scores, key=sil_scores.get) if sil_scores else 2
    best_labels = fcluster(Z, t=best_c, criterion="maxclust")

    # 额外: 在 "gene family 数" 这个解像度上切, 直接对照分子家族
    n_gene = int(meta["gene"].nunique())
    labels_ngene = fcluster(Z, t=min(n_gene, n - 1), criterion="maxclust")

    # 与 gene family 标签对比 (两个解像度)
    ari_gene = float(adjusted_rand_score(meta["gene"].values, best_labels))
    ari_gene_at_ngene = float(adjusted_rand_score(meta["gene"].values, labels_ngene))

    # --- PCoA ---
    coords, explained = pcoa(D, n_axes=10)
    pcoa_df = meta.copy()
    for j in range(coords.shape[1]):
        pcoa_df[f"PCo{j + 1}"] = coords[:, j]
    pcoa_df["cluster_best"] = best_labels
    pcoa_df.to_csv(os.path.join(kdir, "pcoa_coords.tsv"), sep="\t", index=False)

    # --- 聚类表 ---
    clu = meta.copy()
    clu["cluster_best"] = best_labels
    clu["cluster_at_ngenes"] = labels_ngene
    clu.to_csv(os.path.join(kdir, "cluster_assignments.tsv"), sep="\t", index=False)

    # cluster × gene 交叉表
    ct = pd.crosstab(clu["cluster_best"], clu["gene"])
    ct.to_csv(os.path.join(kdir, "cluster_vs_gene.tsv"), sep="\t")

    # --- 可视化 ---
    cmap, top_genes = gene_color_map(meta["gene"].values, args.top_genes)
    _plot_pcoa(coords, explained, meta["gene"].values, cmap, top_genes,
               os.path.join(kdir, "pcoa_by_gene.pdf"), k)
    _plot_dendrogram(Z, best_c, os.path.join(kdir, "dendrogram.pdf"), k)

    # iTOL 颜色条
    write_itol_colorstrip(os.path.join(kdir, "itol_genefamily_colorstrip.txt"),
                          meta["id"].values, meta["gene"].values, cmap)

    summary = {
        "k": k,
        "n_sequences": int(n),
        "n_kmers_total": int(n_total),
        "n_kmers_kept": int(n_kept),
        "ubiquitous_threshold": args.ubiquitous_threshold,
        "best_n_clusters": int(best_c),
        "best_silhouette": float(sil_scores.get(best_c, float("nan"))),
        "ARI_vs_gene_family": ari_gene,
        "ARI_vs_gene_at_ngenes": ari_gene_at_ngene,
        "n_gene_families": int(meta["gene"].nunique()),
        "n_species": int(meta["species"].nunique()),
        "silhouette_scan": sil_scores,
        "runtime_sec": round(time.time() - t0, 1),
    }
    with open(os.path.join(kdir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(f"    最优簇数={best_c}  Silhouette={summary['best_silhouette']:.3f}  "
          f"ARI(vs gene)={ari_gene:.3f}  用时 {summary['runtime_sec']}s", flush=True)
    return summary


def _plot_pcoa(coords, explained, genes, cmap, top_genes, path, k):
    fig, ax = plt.subplots(figsize=(9, 8))
    genes = np.asarray(genes)
    top_set = set(top_genes)
    # 先画 other (灰)
    other = ~np.isin(genes, list(top_set))
    if other.any():
        ax.scatter(coords[other, 0], coords[other, 1], s=8, c="#dddddd",
                   label="other", linewidths=0)
    for g in top_genes:
        m = genes == g
        ax.scatter(coords[m, 0], coords[m, 1], s=12, c=cmap[g], label=g, linewidths=0)
    ax.set_xlabel(f"PCo1 ({explained[0] * 100:.1f}%)")
    ax.set_ylabel(f"PCo2 ({explained[1] * 100:.1f}%)")
    ax.set_title(f"PCoA (k={k}, Mash distance) — colored by gene family")
    ax.legend(markerscale=2, fontsize=7, ncol=2, loc="best")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_dendrogram(Z, color_threshold_clusters, path, k):
    fig, ax = plt.subplots(figsize=(16, 6))
    # 用第 (n-c) 高的合并高度作为着色阈值
    heights = np.sort(Z[:, 2])
    thr = heights[-color_threshold_clusters] if color_threshold_clusters <= len(heights) else 0
    dendrogram(Z, no_labels=True, color_threshold=thr, ax=ax)
    ax.set_title(f"UPGMA dendrogram (k={k}); colored at {color_threshold_clusters} clusters")
    ax.set_ylabel("Mash distance")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="方案A: k-mer + Mash距离 + UPGMA 聚类")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", default="output")
    ap.add_argument("--k_values", default="13,15,17,19")
    ap.add_argument("--ubiquitous_threshold", type=float, default=0.95,
                    help="出现比例 > 此值的 k-mer 视为无信息并删除 (默认 0.95)")
    ap.add_argument("--max_clusters", type=int, default=60,
                    help="Silhouette 扫描的最大簇数")
    ap.add_argument("--top_genes", type=int, default=15,
                    help="PCoA / iTOL 着色的 top gene family 数")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    k_values = [int(x) for x in args.k_values.split(",") if x.strip()]

    print(f"读取 fasta: {args.fasta}", flush=True)
    records = read_fasta(args.fasta)
    print(f"  序列数: {len(records)}", flush=True)

    genes, species = zip(*[parse_header(h) for h, _ in records])
    ids = [sanitize(h) for h, _ in records]
    meta = pd.DataFrame({"id": ids, "header": [h for h, _ in records],
                         "gene": genes, "species": species})
    meta.to_csv(os.path.join(args.outdir, "metadata.tsv"), sep="\t", index=False)
    print(f"  gene family 数: {meta['gene'].nunique()};  物种数: {meta['species'].nunique()}", flush=True)

    summaries = []
    for k in k_values:
        summaries.append(run_one_k(records, k, args.outdir, args, meta))

    # 汇总表
    summ_df = pd.DataFrame([{
        "k": s["k"], "n_kmers_kept": s["n_kmers_kept"],
        "best_n_clusters": s["best_n_clusters"],
        "best_silhouette": round(s["best_silhouette"], 4),
        "ARI_vs_gene_best": round(s["ARI_vs_gene_family"], 4),
        "ARI_vs_gene_at_ngenes": round(s["ARI_vs_gene_at_ngenes"], 4),
        "runtime_sec": s["runtime_sec"],
    } for s in summaries])
    summ_df.to_csv(os.path.join(args.outdir, "summary_all_k.tsv"), sep="\t", index=False)
    print("\n=== 汇总 ===", flush=True)
    print(summ_df.to_string(index=False), flush=True)
    print(f"\n输出目录: {os.path.abspath(args.outdir)}", flush=True)


if __name__ == "__main__":
    main()
