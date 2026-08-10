#!/usr/bin/env python3
"""Infer parental origin of phased child haplotypes around the PGA locus."""

from __future__ import annotations

import argparse
import multiprocessing as mp
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

PARENT_LABELS = ("F_hap1", "F_hap2", "M_hap1", "M_hap2")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ped", required=True, type=Path, help="Pedigree TSV with Sample, Father and Mother columns.")
    parser.add_argument("--root", required=True, type=Path, help="Root directory containing one directory per child.")
    parser.add_argument("--out-prefix", required=True, type=Path, help="Prefix for output TSV files.")
    parser.add_argument("--bin-size", type=int, default=10_000)
    parser.add_argument("--flank-size", type=int, default=10_000_000)
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--minimap2", default="minimap2")
    return parser.parse_args()


def load_pedigree(path: Path) -> Dict[str, Dict[str, str]]:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    mother_col = "Mother" if "Mother" in df.columns else "Mather" if "Mather" in df.columns else None
    required = {"Sample", "Father"}
    if mother_col is None or not required.issubset(df.columns):
        raise ValueError("Pedigree file must contain Sample, Father and Mother columns")

    trios = {}
    for _, row in df.iterrows():
        child = row["Sample"]
        father = row["Father"]
        mother = row[mother_col]
        if father not in {"", "0"} and mother not in {"", "0"}:
            trios[child] = {"Father": father, "Mother": mother}
    return trios


def get_pga_cluster_bounds(gff: Path) -> Tuple[int, int]:
    starts, ends = [], []
    with gff.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if not re.search(r"gene_name=PGA[345](?:;|$)", fields[8]):
                continue
            if "Ensembl_canonical" not in fields[8] and "Ensembl canonical" not in fields[8]:
                continue
            starts.append(int(fields[3]))
            ends.append(int(fields[4]))
    if not starts:
        raise ValueError(f"No canonical PGA3/4/5 genes found in {gff}")
    return min(starts), max(ends)


def read_single_fasta(path: Path) -> str:
    seq = []
    seen_header = False
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if seen_header and seq:
                    break
                seen_header = True
                continue
            seq.append(line.strip())
    if not seq:
        raise ValueError(f"No sequence found in {path}")
    return "".join(seq).upper()


def flanking_bins(
    sequence: str,
    pga_start: int,
    pga_end: int,
    flank_size: int,
    bin_size: int,
) -> List[Tuple[int, str]]:
    """Return 1-based bin-start coordinates and sequences, excluding the PGA core."""
    seq_len = len(sequence)
    region_start = max(1, pga_start - flank_size)
    region_end = min(seq_len, pga_end + flank_size)
    bins: List[Tuple[int, str]] = []

    # GFF positions are 1-based inclusive; Python slices are 0-based half-open.
    for start1, end1 in ((region_start, pga_start - 1), (pga_end + 1, region_end)):
        if end1 < start1:
            continue
        for bin_start in range(start1, end1 + 1, bin_size):
            bin_end = min(end1, bin_start + bin_size - 1)
            seq = sequence[bin_start - 1 : bin_end]
            if seq:
                bins.append((bin_start, seq))
    return bins


def write_query_fasta(bins: Iterable[Tuple[int, str]], path: Path) -> None:
    with path.open("w") as out:
        for coord, seq in bins:
            out.write(f">{coord}\n{seq}\n")


def minimap2_identity_scores(minimap2: str, query: Path, target: Path) -> Dict[int, float]:
    cmd = [minimap2, "-t", "2", "-c", "--eqx", "-x", "asm5", str(target), str(query)]
    proc = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

    scores: Dict[int, float] = {}
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        coord = int(fields[0])
        matches = float(fields[9])
        block_len = float(fields[10])
        identity = matches / block_len if block_len > 0 else 0.0
        if coord not in scores or identity > scores[coord]:
            scores[coord] = identity
    return scores


def four_way_origins(coords: List[int], scores: Dict[str, Dict[int, float]]) -> List[str]:
    """Assign the best parental haplotype; carry upstream origin through exact ties."""
    origins: List[str] = []
    previous = "Unknown"

    for coord in coords:
        values = {label: scores[label].get(coord, 0.0) for label in PARENT_LABELS}
        best_score = max(values.values())
        if best_score <= 0:
            origins.append(previous)
            continue

        winners = [label for label, value in values.items() if abs(value - best_score) <= 1e-6]
        if len(winners) == 1:
            previous = winners[0]
        elif previous == "Unknown":
            previous = "Unknown"
        # On a tie, retain the nearest upstream assignment when one exists.
        origins.append(previous)

    return origins


def transmitting_parent(origins: List[str]) -> Tuple[str, int, int]:
    father = sum(origin.startswith("F_") for origin in origins)
    mother = sum(origin.startswith("M_") for origin in origins)
    if father == mother:
        return "Unknown", father, mother
    return ("Father" if father > mother else "Mother"), father, mother


def constrained_origins(
    coords: List[int],
    scores: Dict[str, Dict[int, float]],
    parent: str,
) -> List[str]:
    if parent == "Father":
        hap1, hap2 = "F_hap1", "F_hap2"
    elif parent == "Mother":
        hap1, hap2 = "M_hap1", "M_hap2"
    else:
        return ["Unknown"] * len(coords)

    out = []
    for coord in coords:
        score1 = scores[hap1].get(coord, 0.0)
        score2 = scores[hap2].get(coord, 0.0)
        if abs(score1 - score2) <= 1e-6:
            out.append("Unknown")
        elif score1 > score2:
            out.append(hap1)
        else:
            out.append(hap2)
    return out


def resolve_parent_fasta(root: Path, child_dir: Path, sample: str, hap: str) -> Path:
    candidates = [
        root / sample / f"{sample}.{hap}.fa",
        child_dir / f"{sample}.{hap}.fa",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(f"No FASTA found for {sample}.{hap}; checked: {candidates}")


def process_child(task):
    child, parents, root, bin_size, flank_size, minimap2 = task
    child_dir = root / child
    if not child_dir.is_dir():
        raise FileNotFoundError(child_dir)

    summary_rows = []
    bin_rows = []

    for child_hap in ("hap1", "hap2"):
        gff = child_dir / f"{child}.{child_hap}.liftoff.gff_polished"
        fasta = child_dir / f"{child}.{child_hap}.fa"
        if not gff.exists() or not fasta.exists():
            raise FileNotFoundError(f"Missing child input for {child}.{child_hap}: {gff} or {fasta}")

        pga_start, pga_end = get_pga_cluster_bounds(gff)
        sequence = read_single_fasta(fasta)
        bins = flanking_bins(sequence, pga_start, pga_end, flank_size, bin_size)
        coords = [coord for coord, _ in bins]

        parent_fastas = {
            "F_hap1": resolve_parent_fasta(root, child_dir, parents["Father"], "hap1"),
            "F_hap2": resolve_parent_fasta(root, child_dir, parents["Father"], "hap2"),
            "M_hap1": resolve_parent_fasta(root, child_dir, parents["Mother"], "hap1"),
            "M_hap2": resolve_parent_fasta(root, child_dir, parents["Mother"], "hap2"),
        }

        with tempfile.TemporaryDirectory(prefix=f"{child}.{child_hap}.") as tmp:
            query = Path(tmp) / "bins.fa"
            write_query_fasta(bins, query)
            scores = {
                label: minimap2_identity_scores(minimap2, query, path)
                for label, path in parent_fastas.items()
            }

        four_way = four_way_origins(coords, scores)
        parent, father_count, mother_count = transmitting_parent(four_way)
        two_way = constrained_origins(coords, scores, parent)

        summary_rows.append(
            {
                "Child": child,
                "Haplotype": child_hap,
                "PGA_Start": pga_start,
                "PGA_End": pga_end,
                "Transmitting_Parent": parent,
                "Father_Assigned_Bins": father_count,
                "Mother_Assigned_Bins": mother_count,
                "N_Bins": len(coords),
            }
        )

        for i, coord in enumerate(coords):
            row = {
                "Child": child,
                "Haplotype": child_hap,
                "Bin_Start": coord,
                "Four_Way_Origin": four_way[i],
                "Transmitting_Parent": parent,
                "Parent_Constrained_Origin": two_way[i],
            }
            for label in PARENT_LABELS:
                row[f"{label}_Identity"] = scores[label].get(coord, 0.0)
            bin_rows.append(row)

    return summary_rows, bin_rows


def main() -> None:
    args = parse_args()
    trios = load_pedigree(args.ped)
    tasks = [
        (child, parents, args.root, args.bin_size, args.flank_size, args.minimap2)
        for child, parents in trios.items()
    ]

    n_proc = max(1, min(args.processes, len(tasks))) if tasks else 1
    if n_proc == 1:
        results = [process_child(task) for task in tasks]
    else:
        with mp.Pool(processes=n_proc) as pool:
            results = pool.map(process_child, tasks)

    summaries = [row for summary_rows, _ in results for row in summary_rows]
    bins = [row for _, bin_rows in results for row in bin_rows]

    summary_path = Path(f"{args.out_prefix}.transmission_summary.tsv")
    bins_path = Path(f"{args.out_prefix}.bin_origins.tsv")
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(summaries).to_csv(summary_path, sep="\t", index=False)
    pd.DataFrame(bins).to_csv(bins_path, sep="\t", index=False)

    print(f"Trios: {len(trios)}")
    print(f"Transmission summary: {summary_path}")
    print(f"Bin origins: {bins_path}")


if __name__ == "__main__":
    main()
