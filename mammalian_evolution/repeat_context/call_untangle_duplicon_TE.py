#!/usr/bin/env python3
"""
Direction-aware PGA duplicon core / junction TE caller from odgi untangle output.

The biological direction is inferred only from TOGA PGA gene BED strands.  The
untangle strand is used for reference traversal segmentation, not for duplicon
orientation.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


PGA_RE = re.compile(r"^PGA_.")


@dataclass
class Gene:
    short: str
    full_name: str
    chrom: str
    gene: str
    raw_start: int
    raw_end: int
    strand: str
    local_start: int
    local_end: int
    local_mid: float
    bio_order: int = 0
    oriented_start: int = 0
    oriented_end: int = 0


@dataclass
class SpeciesMeta:
    short: str
    full_name: str
    directory: Path | None
    length: int | None
    anchor_start: int | None
    anchor_end: int | None
    pga_strand: str
    genes: list[Gene]
    rm_file: Path | None
    status: str


@dataclass
class Segment:
    line_no: int
    query_name: str
    qstart: int
    qend: int
    ref_name: str
    rstart: int
    rend: int
    score: float
    inv: str
    self_cov: float
    nth: str


@dataclass
class CandidateCore:
    candidate_id: int
    segments: list[Segment]
    qstart: int
    qend: int
    ref_min: int
    ref_max: int
    strand: str
    score_min: float
    score_median: float
    self_cov_median: float
    self_cov_max: float


def parse_args() -> argparse.Namespace:
    default_base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Call direction-aware PGA duplicon cores, transitions, and TE from odgi untangle outputs."
    )
    parser.add_argument("--base-dir", type=Path, default=default_base)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--species-fa", type=Path, default=None)
    parser.add_argument("--untangle-dir", type=Path, default=None)
    parser.add_argument("--te-window", type=int, default=0)
    parser.add_argument("--selfcov-ratio", type=float, default=0.75)
    parser.add_argument("--min-selfcov", type=float, default=1.5)
    parser.add_argument("--ref-backtrack-tol", type=int, default=100)
    parser.add_argument("--min-core-len", type=int, default=1000)
    parser.add_argument("--min-gene-overlap-frac", type=float, default=0.5)
    return parser.parse_args()


def read_graph_species(fasta: Path) -> list[str]:
    names: list[str] = []
    with fasta.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                names.append(line[1:].strip().split()[0])
    return names


def read_fai_lengths(fai: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    if not fai.exists():
        return lengths
    with fai.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
    return lengths


def find_species_dirs(base_dir: Path) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for item in base_dir.iterdir():
        if item.is_dir() and item.name != "apes_owms":
            out[item.name.split("__")[0]] = item
    return out


def read_anchor(anchor_bed: Path) -> tuple[str, int, int]:
    with anchor_bed.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.strip():
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    raise ValueError(f"Invalid anchor BED: {anchor_bed}")
                return fields[0], int(fields[1]), int(fields[2])
    raise ValueError(f"Empty anchor BED: {anchor_bed}")


def find_rm_file(species_dir: Path) -> Path | None:
    rms = sorted(species_dir.glob("*.pga.anchor.locus.fa.out"))
    return rms[0] if rms else None


def oriented_coord(pos: int, strand: str, seq_len: int | None) -> int:
    if strand == "+":
        return pos
    if seq_len is None:
        return -pos
    return seq_len - pos


def load_species_meta(base_dir: Path, species_names: list[str]) -> dict[str, SpeciesMeta]:
    dirs = find_species_dirs(base_dir)
    lengths = read_fai_lengths(base_dir / "apes_owms.fa.fai")
    result: dict[str, SpeciesMeta] = {}
    for short in species_names:
        species_dir = dirs.get(short)
        if species_dir is None:
            result[short] = SpeciesMeta(short, short, None, lengths.get(short), None, None, "NA", [], None, "missing_species_dir")
            continue

        anchor_bed = species_dir / "pga.anchor.locus.bed"
        toga_bed = species_dir / "toga.PGA_like.local.v4.assign_candidate_ids.bed"
        rm_file = find_rm_file(species_dir)
        if not anchor_bed.exists() or not toga_bed.exists():
            status = []
            if not anchor_bed.exists():
                status.append("missing_anchor_bed")
            if not toga_bed.exists():
                status.append("missing_toga_bed")
            result[short] = SpeciesMeta(short, species_dir.name, species_dir, lengths.get(short), None, None, "NA", [], rm_file, ";".join(status))
            continue

        chrom, anchor_start, anchor_end = read_anchor(anchor_bed)
        seq_len = anchor_end - anchor_start
        genes: list[Gene] = []
        with toga_bed.open(encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6 or not PGA_RE.match(fields[3]):
                    continue
                raw_start = int(fields[1])
                raw_end = int(fields[2])
                local_start = max(0, raw_start - anchor_start)
                local_end = min(seq_len, raw_end - anchor_start)
                if local_end <= local_start:
                    continue
                genes.append(
                    Gene(
                        short=short,
                        full_name=species_dir.name,
                        chrom=fields[0],
                        gene=fields[3],
                        raw_start=raw_start,
                        raw_end=raw_end,
                        strand=fields[5],
                        local_start=local_start,
                        local_end=local_end,
                        local_mid=(local_start + local_end) / 2,
                    )
                )
        strands = sorted({gene.strand for gene in genes})
        if not genes:
            pga_strand = "NA"
            status = "no_PGA_gene"
        elif len(strands) > 1:
            pga_strand = "mixed"
            status = "mixed_PGA_strand"
        else:
            pga_strand = strands[0]
            status = "OK"

        if pga_strand in {"+", "-"}:
            genes_for_order = sorted(genes, key=lambda g: g.local_start, reverse=(pga_strand == "-"))
            for idx, gene in enumerate(genes_for_order, start=1):
                gene.bio_order = idx
            for gene in genes:
                gene.oriented_start = oriented_coord(gene.local_start, pga_strand, seq_len)
                gene.oriented_end = oriented_coord(gene.local_end, pga_strand, seq_len)

        result[short] = SpeciesMeta(short, species_dir.name, species_dir, lengths.get(short, seq_len), anchor_start, anchor_end, pga_strand, genes, rm_file, status)
    return result


def parse_untangle_name(path: Path) -> tuple[str, str]:
    name = path.name
    m = re.match(r"^(.+?)\.(\d+)\.untangle\.txt$", name)
    if m:
        return m.group(1), m.group(2)
    m = re.match(r"^(.+?)\.untangle\.txt$", name)
    if m:
        return m.group(1), "NA"
    return name.replace(".untangle.txt", ""), "NA"


def read_untangle(path: Path, allowed_species: set[str]) -> dict[str, list[Segment]]:
    by_query: dict[str, list[Segment]] = defaultdict(list)
    with path.open(encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for line_no, row in enumerate(reader, start=2):
            query = row.get("query.name", "")
            if query not in allowed_species:
                continue
            try:
                seg = Segment(
                    line_no=line_no,
                    query_name=query,
                    qstart=int(row["query.start"]),
                    qend=int(row["query.end"]),
                    ref_name=row["ref.name"],
                    rstart=int(row["ref.start"]),
                    rend=int(row["ref.end"]),
                    score=float(row["score"]),
                    inv=row["inv"],
                    self_cov=float(row["self.cov"]),
                    nth=row.get("n.th", ""),
                )
            except (KeyError, ValueError):
                continue
            by_query[query].append(seg)
    for query in by_query:
        by_query[query].sort(key=lambda s: (s.qstart, s.qend, s.rstart, s.rend))
    return by_query


def ref_backtracked(prev: Segment, cur: Segment, tol: int) -> bool:
    inv = prev.inv if prev.inv in {"+", "-"} else cur.inv
    if inv == "+":
        return cur.rstart < prev.rstart - tol
    if inv == "-":
        return cur.rstart > prev.rstart + tol
    return False


def candidate_cores(
    rows: list[Segment],
    selfcov_threshold: float,
    ref_backtrack_tol: int,
    min_core_len: int,
) -> list[CandidateCore]:
    if not rows:
        return []
    chunks: list[list[Segment]] = []
    cur: list[Segment] = []
    for row in rows:
        if cur and ref_backtracked(cur[-1], row, ref_backtrack_tol):
            chunks.append(cur)
            cur = []
        cur.append(row)
    if cur:
        chunks.append(cur)

    candidates: list[CandidateCore] = []
    for chunk in chunks:
        high_idx = [idx for idx, seg in enumerate(chunk) if seg.self_cov >= selfcov_threshold]
        if not high_idx:
            continue
        trimmed = chunk[min(high_idx) : max(high_idx) + 1]
        qstart = min(seg.qstart for seg in trimmed)
        qend = max(seg.qend for seg in trimmed)
        if qend - qstart < min_core_len:
            continue
        scores = [seg.score for seg in trimmed]
        covs = [seg.self_cov for seg in trimmed]
        ref_min = min(min(seg.rstart, seg.rend) for seg in trimmed)
        ref_max = max(max(seg.rstart, seg.rend) for seg in trimmed)
        strands = [seg.inv for seg in trimmed if seg.inv in {"+", "-"}]
        strand = statistics.mode(strands) if strands else "NA"
        candidates.append(
            CandidateCore(
                candidate_id=len(candidates) + 1,
                segments=trimmed,
                qstart=qstart,
                qend=qend,
                ref_min=ref_min,
                ref_max=ref_max,
                strand=strand,
                score_min=min(scores),
                score_median=statistics.median(scores),
                self_cov_median=statistics.median(covs),
                self_cov_max=max(covs),
            )
        )
    return candidates


def assign_cores_to_genes(
    candidates: list[CandidateCore],
    genes: list[Gene],
    min_overlap_frac: float,
) -> tuple[list[dict], list[CandidateCore]]:
    assigned: list[dict] = []
    used: set[int] = set()
    unassigned = set(range(len(candidates)))
    for gene in sorted(genes, key=lambda g: g.local_start):
        gene_len = gene.local_end - gene.local_start
        best: tuple[float, float, int, CandidateCore, int, bool] | None = None
        for idx, cand in enumerate(candidates):
            if idx in used:
                continue
            overlap = max(0, min(cand.qend, gene.local_end) - max(cand.qstart, gene.local_start))
            overlap_frac = overlap / gene_len if gene_len else 0
            midpoint_hit = cand.qstart <= gene.local_mid <= cand.qend
            if not midpoint_hit and overlap_frac < min_overlap_frac:
                continue
            distance = abs(((cand.qstart + cand.qend) / 2) - gene.local_mid)
            key = (overlap_frac, overlap, -distance)
            if best is None or key > best[:3]:
                best = (overlap_frac, overlap, -distance, cand, idx, midpoint_hit)
        if best is None:
            assigned.append({"gene": gene, "candidate": None, "status": "missing_core"})
            continue
        overlap_frac, overlap, _neg_distance, cand, idx, midpoint_hit = best
        used.add(idx)
        unassigned.discard(idx)
        assigned.append(
            {
                "gene": gene,
                "candidate": cand,
                "status": "OK",
                "gene_overlap_bp": int(overlap),
                "gene_overlap_frac": overlap_frac,
                "gene_midpoint_in_core": midpoint_hit,
            }
        )
    return assigned, [candidates[idx] for idx in sorted(unassigned)]


def parse_repeatmasker(path: Path | None) -> list[dict]:
    if path is None or not path.exists():
        return []
    repeats: list[dict] = []
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip():
                continue
            stripped = line.strip()
            if stripped.startswith("SW") or stripped.startswith("score") or stripped.startswith("There were"):
                continue
            fields = stripped.split()
            if len(fields) < 14:
                continue
            try:
                sw_score = int(fields[0])
                perc_div = float(fields[1])
                seqid = fields[4]
                start_1 = int(fields[5])
                end_1 = int(fields[6])
                strand = fields[8]
                repeat_name = fields[9]
                repeat_class = fields[10]
                rm_id = fields[-1].rstrip("*")
            except ValueError:
                continue
            repeats.append(
                {
                    "seqid": seqid,
                    "start_1based": start_1,
                    "end_1based": end_1,
                    "start0": start_1 - 1,
                    "end0": end_1,
                    "strand": strand,
                    "repeat_name": repeat_name,
                    "repeat_class": repeat_class,
                    "repeat_group": simplify_te_class(repeat_name, repeat_class),
                    "sw_score": sw_score,
                    "perc_div": perc_div,
                    "rm_id": rm_id,
                }
            )
    return repeats


def simplify_te_class(name: str, repeat_class: str) -> str:
    rc = repeat_class.upper()
    nm = name.upper()
    if "ALU" in nm or rc.startswith("SINE/ALU"):
        return "SINE/Alu"
    if rc.startswith("LINE/L1"):
        return "LINE/L1"
    if rc.startswith("LTR") or "ERV" in rc:
        return "LTR/ERV"
    if rc.startswith("DNA"):
        return "DNA"
    if rc.startswith("SINE"):
        return "Other SINE"
    if rc.startswith("LINE"):
        return "Other LINE"
    return "Other TE"


def endpoint_1based(raw_boundary: int, boundary_kind: str) -> int:
    if boundary_kind == "raw_start":
        return raw_boundary + 1
    if boundary_kind == "raw_end":
        return raw_boundary
    raise ValueError(f"Unknown boundary kind: {boundary_kind}")


def endpoint_hits(repeats: list[dict], endpoint: int, boundary_kind: str, window: int) -> list[dict]:
    pos_1 = endpoint_1based(endpoint, boundary_kind)
    hits = []
    for rep in repeats:
        if rep["start_1based"] - window <= pos_1 <= rep["end_1based"] + window:
            dist = 0
            if pos_1 < rep["start_1based"]:
                dist = rep["start_1based"] - pos_1
            elif pos_1 > rep["end_1based"]:
                dist = pos_1 - rep["end_1based"]
            hit = dict(rep)
            hit["endpoint_1based"] = pos_1
            hit["dist_to_endpoint"] = dist
            hits.append(hit)
    return hits


def interval_hits(repeats: list[dict], raw_start: int, raw_end: int) -> list[dict]:
    if raw_end < raw_start:
        raw_start, raw_end = raw_end, raw_start
    start_1 = raw_start + 1
    end_1 = raw_end
    hits = []
    for rep in repeats:
        ov_start = max(start_1, rep["start_1based"])
        ov_end = min(end_1, rep["end_1based"])
        if ov_start <= ov_end:
            hit = dict(rep)
            hit["overlap_bp"] = ov_end - ov_start + 1
            hit["interval_start_1based"] = start_1
            hit["interval_end_1based"] = end_1
            hits.append(hit)
    return hits


def boundary_for_role(raw_start: int, raw_end: int, pga_strand: str, role: str) -> tuple[int, str]:
    if pga_strand == "+":
        return (raw_start, "raw_start") if role == "biological_start" else (raw_end, "raw_end")
    return (raw_end, "raw_end") if role == "biological_start" else (raw_start, "raw_start")


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def add_endpoint_te_rows(
    out_rows: list[dict],
    repeats: list[dict],
    base: dict,
    raw_start: int,
    raw_end: int,
    pga_strand: str,
    window: int,
) -> None:
    for role in ("biological_start", "biological_end"):
        endpoint, boundary_kind = boundary_for_role(raw_start, raw_end, pga_strand, role)
        hits = endpoint_hits(repeats, endpoint, boundary_kind, window)
        if not hits:
            row = dict(base)
            row.update(
                {
                    "endpoint_role": role,
                    "raw_boundary": boundary_kind,
                    "endpoint_raw": endpoint,
                    "endpoint_1based": endpoint_1based(endpoint, boundary_kind),
                    "te_hit": 0,
                }
            )
            out_rows.append(row)
            continue
        for hit in hits:
            row = dict(base)
            row.update(
                {
                    "endpoint_role": role,
                    "raw_boundary": boundary_kind,
                    "endpoint_raw": endpoint,
                    "endpoint_1based": hit["endpoint_1based"],
                    "te_hit": 1,
                    "te_name": hit["repeat_name"],
                    "te_class": hit["repeat_class"],
                    "te_group": hit["repeat_group"],
                    "te_start_1based": hit["start_1based"],
                    "te_end_1based": hit["end_1based"],
                    "te_start0": hit["start0"],
                    "te_end0": hit["end0"],
                    "rm_id": hit["rm_id"],
                    "dist_to_endpoint": hit["dist_to_endpoint"],
                }
            )
            out_rows.append(row)


def median_or_empty(vals: Iterable[float]) -> str:
    vals = list(vals)
    return "" if not vals else f"{statistics.median(vals):.6g}"


def main() -> None:
    args = parse_args()
    base_dir: Path = args.base_dir
    species_fa = args.species_fa or (base_dir / "apes_owms.fa")
    untangle_dir = args.untangle_dir or (base_dir / "apes_owms")
    output_dir = args.output_dir or (base_dir / "untangle_duplicon_TE")
    species_names = read_graph_species(species_fa)
    species_set = set(species_names)
    species_meta = load_species_meta(base_dir, species_names)
    repeat_cache = {short: parse_repeatmasker(meta.rm_file) for short, meta in species_meta.items()}

    gene_rows: list[dict] = []
    for short in species_names:
        meta = species_meta[short]
        for gene in sorted(meta.genes, key=lambda g: g.local_start):
            gene_rows.append(
                {
                    "species": short,
                    "full_name": meta.full_name,
                    "pga_cn": len(meta.genes),
                    "pga_strand": meta.pga_strand,
                    "gene": gene.gene,
                    "bio_order": gene.bio_order,
                    "chrom": gene.chrom,
                    "raw_genomic_start": gene.raw_start,
                    "raw_genomic_end": gene.raw_end,
                    "local_start": gene.local_start,
                    "local_end": gene.local_end,
                    "local_mid": f"{gene.local_mid:.3f}",
                    "oriented_start": gene.oriented_start,
                    "oriented_end": gene.oriented_end,
                    "anchor_start": meta.anchor_start,
                    "anchor_end": meta.anchor_end,
                    "status": meta.status,
                }
            )

    core_rows: list[dict] = []
    transition_rows: list[dict] = []
    duplicon_rows: list[dict] = []
    core_te_rows: list[dict] = []
    transition_te_rows: list[dict] = []
    duplicon_te_rows: list[dict] = []
    summary_rows: list[dict] = []

    untangle_files = sorted(untangle_dir.glob("*untangle.txt"))
    for untangle_file in untangle_files:
        reference_name, m_value = parse_untangle_name(untangle_file)
        by_query = read_untangle(untangle_file, species_set)
        for short in species_names:
            meta = species_meta[short]
            rows = by_query.get(short, [])
            pga_cn = len(meta.genes)
            summary = {
                "untangle_file": untangle_file.name,
                "reference_name": reference_name,
                "m": m_value,
                "query_species": short,
                "pga_strand": meta.pga_strand,
                "expected_pga_cn": pga_cn,
                "untangle_segments": len(rows),
                "called_core_count": 0,
                "called_duplicon_count": 0,
                "called_transition_count": 0,
                "status": meta.status,
            }
            if meta.status != "OK" or pga_cn < 2:
                if pga_cn < 2 and meta.status == "OK":
                    summary["status"] = "single_copy_or_no_dup_call"
                summary_rows.append(summary)
                continue
            if not rows:
                summary["status"] = "no_untangle_rows"
                summary_rows.append(summary)
                continue

            threshold = max(args.min_selfcov, pga_cn * args.selfcov_ratio)
            candidates = candidate_cores(rows, threshold, args.ref_backtrack_tol, args.min_core_len)
            assigned, unassigned = assign_cores_to_genes(candidates, meta.genes, args.min_gene_overlap_frac)
            called = []
            for item in assigned:
                gene = item["gene"]
                cand = item.get("candidate")
                if cand is None:
                    continue
                oriented_start = oriented_coord(cand.qstart, meta.pga_strand, meta.length)
                oriented_end = oriented_coord(cand.qend, meta.pga_strand, meta.length)
                support_ids = ",".join(str(seg.line_no) for seg in cand.segments)
                original_lines = ",".join(f"{seg.qstart}-{seg.qend}:{seg.rstart}-{seg.rend}" for seg in cand.segments)
                core_id = f"{reference_name}|m{m_value}|{short}|{gene.gene}"
                row = {
                    "core_id": core_id,
                    "untangle_file": untangle_file.name,
                    "reference_name": reference_name,
                    "m": m_value,
                    "query_species": short,
                    "gene": gene.gene,
                    "gene_bio_order": gene.bio_order,
                    "pga_strand": meta.pga_strand,
                    "expected_pga_cn": pga_cn,
                    "selfcov_threshold": f"{threshold:.6g}",
                    "core_raw_start": cand.qstart,
                    "core_raw_end": cand.qend,
                    "core_oriented_start": oriented_start,
                    "core_oriented_end": oriented_end,
                    "core_len": cand.qend - cand.qstart,
                    "ref_min": cand.ref_min,
                    "ref_max": cand.ref_max,
                    "untangle_inv": cand.strand,
                    "n_segments": len(cand.segments),
                    "segment_line_ids": support_ids,
                    "segment_intervals": original_lines,
                    "score_min": f"{cand.score_min:.6g}",
                    "score_median": f"{cand.score_median:.6g}",
                    "self_cov_median": f"{cand.self_cov_median:.6g}",
                    "self_cov_max": f"{cand.self_cov_max:.6g}",
                    "gene_local_start": gene.local_start,
                    "gene_local_end": gene.local_end,
                    "gene_overlap_bp": item.get("gene_overlap_bp", ""),
                    "gene_overlap_frac": f"{item.get('gene_overlap_frac', 0):.6g}",
                    "gene_midpoint_in_core": item.get("gene_midpoint_in_core", ""),
                    "status": item.get("status", "OK"),
                }
                core_rows.append(row)
                called.append({"gene": gene, "candidate": cand, "row": row})
                add_endpoint_te_rows(
                    core_te_rows,
                    repeat_cache.get(short, []),
                    {
                        "core_id": core_id,
                        "untangle_file": untangle_file.name,
                        "reference_name": reference_name,
                        "m": m_value,
                        "query_species": short,
                        "gene": gene.gene,
                        "pga_strand": meta.pga_strand,
                    },
                    cand.qstart,
                    cand.qend,
                    meta.pga_strand,
                    args.te_window,
                )

            called.sort(key=lambda x: x["candidate"].qstart, reverse=(meta.pga_strand == "-"))
            for bio_idx, item in enumerate(called, start=1):
                item["bio_idx"] = bio_idx
            for idx, item in enumerate(called):
                cand = item["candidate"]
                gene = item["gene"]
                next_item = called[idx + 1] if idx + 1 < len(called) else None
                transition_id = ""
                transition_raw_start = cand.qend
                transition_raw_end = cand.qend
                transition_len = 0
                transition_type = "terminal_no_following_core"
                if next_item is not None:
                    next_cand = next_item["candidate"]
                    transition_id = f"{reference_name}|m{m_value}|{short}|{gene.gene}|to|{next_item['gene'].gene}"
                    if meta.pga_strand == "+":
                        transition_raw_start = cand.qend
                        transition_raw_end = next_cand.qstart
                    else:
                        transition_raw_start = next_cand.qend
                        transition_raw_end = cand.qstart
                    if transition_raw_end < transition_raw_start:
                        transition_raw_start, transition_raw_end = transition_raw_end, transition_raw_start
                        transition_type = "transition_overlap_or_reversed"
                    else:
                        transition_type = "transition_gap" if transition_raw_end > transition_raw_start else "transition_adjacent"
                    transition_len = transition_raw_end - transition_raw_start
                    transition_oriented_start = oriented_coord(transition_raw_start, meta.pga_strand, meta.length)
                    transition_oriented_end = oriented_coord(transition_raw_end, meta.pga_strand, meta.length)
                    trow = {
                        "transition_id": transition_id,
                        "untangle_file": untangle_file.name,
                        "reference_name": reference_name,
                        "m": m_value,
                        "query_species": short,
                        "pga_strand": meta.pga_strand,
                        "upstream_core_gene": gene.gene,
                        "downstream_core_gene": next_item["gene"].gene,
                        "upstream_core_bio_order": item["bio_idx"],
                        "transition_type": transition_type,
                        "transition_raw_start": transition_raw_start,
                        "transition_raw_end": transition_raw_end,
                        "transition_oriented_start": transition_oriented_start,
                        "transition_oriented_end": transition_oriented_end,
                        "transition_len": transition_len,
                    }
                    transition_rows.append(trow)
                    hits = interval_hits(repeat_cache.get(short, []), transition_raw_start, transition_raw_end)
                    if not hits:
                        transition_te_rows.append({**trow, "te_hit": 0})
                    else:
                        for hit in hits:
                            transition_te_rows.append(
                                {
                                    **trow,
                                    "te_hit": 1,
                                    "te_name": hit["repeat_name"],
                                    "te_class": hit["repeat_class"],
                                    "te_group": hit["repeat_group"],
                                    "te_start_1based": hit["start_1based"],
                                    "te_end_1based": hit["end_1based"],
                                    "te_start0": hit["start0"],
                                    "te_end0": hit["end0"],
                                    "rm_id": hit["rm_id"],
                                    "overlap_bp": hit["overlap_bp"],
                                }
                            )

                if next_item is not None and transition_len >= 0:
                    if meta.pga_strand == "+":
                        dup_start = cand.qstart
                        dup_end = transition_raw_end
                    else:
                        dup_start = transition_raw_start
                        dup_end = cand.qend
                else:
                    dup_start = cand.qstart
                    dup_end = cand.qend
                if dup_end < dup_start:
                    dup_start, dup_end = dup_end, dup_start
                duplicon_id = f"{reference_name}|m{m_value}|{short}|{gene.gene}|duplicon"
                drow = {
                    "duplicon_id": duplicon_id,
                    "untangle_file": untangle_file.name,
                    "reference_name": reference_name,
                    "m": m_value,
                    "query_species": short,
                    "gene": gene.gene,
                    "core_id": item["row"]["core_id"],
                    "transition_id": transition_id,
                    "pga_strand": meta.pga_strand,
                    "bio_order": item["bio_idx"],
                    "duplicon_raw_start": dup_start,
                    "duplicon_raw_end": dup_end,
                    "duplicon_oriented_start": oriented_coord(dup_start, meta.pga_strand, meta.length),
                    "duplicon_oriented_end": oriented_coord(dup_end, meta.pga_strand, meta.length),
                    "duplicon_len": dup_end - dup_start,
                    "core_raw_start": cand.qstart,
                    "core_raw_end": cand.qend,
                    "transition_raw_start": transition_raw_start if transition_id else "",
                    "transition_raw_end": transition_raw_end if transition_id else "",
                    "transition_len": transition_len if transition_id else "",
                }
                duplicon_rows.append(drow)
                add_endpoint_te_rows(
                    duplicon_te_rows,
                    repeat_cache.get(short, []),
                    {
                        "duplicon_id": duplicon_id,
                        "untangle_file": untangle_file.name,
                        "reference_name": reference_name,
                        "m": m_value,
                        "query_species": short,
                        "gene": gene.gene,
                        "pga_strand": meta.pga_strand,
                    },
                    dup_start,
                    dup_end,
                    meta.pga_strand,
                    args.te_window,
                )

            summary["called_core_count"] = len(called)
            summary["called_duplicon_count"] = len(called)
            summary["called_transition_count"] = max(0, len(called) - 1)
            if len(called) != pga_cn:
                summary["status"] = f"core_count_mismatch:{len(called)}"
            summary_rows.append(summary)

    consensus_rows: list[dict] = []
    grouped: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    dup_by_core: dict[str, dict] = {row["core_id"]: row for row in duplicon_rows}
    for row in core_rows:
        grouped[(row["reference_name"], row["query_species"], row["gene"])].append(row)
    for (reference_name, query_species, gene), rows in sorted(grouped.items()):
        starts = [int(row["core_raw_start"]) for row in rows]
        ends = [int(row["core_raw_end"]) for row in rows]
        dstarts = []
        dends = []
        for row in rows:
            dup = dup_by_core.get(row["core_id"])
            if dup:
                dstarts.append(int(dup["duplicon_raw_start"]))
                dends.append(int(dup["duplicon_raw_end"]))
        consensus_rows.append(
            {
                "reference_name": reference_name,
                "query_species": query_species,
                "gene": gene,
                "n_m_supported": len(rows),
                "m_values": ",".join(str(row["m"]) for row in rows),
                "core_start_min": min(starts),
                "core_start_max": max(starts),
                "core_end_min": min(ends),
                "core_end_max": max(ends),
                "core_start_range": max(starts) - min(starts),
                "core_end_range": max(ends) - min(ends),
                "duplicon_start_min": min(dstarts) if dstarts else "",
                "duplicon_start_max": max(dstarts) if dstarts else "",
                "duplicon_end_min": min(dends) if dends else "",
                "duplicon_end_max": max(dends) if dends else "",
            }
        )

    write_tsv(
        output_dir / "untangle_PGA_genes.local.tsv",
        gene_rows,
        [
            "species",
            "full_name",
            "pga_cn",
            "pga_strand",
            "gene",
            "bio_order",
            "chrom",
            "raw_genomic_start",
            "raw_genomic_end",
            "local_start",
            "local_end",
            "local_mid",
            "oriented_start",
            "oriented_end",
            "anchor_start",
            "anchor_end",
            "status",
        ],
    )
    write_tsv(
        output_dir / "untangle_core_intervals.tsv",
        core_rows,
        [
            "core_id",
            "untangle_file",
            "reference_name",
            "m",
            "query_species",
            "gene",
            "gene_bio_order",
            "pga_strand",
            "expected_pga_cn",
            "selfcov_threshold",
            "core_raw_start",
            "core_raw_end",
            "core_oriented_start",
            "core_oriented_end",
            "core_len",
            "ref_min",
            "ref_max",
            "untangle_inv",
            "n_segments",
            "segment_line_ids",
            "segment_intervals",
            "score_min",
            "score_median",
            "self_cov_median",
            "self_cov_max",
            "gene_local_start",
            "gene_local_end",
            "gene_overlap_bp",
            "gene_overlap_frac",
            "gene_midpoint_in_core",
            "status",
        ],
    )
    write_tsv(
        output_dir / "untangle_transition_regions.tsv",
        transition_rows,
        [
            "transition_id",
            "untangle_file",
            "reference_name",
            "m",
            "query_species",
            "pga_strand",
            "upstream_core_gene",
            "downstream_core_gene",
            "upstream_core_bio_order",
            "transition_type",
            "transition_raw_start",
            "transition_raw_end",
            "transition_oriented_start",
            "transition_oriented_end",
            "transition_len",
        ],
    )
    write_tsv(
        output_dir / "untangle_duplicon_intervals.tsv",
        duplicon_rows,
        [
            "duplicon_id",
            "untangle_file",
            "reference_name",
            "m",
            "query_species",
            "gene",
            "core_id",
            "transition_id",
            "pga_strand",
            "bio_order",
            "duplicon_raw_start",
            "duplicon_raw_end",
            "duplicon_oriented_start",
            "duplicon_oriented_end",
            "duplicon_len",
            "core_raw_start",
            "core_raw_end",
            "transition_raw_start",
            "transition_raw_end",
            "transition_len",
        ],
    )
    endpoint_fields = [
        "untangle_file",
        "reference_name",
        "m",
        "query_species",
        "gene",
        "pga_strand",
        "core_id",
        "duplicon_id",
        "endpoint_role",
        "raw_boundary",
        "endpoint_raw",
        "endpoint_1based",
        "te_hit",
        "te_name",
        "te_class",
        "te_group",
        "te_start_1based",
        "te_end_1based",
        "te_start0",
        "te_end0",
        "rm_id",
        "dist_to_endpoint",
    ]
    write_tsv(output_dir / "untangle_core_endpoint_TE.tsv", core_te_rows, endpoint_fields)
    write_tsv(
        output_dir / "untangle_transition_region_TE.tsv",
        transition_te_rows,
        [
            "transition_id",
            "untangle_file",
            "reference_name",
            "m",
            "query_species",
            "pga_strand",
            "upstream_core_gene",
            "downstream_core_gene",
            "upstream_core_bio_order",
            "transition_type",
            "transition_raw_start",
            "transition_raw_end",
            "transition_oriented_start",
            "transition_oriented_end",
            "transition_len",
            "te_hit",
            "te_name",
            "te_class",
            "te_group",
            "te_start_1based",
            "te_end_1based",
            "te_start0",
            "te_end0",
            "rm_id",
            "overlap_bp",
        ],
    )
    write_tsv(output_dir / "untangle_duplicon_endpoint_TE.tsv", duplicon_te_rows, endpoint_fields)
    write_tsv(
        output_dir / "untangle_m_consensus.tsv",
        consensus_rows,
        [
            "reference_name",
            "query_species",
            "gene",
            "n_m_supported",
            "m_values",
            "core_start_min",
            "core_start_max",
            "core_end_min",
            "core_end_max",
            "core_start_range",
            "core_end_range",
            "duplicon_start_min",
            "duplicon_start_max",
            "duplicon_end_min",
            "duplicon_end_max",
        ],
    )
    write_tsv(
        output_dir / "untangle_run_summary.tsv",
        summary_rows,
        [
            "untangle_file",
            "reference_name",
            "m",
            "query_species",
            "pga_strand",
            "expected_pga_cn",
            "untangle_segments",
            "called_core_count",
            "called_duplicon_count",
            "called_transition_count",
            "status",
        ],
    )
    print(f"Wrote output directory: {output_dir}")
    print(f"Species: {len(species_names)}")
    print(f"Untangle files: {len(untangle_files)}")
    print(f"Core rows: {len(core_rows)}")
    print(f"Duplicon rows: {len(duplicon_rows)}")


if __name__ == "__main__":
    main()
