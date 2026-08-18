#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
02_refine_toga_local_units.py

Purpose
-------
Refine TOGA annotations inside a VPS37C--VWCE anchor interval into local
PGA/PGA-like/PAG-like gene-structure units.

Compared with v2, v3 adds target-whitelist filtering to avoid counting
non-PGA genes introduced by local rearrangements, which is common in some
rodent assemblies.

Main logic
----------
1. Parse BED12 TOGA models.
2. Classify each model as target/non-target by source transcript/name/gene-label
   whitelist generated from canonical VPS37C--VWCE intervals of reference genomes.
3. Only target models are allowed to become local-unit seeds.
4. Compact complete target models define primary local units.
5. Long-intron target models sharing exons with compact units are marked as
   stretched_projection and are not counted as independent copies.
6. Very long target models are recorded as fused/bridging projections and are not
   counted as independent copies.
7. Non-target models are exported as rearrangement evidence and never enter CN.

Recommended output for CN / assignment
--------------------------------------
- *.local_units.tsv: one row per local target unit.
- *.assign_candidate_ids.txt: complete + partial target local units for downstream
  nucleotide/protein PGA/PAG assignment.
- *.rearrangement_qc.tsv: assembly-level flag for local rearrangement.
"""

import argparse
import gzip
import os
import re
import sys
from collections import defaultdict, Counter
from statistics import median

BED12_COLS = [
    "chrom", "start", "end", "name", "score", "strand",
    "thickStart", "thickEnd", "rgb", "blockCount", "blockSizes", "blockStarts"
]


def open_maybe_gzip(path, mode="rt"):
    if path is None or path == "":
        return None
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_bed12(path):
    rows = []
    with open_maybe_gzip(path, "rt") as fh:
        for line_no, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 12:
                sys.stderr.write(f"[WARN] skip non-BED12 line {line_no}: {line}\n")
                continue
            rec = dict(zip(BED12_COLS, parts[:12]))
            try:
                for k in ["start", "end", "thickStart", "thickEnd", "blockCount"]:
                    rec[k] = int(rec[k])
            except ValueError:
                sys.stderr.write(f"[WARN] skip line with non-integer coordinates {line_no}: {line}\n")
                continue
            rec["span"] = rec["end"] - rec["start"]
            rec["source_transcript_id"], rec["gene_label"], rec["source_name"] = parse_name(rec["name"])
            rec["exons"] = bed12_exons(rec)
            rec["exonic_bp"] = sum(e2 - e1 for e1, e2 in rec["exons"])
            rows.append(rec)
    return rows


def parse_name(name):
    """
    Examples:
      hg38.ENST00000312403.10#PGA5#8382
      mm10.ENSMUST00000025647#Pga5#213
      HLbosTau10.ENSBTAT00000013786#PAG7#27
      hg38_PGA5  (seed-like IDs; not expected in BED but handled)
    """
    s = str(name)
    parts = s.split("#")
    source_transcript_id = parts[0]
    gene_label = parts[1] if len(parts) >= 2 else "NA"
    source_name = source_transcript_id + "#" + gene_label
    return source_transcript_id, gene_label, source_name


def bed12_exons(rec):
    sizes = [x for x in str(rec["blockSizes"]).rstrip(",").split(",") if x != ""]
    starts = [x for x in str(rec["blockStarts"]).rstrip(",").split(",") if x != ""]
    exons = []
    if len(sizes) != len(starts):
        return exons
    for size, rel_start in zip(sizes, starts):
        try:
            s = rec["start"] + int(rel_start)
            e = s + int(size)
        except ValueError:
            continue
        if e > s:
            exons.append((s, e))
    exons.sort()
    return exons


def read_fasta_info(path, protein=False):
    """
    Return dict: id -> {length, internal_stop}
    FASTA ID is the first token after '>'.
    """
    info = {}
    if not path:
        return info
    if not os.path.exists(path):
        sys.stderr.write(f"[WARN] FASTA does not exist: {path}\n")
        return info

    fh = open_maybe_gzip(path, "rt")
    name = None
    seq_chunks = []

    def flush():
        if name is None:
            return
        seq = "".join(seq_chunks).replace(" ", "").replace("\t", "")
        if protein:
            internal = "*" in seq[:-1] if len(seq) > 0 else False
        else:
            internal = False
        info[name] = {"length": len(seq), "internal_stop": internal}

    with fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        flush()

    return info


def load_simple_list(path):
    vals = set()
    if not path:
        return vals
    with open_maybe_gzip(path, "rt") as fh:
        for line in fh:
            x = line.strip()
            if not x or x.startswith("#"):
                continue
            vals.add(x)
    return vals


def load_whitelists(args):
    tx = set()
    srcname = set()
    gene = set()

    tx |= load_simple_list(args.source_transcript_whitelist)
    srcname |= load_simple_list(args.source_name_whitelist)
    gene |= load_simple_list(args.gene_label_whitelist)

    if args.whitelist_tsv:
        with open_maybe_gzip(args.whitelist_tsv, "rt") as fh:
            header = None
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if header is None:
                    header = {name: i for i, name in enumerate(parts)}
                    required_any = {"source_transcript_id", "source_name", "gene_label"}
                    if len(required_any & set(header)) == 0:
                        raise ValueError(
                            "whitelist-tsv must contain at least one of columns: "
                            "source_transcript_id, source_name, gene_label"
                        )
                    continue

                status = parts[header["target_status"]] if "target_status" in header and header["target_status"] < len(parts) else ""
                if (not args.include_anchor_whitelist) and status.lower() == "anchor":
                    continue

                if "source_transcript_id" in header and header["source_transcript_id"] < len(parts):
                    val = parts[header["source_transcript_id"]].strip()
                    if val:
                        tx.add(val)
                if "source_name" in header and header["source_name"] < len(parts):
                    val = parts[header["source_name"]].strip()
                    if val:
                        srcname.add(val)
                if "gene_label" in header and header["gene_label"] < len(parts):
                    val = parts[header["gene_label"]].strip()
                    if val:
                        gene.add(val)

    return tx, srcname, gene


def comma_lower_set(x):
    return set(i.strip().lower() for i in str(x).split(",") if i.strip())


def classify_target(rec, tx_whitelist, source_name_whitelist, gene_whitelist, anchor_labels, gene_regex=None):
    gene = rec["gene_label"]
    gene_lower = str(gene).lower()
    if gene_lower in anchor_labels:
        return "anchor"

    if rec["source_transcript_id"] in tx_whitelist:
        return "target_by_source_transcript"
    if rec["source_name"] in source_name_whitelist:
        return "target_by_source_name"
    if rec["gene_label"] in gene_whitelist:
        return "target_by_gene_label"
    if gene_regex is not None and gene_regex.search(str(rec["gene_label"])):
        return "target_by_gene_regex"
    return "non_target"


def interval_overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def exon_overlap(exons1, exons2):
    i = j = 0
    ov = 0
    shared = 0
    while i < len(exons1) and j < len(exons2):
        a1, a2 = exons1[i]
        b1, b2 = exons2[j]
        x = interval_overlap(a1, a2, b1, b2)
        if x > 0:
            ov += x
            shared += 1
        if a2 <= b2:
            i += 1
        else:
            j += 1
    return ov, shared


def exon_union(exons):
    if not exons:
        return []
    xs = sorted(exons)
    out = [xs[0]]
    for s, e in xs[1:]:
        ps, pe = out[-1]
        if s <= pe:
            out[-1] = (ps, max(pe, e))
        else:
            out.append((s, e))
    return out


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1


def should_cluster_seed(a, b, args):
    if a["chrom"] != b["chrom"] or a["strand"] != b["strand"]:
        return False
    iov = interval_overlap(a["start"], a["end"], b["start"], b["end"])
    min_span = max(1, min(a["span"], b["span"]))
    if iov / min_span >= args.seed_cluster_overlap:
        return True
    eov, shared = exon_overlap(a["exons"], b["exons"])
    min_exonic = max(1, min(a["exonic_bp"], b["exonic_bp"]))
    if eov / min_exonic >= args.seed_exonic_overlap:
        return True
    if shared >= args.seed_shared_exons:
        return True
    return False


def make_unit(unit_id, seed_records, status="complete_local_structure"):
    chroms = sorted(set(r["chrom"] for r in seed_records))
    strands = sorted(set(r["strand"] for r in seed_records))
    chrom = chroms[0]
    strand = strands[0]
    start = min(r["start"] for r in seed_records)
    end = max(r["end"] for r in seed_records)
    exons = exon_union([e for r in seed_records for e in r["exons"]])
    return {
        "unit_id": unit_id,
        "chrom": chrom,
        "start": start,
        "end": end,
        "span": end - start,
        "strand": strand,
        "exons": exons,
        "status": status,
        "complete_seed": [],
        "partial_seed": [],
        "weak_seed": [],
        "attached_fragment": [],
        "stretched_projection": [],
        "fused_projection": [],
    }


def unit_exonic_bp(unit):
    return sum(e2 - e1 for e1, e2 in unit["exons"])


def find_best_unit_for_record(rec, units, args):
    best = None
    best_score = -1
    best_metrics = None
    for u in units:
        if rec["chrom"] != u["chrom"] or rec["strand"] != u["strand"]:
            continue
        iov = interval_overlap(rec["start"], rec["end"], u["start"], u["end"])
        interval_frac_rec = iov / max(1, rec["span"])
        interval_frac_unit = iov / max(1, u["span"])
        eov, shared = exon_overlap(rec["exons"], u["exons"])
        exon_frac_rec = eov / max(1, rec["exonic_bp"])
        exon_frac_unit = eov / max(1, unit_exonic_bp(u))
        score = max(interval_frac_rec, interval_frac_unit, exon_frac_rec, exon_frac_unit) + 0.05 * shared
        if score > best_score:
            best_score = score
            best = u
            best_metrics = {
                "interval_overlap_bp": iov,
                "interval_frac_rec": interval_frac_rec,
                "interval_frac_unit": interval_frac_unit,
                "exon_overlap_bp": eov,
                "exon_frac_rec": exon_frac_rec,
                "exon_frac_unit": exon_frac_unit,
                "shared_exons": shared,
            }
    return best, best_metrics


def is_attachable(metrics, args):
    if metrics is None:
        return False
    if metrics["interval_overlap_bp"] >= args.attach_min_bp and metrics["interval_frac_rec"] >= args.attach_overlap:
        return True
    if metrics["exon_frac_rec"] >= args.attach_exonic_overlap:
        return True
    if metrics["shared_exons"] >= args.attach_shared_exons:
        return True
    return False


def is_stretched_projection(rec, unit, metrics, args):
    if unit is None or metrics is None:
        return False
    if rec["span"] < max(1, unit["span"]) * args.stretch_fold:
        return False
    if metrics["exon_frac_rec"] >= args.stretch_exon_ov:
        return True
    if metrics["exon_frac_unit"] >= args.stretch_exon_ov_unit:
        return True
    if metrics["shared_exons"] >= args.stretch_shared_exons:
        return True
    return False


def representative_score(rec, median_compact_span, args, nuc_info, prot_info):
    score = 0
    tid = rec["source_transcript_id"]
    gene = rec["gene_label"]
    name = rec["name"]
    status = rec.get("target_status", "")

    if status == "target_by_source_transcript":
        score += 500
    elif status == "target_by_source_name":
        score += 400
    elif status == "target_by_gene_label":
        score += 200
    elif status == "target_by_gene_regex":
        score += 50

    if args.prefer_id_prefix and tid.startswith(args.prefer_id_prefix):
        score += 300
    if rec["blockCount"] == args.target_exons:
        score += 250
    if re.search(args.prefer_gene_regex, str(gene)):
        score += 120

    if args.min_seed_span <= rec["span"] <= args.max_compact_seed_span:
        score += 150
    if median_compact_span:
        score -= min(200, abs(rec["span"] - median_compact_span) / max(1, median_compact_span) * 100)

    if name in nuc_info:
        score += 40
    # FASTA IDs normally match name. Also try source transcript id for reference FASTAs.
    pinfo = prot_info.get(name) or prot_info.get(tid)
    if pinfo:
        score += 60
        plen = pinfo["length"]
        if args.min_protein_len <= plen <= args.max_protein_len:
            score += 80
        if not pinfo["internal_stop"]:
            score += 60
        else:
            score -= 80

    if rec["span"] >= args.fused_hard_span:
        score -= 1000
    return score


def choose_representative(records, median_compact_span, args, nuc_info, prot_info):
    if not records:
        return None
    return sorted(records, key=lambda r: representative_score(r, median_compact_span, args, nuc_info, prot_info), reverse=True)[0]


def format_list(records, field="name"):
    if not records:
        return "NA"
    vals = []
    seen = set()
    for r in records:
        v = str(r[field])
        if v not in seen:
            vals.append(v)
            seen.add(v)
    return ";".join(vals) if vals else "NA"


def write_lines(path, vals):
    with open(path, "w") as out:
        for v in vals:
            if v and v != "NA":
                out.write(str(v) + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Refine TOGA PGA/PGA-like local units with canonical-source whitelist filtering."
    )
    ap.add_argument("--bed", required=True, help="TOGA BED12 annotations inside VPS37C--VWCE query interval")
    ap.add_argument("--assembly", default="NA", help="Assembly/sample name")
    ap.add_argument("--out-prefix", required=True, help="Output prefix")
    ap.add_argument("--nucleotide", default=None, help="Optional nucleotide FASTA(.gz) for representative QC")
    ap.add_argument("--protein", default=None, help="Optional protein FASTA(.gz) for representative QC")

    # Whitelists
    ap.add_argument("--whitelist-tsv", default=None, help="Canonical interval whitelist TSV with source_transcript_id/source_name/gene_label columns")
    ap.add_argument("--source-transcript-whitelist", default=None, help="One source transcript ID per line")
    ap.add_argument("--source-name-whitelist", default=None, help="One source_name transcript#gene per line")
    ap.add_argument("--gene-label-whitelist", default=None, help="One gene label per line")
    ap.add_argument("--include-anchor-whitelist", action="store_true", help="Do not drop target_status=anchor rows from whitelist TSV")
    ap.add_argument("--use-gene-regex-fallback", action="store_true", help="Allow gene-label regex fallback as target evidence")
    ap.add_argument("--target-gene-regex", default=r"^(PGA[0-9A-Za-z_.-]*|Pga[0-9A-Za-z_.-]*|PAG[0-9A-Za-z_.-]*|Pag[0-9A-Za-z_.-]*)$", help="Fallback target gene regex")
    ap.add_argument("--anchor-labels", default="VPS37C,VWCE,Vps37c,Vwce", help="Comma-separated anchor gene labels to exclude from target units")

    # Local-unit parameters
    ap.add_argument("--target-exons", type=int, default=9)
    ap.add_argument("--min-seed-span", type=int, default=4000)
    ap.add_argument("--max-compact-seed-span", type=int, default=15000)
    ap.add_argument("--max-seed-span", type=int, default=18000)
    ap.add_argument("--partial-min-exons", type=int, default=5)
    ap.add_argument("--weak-min-exons", type=int, default=2)

    ap.add_argument("--seed-cluster-overlap", type=float, default=0.35)
    ap.add_argument("--seed-exonic-overlap", type=float, default=0.50)
    ap.add_argument("--seed-shared-exons", type=int, default=3)

    ap.add_argument("--attach-overlap", type=float, default=0.35)
    ap.add_argument("--attach-exonic-overlap", type=float, default=0.35)
    ap.add_argument("--attach-shared-exons", type=int, default=2)
    ap.add_argument("--attach-min-bp", type=int, default=500)

    ap.add_argument("--fused-hard-span", type=int, default=35000)
    ap.add_argument("--fused-fold", type=float, default=2.0)
    ap.add_argument("--min-units-for-bridging", type=int, default=2)

    ap.add_argument("--stretch-fold", type=float, default=1.35)
    ap.add_argument("--stretch-exon-ov", type=float, default=0.30)
    ap.add_argument("--stretch-exon-ov-unit", type=float, default=0.50)
    ap.add_argument("--stretch-shared-exons", type=int, default=3)

    # Representative preferences
    ap.add_argument("--prefer-id-prefix", default="hg38.ENST")
    ap.add_argument("--prefer-gene-regex", default=r"^(PGA|Pga|PAG|Pag)")
    ap.add_argument("--min-protein-len", type=int, default=250)
    ap.add_argument("--max-protein-len", type=int, default=500)

    # Rearrangement QC
    ap.add_argument("--rearrange-min-non-target-genes", type=int, default=3)
    ap.add_argument("--rearrange-min-non-target-frac", type=float, default=0.30)

    args = ap.parse_args()

    rows = parse_bed12(args.bed)
    if len(rows) == 0:
        raise ValueError("No valid BED12 rows were read.")

    nuc_info = read_fasta_info(args.nucleotide, protein=False) if args.nucleotide else {}
    prot_info = read_fasta_info(args.protein, protein=True) if args.protein else {}

    tx_w, source_name_w, gene_w = load_whitelists(args)
    anchor_labels = comma_lower_set(args.anchor_labels)
    gene_regex = re.compile(args.target_gene_regex) if args.use_gene_regex_fallback else None

    for r in rows:
        r["target_status"] = classify_target(r, tx_w, source_name_w, gene_w, anchor_labels, gene_regex)

    target_rows = [r for r in rows if r["target_status"].startswith("target_by")]
    anchor_rows = [r for r in rows if r["target_status"] == "anchor"]
    non_target_rows = [r for r in rows if r["target_status"] == "non_target"]

    compact_seeds = [
        r for r in target_rows
        if r["blockCount"] == args.target_exons
        and args.min_seed_span <= r["span"] <= args.max_compact_seed_span
    ]

    median_compact_span = median([r["span"] for r in compact_seeds]) if compact_seeds else None

    # Build primary units from compact complete seeds.
    units = []
    used_ids = set()
    if compact_seeds:
        uf = UnionFind(len(compact_seeds))
        for i in range(len(compact_seeds)):
            for j in range(i + 1, len(compact_seeds)):
                if should_cluster_seed(compact_seeds[i], compact_seeds[j], args):
                    uf.union(i, j)
        groups = defaultdict(list)
        for i, r in enumerate(compact_seeds):
            groups[uf.find(i)].append(r)
        for idx, group in enumerate(sorted(groups.values(), key=lambda g: (g[0]["chrom"], min(x["start"] for x in g))), 1):
            u = make_unit(f"unit{idx:04d}", group, status="complete_local_structure")
            u["complete_seed"].extend(group)
            units.append(u)
            for r in group:
                used_ids.add(id(r))

    event_rows = []

    def add_event(rec, event_type, unit_ids, metrics=None):
        event_rows.append({
            "assembly": args.assembly,
            "chrom": rec["chrom"],
            "start": rec["start"],
            "end": rec["end"],
            "span": rec["span"],
            "strand": rec["strand"],
            "name": rec["name"],
            "source_transcript_id": rec["source_transcript_id"],
            "gene_label": rec["gene_label"],
            "target_status": rec["target_status"],
            "event_type": event_type,
            "overlapped_units": ";".join(unit_ids) if unit_ids else "NA",
            "interval_overlap_bp": metrics.get("interval_overlap_bp", "NA") if metrics else "NA",
            "exon_overlap_bp": metrics.get("exon_overlap_bp", "NA") if metrics else "NA",
            "shared_exons": metrics.get("shared_exons", "NA") if metrics else "NA",
        })

    # Process remaining target rows in coordinate order.
    remaining = [r for r in target_rows if id(r) not in used_ids]
    remaining.sort(key=lambda r: (r["chrom"], r["start"], r["end"], r["name"]))

    unit_counter = len(units)

    for rec in remaining:
        # Determine best overlap to existing unit.
        best_unit, metrics = find_best_unit_for_record(rec, units, args) if units else (None, None)

        # Large fused/bridging target projections should not define copy boundaries.
        is_hard_fused = rec["span"] >= args.fused_hard_span
        is_relative_fused = False
        if median_compact_span is not None and rec["span"] >= median_compact_span * args.fused_fold and rec["span"] > args.max_seed_span:
            is_relative_fused = True

        overlapped_units = []
        for u in units:
            if rec["chrom"] != u["chrom"] or rec["strand"] != u["strand"]:
                continue
            iov = interval_overlap(rec["start"], rec["end"], u["start"], u["end"])
            eov, shared = exon_overlap(rec["exons"], u["exons"])
            if iov > 0 or eov > 0 or shared > 0:
                overlapped_units.append(u)

        is_bridging = len(overlapped_units) >= args.min_units_for_bridging and (is_hard_fused or is_relative_fused or rec["span"] > args.max_seed_span)

        if is_hard_fused or is_relative_fused or is_bridging:
            unit_ids = []
            for u in overlapped_units:
                u["fused_projection"].append(rec)
                unit_ids.append(u["unit_id"])
            add_event(rec, "large_fused_or_bridging_projection", unit_ids, metrics or {})
            continue

        # Stretched long-intron target projection; attach to best compact/partial unit.
        if best_unit is not None and is_stretched_projection(rec, best_unit, metrics, args):
            best_unit["stretched_projection"].append(rec)
            add_event(rec, "stretched_projection", [best_unit["unit_id"]], metrics)
            continue

        # Attach fragments/alternative projections to existing units.
        if best_unit is not None and is_attachable(metrics, args):
            best_unit["attached_fragment"].append(rec)
            add_event(rec, "attached_target_fragment", [best_unit["unit_id"]], metrics)
            continue

        # If not attached, allow target models to seed new units.
        if rec["blockCount"] == args.target_exons and args.min_seed_span <= rec["span"] <= args.max_seed_span:
            unit_counter += 1
            u = make_unit(f"unit{unit_counter:04d}", [rec], status="complete_local_structure")
            u["complete_seed"].append(rec)
            units.append(u)
        elif rec["blockCount"] >= args.partial_min_exons and args.min_seed_span <= rec["span"] <= args.max_seed_span:
            unit_counter += 1
            u = make_unit(f"unit{unit_counter:04d}", [rec], status="partial_local_structure")
            u["partial_seed"].append(rec)
            units.append(u)
        elif rec["blockCount"] >= args.weak_min_exons:
            unit_counter += 1
            u = make_unit(f"unit{unit_counter:04d}", [rec], status="weak_fragment_structure")
            u["weak_seed"].append(rec)
            units.append(u)
        else:
            add_event(rec, "tiny_target_fragment_unassigned", [], metrics or {})

    # Re-sort units and rename by coordinate for stable output.
    units.sort(key=lambda u: (u["chrom"], u["start"], u["end"], u["unit_id"]))
    old_to_new = {}
    for i, u in enumerate(units, 1):
        old = u["unit_id"]
        new = f"unit{i:04d}"
        old_to_new[old] = new
        u["unit_id"] = new

    # Prepare outputs.
    local_unit_rows = []
    all_reps = []
    complete_reps = []
    assign_reps = []

    for u in units:
        core_records = u["complete_seed"] + u["partial_seed"] + u["weak_seed"]
        rep = choose_representative(core_records, median_compact_span, args, nuc_info, prot_info)
        if rep is None:
            continue

        unit_status = u["status"]
        if len(u["complete_seed"]) > 0:
            unit_status = "complete_local_structure"
        elif len(u["partial_seed"]) > 0:
            unit_status = "partial_local_structure"
        else:
            unit_status = "weak_fragment_structure"

        pinfo = prot_info.get(rep["name"]) or prot_info.get(rep["source_transcript_id"]) or {}
        ninfo = nuc_info.get(rep["name"]) or nuc_info.get(rep["source_transcript_id"]) or {}

        all_reps.append(rep["name"])
        if unit_status == "complete_local_structure":
            complete_reps.append(rep["name"])
            assign_reps.append(rep["name"])
        elif unit_status == "partial_local_structure":
            assign_reps.append(rep["name"])

        local_unit_rows.append({
            "assembly": args.assembly,
            "unit_id": u["unit_id"],
            "chrom": u["chrom"],
            "start": u["start"],
            "end": u["end"],
            "span": u["span"],
            "strand": u["strand"],
            "unit_status": unit_status,
            "representative_id": rep["name"],
            "representative_source_transcript_id": rep["source_transcript_id"],
            "representative_gene_label": rep["gene_label"],
            "representative_target_status": rep["target_status"],
            "representative_blockCount": rep["blockCount"],
            "representative_span": rep["span"],
            "representative_exonic_bp": rep["exonic_bp"],
            "nucleotide_length": ninfo.get("length", "NA"),
            "protein_length": pinfo.get("length", "NA"),
            "protein_internal_stop": pinfo.get("internal_stop", "NA"),
            "n_complete_seed_models": len(u["complete_seed"]),
            "n_partial_seed_models": len(u["partial_seed"]),
            "n_weak_seed_models": len(u["weak_seed"]),
            "n_attached_fragment_models": len(u["attached_fragment"]),
            "n_stretched_projection_models": len(u["stretched_projection"]),
            "n_fused_models_overlapping": len(u["fused_projection"]),
            "complete_seed_ids": format_list(u["complete_seed"]),
            "partial_seed_ids": format_list(u["partial_seed"]),
            "weak_seed_ids": format_list(u["weak_seed"]),
            "attached_fragment_ids": format_list(u["attached_fragment"]),
            "stretched_projection_ids": format_list(u["stretched_projection"]),
            "overlapping_fused_ids": format_list(u["fused_projection"]),
        })

    # Write table helper.
    def write_tsv(path, records, fields):
        with open(path, "w") as out:
            out.write("\t".join(fields) + "\n")
            for r in records:
                out.write("\t".join(str(r.get(f, "NA")) for f in fields) + "\n")

    unit_fields = [
        "assembly", "unit_id", "chrom", "start", "end", "span", "strand", "unit_status",
        "representative_id", "representative_source_transcript_id", "representative_gene_label",
        "representative_target_status", "representative_blockCount", "representative_span",
        "representative_exonic_bp", "nucleotide_length", "protein_length", "protein_internal_stop",
        "n_complete_seed_models", "n_partial_seed_models", "n_weak_seed_models",
        "n_attached_fragment_models", "n_stretched_projection_models", "n_fused_models_overlapping",
        "complete_seed_ids", "partial_seed_ids", "weak_seed_ids", "attached_fragment_ids",
        "stretched_projection_ids", "overlapping_fused_ids",
    ]
    write_tsv(args.out_prefix + ".local_units.tsv", local_unit_rows, unit_fields)

    event_fields = [
        "assembly", "chrom", "start", "end", "span", "strand", "name",
        "source_transcript_id", "gene_label", "target_status", "event_type", "overlapped_units",
        "interval_overlap_bp", "exon_overlap_bp", "shared_exons",
    ]
    write_tsv(args.out_prefix + ".fused_or_stretched.tsv", event_rows, event_fields)

    # All model target status table.
    model_rows = []
    for r in rows:
        model_rows.append({
            "assembly": args.assembly,
            "chrom": r["chrom"],
            "start": r["start"],
            "end": r["end"],
            "span": r["span"],
            "strand": r["strand"],
            "name": r["name"],
            "source_transcript_id": r["source_transcript_id"],
            "gene_label": r["gene_label"],
            "source_name": r["source_name"],
            "blockCount": r["blockCount"],
            "exonic_bp": r["exonic_bp"],
            "target_status": r["target_status"],
        })
    model_fields = [
        "assembly", "chrom", "start", "end", "span", "strand", "name", "source_transcript_id",
        "gene_label", "source_name", "blockCount", "exonic_bp", "target_status",
    ]
    write_tsv(args.out_prefix + ".all_model_target_status.tsv", model_rows, model_fields)

    # Non-target models table for rearrangement inspection.
    non_target_records = [r for r in model_rows if r["target_status"] == "non_target"]
    write_tsv(args.out_prefix + ".non_target_models.tsv", non_target_records, model_fields)

    write_lines(args.out_prefix + ".all_unit_representative_ids.txt", all_reps)
    write_lines(args.out_prefix + ".complete_unit_representative_ids.txt", complete_reps)
    write_lines(args.out_prefix + ".assign_candidate_ids.txt", assign_reps)

    # Summary and rearrangement QC.
    status_count = Counter(r["unit_status"] for r in local_unit_rows)
    target_status_count = Counter(r["target_status"] for r in model_rows)
    non_target_gene_labels = sorted(set(r["gene_label"] for r in non_target_records if str(r["gene_label"]).lower() not in anchor_labels))
    n_total = len(rows)
    n_non_target = len(non_target_records)
    non_target_frac = n_non_target / max(1, n_total)
    rearr_flag = (len(non_target_gene_labels) >= args.rearrange_min_non_target_genes) or (non_target_frac >= args.rearrange_min_non_target_frac)
    reason = []
    if len(non_target_gene_labels) >= args.rearrange_min_non_target_genes:
        reason.append(f"n_non_target_gene_labels>={args.rearrange_min_non_target_genes}")
    if non_target_frac >= args.rearrange_min_non_target_frac:
        reason.append(f"non_target_frac>={args.rearrange_min_non_target_frac}")
    if not reason:
        reason.append("NA")

    summary_rows = []
    summary_rows.append(("assembly", args.assembly))
    summary_rows.append(("total_bed_models", n_total))
    summary_rows.append(("target_models", len(target_rows)))
    summary_rows.append(("anchor_models", len(anchor_rows)))
    summary_rows.append(("non_target_models", n_non_target))
    summary_rows.append(("non_target_fraction", f"{non_target_frac:.4f}"))
    summary_rows.append(("non_target_gene_labels", ";".join(non_target_gene_labels) if non_target_gene_labels else "NA"))
    summary_rows.append(("local_rearrangement_flag", str(bool(rearr_flag))))
    summary_rows.append(("local_rearrangement_reason", ";".join(reason)))
    summary_rows.append(("total_local_units", len(local_unit_rows)))
    for k in sorted(status_count):
        summary_rows.append((k, status_count[k]))
    for k in sorted(target_status_count):
        summary_rows.append(("model_" + k, target_status_count[k]))
    summary_rows.append(("large_fused_like_rows", sum(1 for e in event_rows if e["event_type"] == "large_fused_or_bridging_projection")))
    summary_rows.append(("stretched_projection_rows", sum(1 for e in event_rows if e["event_type"] == "stretched_projection")))
    summary_rows.append(("median_compact_seed_span", median_compact_span if median_compact_span is not None else "NA"))
    summary_rows.append(("whitelist_source_transcript_n", len(tx_w)))
    summary_rows.append(("whitelist_source_name_n", len(source_name_w)))
    summary_rows.append(("whitelist_gene_label_n", len(gene_w)))

    with open(args.out_prefix + ".local_units.summary.tsv", "w") as out:
        out.write("metric\tvalue\n")
        for k, v in summary_rows:
            out.write(f"{k}\t{v}\n")

    rearr_fields = [
        "assembly", "total_bed_models", "target_models", "anchor_models", "non_target_models",
        "non_target_fraction", "n_non_target_gene_labels", "non_target_gene_labels",
        "local_rearrangement_flag", "local_rearrangement_reason",
        "total_local_units", "complete_local_structure", "partial_local_structure", "weak_fragment_structure",
    ]
    rearr_record = {
        "assembly": args.assembly,
        "total_bed_models": n_total,
        "target_models": len(target_rows),
        "anchor_models": len(anchor_rows),
        "non_target_models": n_non_target,
        "non_target_fraction": f"{non_target_frac:.4f}",
        "n_non_target_gene_labels": len(non_target_gene_labels),
        "non_target_gene_labels": ";".join(non_target_gene_labels) if non_target_gene_labels else "NA",
        "local_rearrangement_flag": str(bool(rearr_flag)),
        "local_rearrangement_reason": ";".join(reason),
        "total_local_units": len(local_unit_rows),
        "complete_local_structure": status_count.get("complete_local_structure", 0),
        "partial_local_structure": status_count.get("partial_local_structure", 0),
        "weak_fragment_structure": status_count.get("weak_fragment_structure", 0),
    }
    write_tsv(args.out_prefix + ".rearrangement_qc.tsv", [rearr_record], rearr_fields)

    print(f"[INFO] assembly: {args.assembly}")
    print(f"[INFO] total BED models: {n_total}")
    print(f"[INFO] target models: {len(target_rows)}")
    print(f"[INFO] non-target models: {n_non_target}")
    print(f"[INFO] local rearrangement flag: {bool(rearr_flag)}")
    print(f"[INFO] total local units: {len(local_unit_rows)}")
    print(f"[INFO] complete local units: {status_count.get('complete_local_structure', 0)}")
    print(f"[INFO] partial local units: {status_count.get('partial_local_structure', 0)}")
    print(f"[INFO] weak fragment units: {status_count.get('weak_fragment_structure', 0)}")
    print(f"[INFO] output prefix: {args.out_prefix}")


if __name__ == "__main__":
    main()
