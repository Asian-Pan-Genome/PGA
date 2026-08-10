#!/usr/bin/env python3
"""
Draw compact Fig.6b case panels for primate and sirenian PGA duplicons.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

import call_untangle_duplicon_TE as call
import plot_untangle_duplicon_TE as base


PRIMATE_CASE_SPECIES = [
    "Ateles_hybridus",
    "Pan_troglodytes",
    "Symphalangus_syndactylus",
    "Rhinopithecus_roxellana",
    "Macaca_mulatta",
]
ATELES_REFERENCE_INTERVAL = (41761, 63116)
ATELES_DIRNAME = "Ateles_hybridus__brown_spider_monkey__HLateHyb1__GCA_916098195.1"

TRICHECHUS_DIRNAME = "Trichechus_inunguis__Amazon_manatee__HLtriInu1A__GCA_046562895.1"
TRICHECHUS_LABEL = "Trichechus_inunguis"


def parse_args() -> argparse.Namespace:
    base_dir = Path(__file__).resolve().parent
    input_dir = base_dir / "untangle_duplicon_TE"
    parser = argparse.ArgumentParser(description="Draw compact Fig.6b duplicon case panels.")
    parser.add_argument("--base-dir", type=Path, default=base_dir)
    parser.add_argument("--input-dir", type=Path, default=input_dir)
    parser.add_argument("--tree", type=Path, default=base_dir.parent / "295.sp.tree")
    parser.add_argument("--reference", default="Ateles_hybridus")
    parser.add_argument("--m", default="256")
    parser.add_argument("--output-dir", type=Path, default=input_dir)
    parser.add_argument("--max-join-transition-len", type=int, default=5000)
    parser.add_argument("--flank-kb", type=float, default=10)
    parser.add_argument("--track-height-ratio", default="2:1:1")
    parser.add_argument("--primates-width", type=float, default=4.8)
    parser.add_argument("--primates-height", type=float, default=2.15)
    parser.add_argument("--sirenia-width", type=float, default=4.8)
    parser.add_argument("--sirenia-height", type=float, default=0.9)
    return parser.parse_args()


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def read_fasta_len(path: Path) -> int:
    total = 0
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def cluster_points(points: list[int], tolerance: int = 500) -> list[int]:
    if not points:
        return []
    points = sorted(points)
    clusters: list[list[int]] = [[points[0]]]
    for point in points[1:]:
        if point - clusters[-1][-1] <= tolerance:
            clusters[-1].append(point)
        else:
            clusters.append([point])
    return [min(cluster) for cluster in clusters]


def ateles_reference_interval_rows(base_dir: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    start, end = ATELES_REFERENCE_INTERVAL
    intervals = [
        {
            "duplicon_id": "Ateles_hybridus|reference|41761-63116",
            "query_species": "Ateles_hybridus",
            "gene": "PGA_1",
            "pga_strand": "+",
            "bio_order": "1",
            "effective_duplicon_raw_start": str(start),
            "effective_duplicon_raw_end": str(end),
            "effective_duplicon_len": str(end - start),
            "join_status": "reference_interval",
            "endpoint_source": "Ateles_reference_interval",
        }
    ]
    rm_path = base_dir / ATELES_DIRNAME / f"{ATELES_DIRNAME}.pga.anchor.locus.fa.out"
    endpoints: list[dict[str, str]] = []
    if rm_path.exists():
        repeats = call.parse_repeatmasker(rm_path)
        for role, endpoint, boundary_kind in (
            ("biological_start", start, "raw_start"),
            ("biological_end", end, "raw_end"),
        ):
            hits = call.endpoint_hits(repeats, endpoint, boundary_kind, 0)
            if not hits:
                continue
            for hit in hits:
                endpoints.append(
                    {
                        "query_species": "Ateles_hybridus",
                        "gene": "PGA_1",
                        "pga_strand": "+",
                        "duplicon_id": "Ateles_hybridus|reference|41761-63116",
                        "endpoint_role": role,
                        "endpoint_source": "Ateles_reference_interval",
                        "raw_boundary": boundary_kind,
                        "endpoint_raw": str(endpoint),
                        "endpoint_1based": str(hit["endpoint_1based"]),
                        "te_hit": "1",
                        "te_name": hit["repeat_name"],
                        "te_class": hit["repeat_class"],
                        "te_group": hit["repeat_group"],
                        "te_start_1based": str(hit["start_1based"]),
                        "te_end_1based": str(hit["end_1based"]),
                        "te_start0": str(hit["start0"]),
                        "te_end0": str(hit["end0"]),
                        "rm_id": hit["rm_id"],
                        "dist_to_endpoint": str(hit["dist_to_endpoint"]),
                    }
                )
    return intervals, endpoints


def plot_primate_tracks_compact(
    out_pdf: Path,
    out_png: Path,
    tree,
    species_order: list[str],
    genes_by_species: dict[str, list[dict[str, str]]],
    length_by_species: dict[str, int],
    pga_strand_by_species: dict[str, str],
    crop_by_species: dict[str, tuple[float, float]],
    track_heights: tuple[float, float, float],
    intervals: list[dict[str, str]],
    endpoint_rows: list[dict[str, str]],
    plot_width: float,
    plot_height: float,
) -> None:
    n = len(species_order)
    y_by_species = {species: n - 1 - idx for idx, species in enumerate(species_order)}
    gene_height, duplicon_height, te_height = track_heights
    gene_y_offset = 0.18
    duplicon_y_offset = -0.10
    fig = plt.figure(figsize=(plot_width, plot_height))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.72, 1.95, 5.7], wspace=0.0)
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_label = fig.add_subplot(gs[0, 1], sharey=ax_tree)
    ax = fig.add_subplot(gs[0, 2], sharey=ax_tree)

    base.draw_tree(ax_tree, tree, species_order, y_by_species)
    ax_label.set_xlim(0, 1)
    ax_label.set_ylim(-0.7, n - 0.3)
    ax_label.axis("off")
    for species in species_order:
        ax_label.text(0.03, y_by_species[species], species, va="center", ha="left", fontsize=5.2, clip_on=False)

    max_x = max(
        (crop_by_species.get(sp, (0.0, 1.0))[1] - crop_by_species.get(sp, (0.0, 1.0))[0] for sp in species_order),
        default=1.0,
    )
    for species in species_order:
        y = y_by_species[species]
        strand = pga_strand_by_species.get(species, "+")
        seq_len = length_by_species.get(species, 0)
        crop_start = crop_by_species.get(species, (0.0, max_x))[0]
        for gene in genes_by_species.get(species, []):
            start, end = base.interval_plot(base.as_int(gene["local_start"]), base.as_int(gene["local_end"]), strand, seq_len, crop_start)
            base.draw_gene_arrow(ax, start, end, y + gene_y_offset, "+", "#2b8f7b", gene_height)

    for row in intervals:
        species = row["query_species"]
        if species not in y_by_species or species not in crop_by_species:
            continue
        seq_len = length_by_species.get(species, 0)
        strand = row.get("pga_strand", "+")
        start, end = base.interval_plot(
            base.as_int(row["effective_duplicon_raw_start"]),
            base.as_int(row["effective_duplicon_raw_end"]),
            strand,
            seq_len,
            crop_by_species[species][0],
        )
        color = "#6f6f6f" if base.as_int(row.get("bio_order", 0)) % 2 else "#a8a8a8"
        base.draw_interval(ax, start, end, y_by_species[species] + duplicon_y_offset, duplicon_height, color, 0.38)

    base.draw_endpoint_ticks(
        ax,
        endpoint_rows,
        y_by_species,
        length_by_species,
        crop_by_species,
        duplicon_y_offset,
        te_height,
    )

    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(0, max_x * 1.02)
    ax.set_xlabel("Relative PGA locus coordinate (kb)", fontsize=8)
    ax.set_yticks([y_by_species[sp] for sp in species_order])
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=6)
    tick_step = 50 if max_x > 120 else 25
    tick_max = int(max_x // tick_step + 1) * tick_step
    ax.set_xticks(list(range(0, tick_max + 1, tick_step)))
    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_primates_case(args: argparse.Namespace, track_heights: tuple[float, float, float]) -> None:
    selected = set(PRIMATE_CASE_SPECIES)
    gene_rows = base.read_tsv(args.input_dir / "untangle_PGA_genes.local.tsv")
    duplicon_rows = base.read_tsv(args.input_dir / "untangle_duplicon_intervals.tsv")
    core_endpoint_rows = base.read_tsv(args.input_dir / "untangle_core_endpoint_TE.tsv")
    duplicon_endpoint_rows = base.read_tsv(args.input_dir / "untangle_duplicon_endpoint_TE.tsv")
    effective_dups, effective_endpoint = base.build_effective_duplicons(
        duplicon_rows,
        core_endpoint_rows,
        duplicon_endpoint_rows,
        args.reference,
        args.m,
        args.max_join_transition_len,
    )
    tree, tree_order = base.load_pruned_tree(args.tree, PRIMATE_CASE_SPECIES)
    species_order = [species for species in tree_order if species in selected]
    for species in PRIMATE_CASE_SPECIES:
        if species not in species_order:
            species_order.append(species)

    genes_by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    pga_strand_by_species: dict[str, str] = {}
    length_by_species: dict[str, int] = {}
    for row in gene_rows:
        species = row["species"]
        if species not in selected:
            continue
        genes_by_species[species].append(row)
        pga_strand_by_species[species] = row["pga_strand"]
        length_by_species[species] = max(
            length_by_species.get(species, 0),
            base.as_int(row["anchor_end"]) - base.as_int(row["anchor_start"]),
        )
    for species in genes_by_species:
        genes_by_species[species].sort(key=lambda row: base.as_float(row["oriented_start"]))

    intervals = [
        row
        for row in effective_dups
        if row["reference_name"] == args.reference and row["m"] == args.m and row["query_species"] in selected
    ]
    endpoints = [
        row
        for row in effective_endpoint
        if row["reference_name"] == args.reference and row["m"] == args.m and row["query_species"] in selected
    ]
    ateles_intervals, ateles_endpoints = ateles_reference_interval_rows(args.base_dir)
    intervals.extend(ateles_intervals)
    endpoints.extend(ateles_endpoints)
    crop_by_species = base.build_crop_windows(
        species_order,
        intervals,
        genes_by_species,
        length_by_species,
        pga_strand_by_species,
        args.flank_kb,
    )

    out_prefix = args.output_dir / "fig6b_case_primates"
    plot_primate_tracks_compact(
        out_prefix.with_suffix(".pdf"),
        out_prefix.with_suffix(".png"),
        tree,
        species_order,
        genes_by_species,
        length_by_species,
        pga_strand_by_species,
        crop_by_species,
        track_heights,
        intervals,
        endpoints,
        args.primates_width,
        args.primates_height,
    )
    write_tsv(
        out_prefix.with_suffix(".species_order.tsv"),
        [{"plot_order": idx + 1, "species": species} for idx, species in enumerate(species_order)],
        ["plot_order", "species"],
    )
    write_tsv(
        out_prefix.with_suffix(".Ateles_reference_interval.tsv"),
        ateles_intervals,
        [
            "duplicon_id",
            "query_species",
            "gene",
            "pga_strand",
            "bio_order",
            "effective_duplicon_raw_start",
            "effective_duplicon_raw_end",
            "effective_duplicon_len",
            "join_status",
            "endpoint_source",
        ],
    )
    write_tsv(
        out_prefix.with_suffix(".Ateles_reference_interval_endpoint_TE.tsv"),
        ateles_endpoints,
        [
            "query_species",
            "gene",
            "pga_strand",
            "duplicon_id",
            "endpoint_role",
            "endpoint_source",
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
        ],
    )


def build_trichechus_intervals(
    trichechus_dir: Path,
    min_effective_span: int = 20000,
    diagonal_buffer: int = 1000,
) -> tuple[list[dict], list[int]]:
    blocks = base.read_tsv(
        trichechus_dir
        / "Trichechus_inunguis__Amazon_manatee__HLtriInu1A__GCA_046562895.1.pga.self_aln_merged_blocks.tsv"
    )
    non_diag = []
    points: list[int] = []
    for row in blocks:
        ref_start = base.as_int(row["refStart"])
        query_start = base.as_int(row["queryStart"])
        if abs(ref_start - query_start) < diagonal_buffer:
            continue
        if base.as_int(row["effective_span"]) < min_effective_span:
            continue
        non_diag.append(row)
        for key in ("refStart", "refEnd", "queryStart", "queryEnd"):
            points.append(base.as_int(row[key]))
    boundaries = cluster_points(points, tolerance=500)
    if len(boundaries) < 4:
        raise SystemExit("Could not recover at least four Trichechus self-aln duplicon boundaries.")

    intervals: list[dict] = []
    for idx, (start, end) in enumerate(zip(boundaries[:-1], boundaries[1:]), start=1):
        intervals.append(
            {
                "duplicon_id": f"Trichechus_inunguis|self_aln|duplicon_{idx}",
                "query_species": TRICHECHUS_LABEL,
                "gene": f"duplicon_{idx}",
                "pga_strand": "-",
                "bio_order": idx,
                "effective_duplicon_raw_start": start,
                "effective_duplicon_raw_end": end,
                "effective_duplicon_len": end - start,
                "support_aln_ids": ",".join(
                    sorted(
                        {
                            row["aln_id"]
                            for row in non_diag
                            if start in {base.as_int(row["refStart"]), base.as_int(row["queryStart"])}
                            or end in {base.as_int(row["refEnd"]), base.as_int(row["queryEnd"])}
                            or abs(start - base.as_int(row["refStart"])) <= 5
                            or abs(start - base.as_int(row["queryStart"])) <= 5
                            or abs(end - base.as_int(row["refEnd"])) <= 5
                            or abs(end - base.as_int(row["queryEnd"])) <= 5
                        },
                        key=lambda value: int(value),
                    )
                ),
            }
        )
    return intervals, boundaries


def build_trichechus_endpoint_rows(trichechus_dir: Path, boundaries: list[int]) -> list[dict]:
    rows = base.read_tsv(
        trichechus_dir
        / "Trichechus_inunguis__Amazon_manatee__HLtriInu1A__GCA_046562895.1.pga.self_aln_breakpoints_TE.tsv"
    )
    out: dict[tuple[int, str, str], dict] = {}
    for row in rows:
        if row.get("is_diagonal") != "FALSE":
            continue
        bp = base.as_int(row["bp_pos"])
        nearest = min(boundaries, key=lambda point: abs(point - bp))
        if abs(nearest - bp) > 5:
            continue
        key = (nearest, row.get("TE_name", ""), row.get("TE_group", ""))
        out[key] = {
            "query_species": TRICHECHUS_LABEL,
            "endpoint_raw": nearest,
            "te_hit": "1",
            "te_name": row.get("TE_name", ""),
            "te_group": row.get("TE_group", ""),
            "te_class": row.get("TE_class", ""),
            "dist_to_endpoint": abs(nearest - bp),
            "source_bp_pos": bp,
        }
    return list(out.values())


def write_trichechus_consistency(trichechus_dir: Path, boundaries: list[int], out_path: Path) -> None:
    pairwise_rows = base.read_tsv(trichechus_dir / "Dugong_dugon.to.Trichechus_inunguis_breakpoints_TE.tsv")
    pairwise_target_bps = [
        base.as_int(row["bp_pos"])
        for row in pairwise_rows
        if row.get("side") == "target" and row.get("species_label", "").startswith("Trichechus_inunguis")
    ]
    rows = []
    for idx, boundary in enumerate(boundaries[1:-1], start=1):
        nearest = min(pairwise_target_bps, key=lambda point: abs(point - boundary)) if pairwise_target_bps else ""
        delta = abs(int(nearest) - boundary) if nearest != "" else ""
        rows.append(
            {
                "boundary_id": f"internal_boundary_{idx}",
                "self_aln_boundary_raw": boundary,
                "nearest_pairwise_target_boundary_raw": nearest,
                "delta_bp": delta,
                "consistent_within_5bp": int(delta != "" and int(delta) <= 5),
            }
        )
    write_tsv(
        out_path,
        rows,
        [
            "boundary_id",
            "self_aln_boundary_raw",
            "nearest_pairwise_target_boundary_raw",
            "delta_bp",
            "consistent_within_5bp",
        ],
    )


def plot_sirenia_case(args: argparse.Namespace, track_heights: tuple[float, float, float]) -> None:
    trichechus_dir = args.base_dir / TRICHECHUS_DIRNAME
    fasta = trichechus_dir / f"{TRICHECHUS_DIRNAME}.pga.anchor.locus.fa"
    seq_len = read_fasta_len(fasta)
    genes = base.read_tsv(trichechus_dir / f"{TRICHECHUS_DIRNAME}.pga.self_aln_genes.local.tsv")
    intervals, boundaries = build_trichechus_intervals(trichechus_dir)
    endpoints = build_trichechus_endpoint_rows(trichechus_dir, boundaries)
    out_prefix = args.output_dir / "fig6b_case_sirenia_Trichechus_inunguis"
    write_trichechus_consistency(trichechus_dir, boundaries, out_prefix.with_suffix(".self_pairwise_boundary_consistency.tsv"))
    write_tsv(
        out_prefix.with_suffix(".self_aln_duplicon_intervals.tsv"),
        intervals,
        [
            "duplicon_id",
            "query_species",
            "gene",
            "pga_strand",
            "bio_order",
            "effective_duplicon_raw_start",
            "effective_duplicon_raw_end",
            "effective_duplicon_len",
            "support_aln_ids",
        ],
    )

    gene_height, duplicon_height, te_height = track_heights
    oriented_intervals = [
        base.interval_oriented(
            base.as_int(row["effective_duplicon_raw_start"]),
            base.as_int(row["effective_duplicon_raw_end"]),
            "-",
            seq_len,
        )
        for row in intervals
    ]
    crop_start = max(0.0, min(start for start, _ in oriented_intervals) - args.flank_kb)
    crop_end = min(seq_len / 1000.0, max(end for _, end in oriented_intervals) + args.flank_kb)
    crop_width = crop_end - crop_start

    fig, ax = plt.subplots(figsize=(args.sirenia_width, args.sirenia_height))
    y = 0.0
    gene_y_offset = 0.18
    duplicon_y_offset = -0.10

    for gene in genes:
        if gene.get("gene_type") != "PGA":
            continue
        start, end = base.interval_plot(base.as_int(gene["local_start"]), base.as_int(gene["local_end"]), "-", seq_len, crop_start)
        base.draw_gene_arrow(ax, start, end, y + gene_y_offset, "+", "#2b8f7b", gene_height)

    for row in intervals:
        start, end = base.interval_plot(
            base.as_int(row["effective_duplicon_raw_start"]),
            base.as_int(row["effective_duplicon_raw_end"]),
            "-",
            seq_len,
            crop_start,
        )
        color = "#6f6f6f" if base.as_int(row["bio_order"]) % 2 else "#a8a8a8"
        base.draw_interval(ax, start, end, y + duplicon_y_offset, duplicon_height, color, 0.38)

    for row in endpoints:
        x = base.point_plot(base.as_int(row["endpoint_raw"]), "-", seq_len, crop_start)
        color = base.TE_COLORS[base.plot_te_group(row.get("te_group", ""))]
        ax.plot([x, x], [y + duplicon_y_offset - te_height / 2, y + duplicon_y_offset + te_height / 2], color=color, lw=1.2)

    ax.text(-crop_width * 0.01, y, TRICHECHUS_LABEL, va="center", ha="right", fontsize=5.5)
    ax.set_xlim(0, crop_width * 1.02)
    ax.set_ylim(-0.55, 0.55)
    ax.set_yticks([])
    ax.set_xlabel("Relative PGA locus coordinate (kb)", fontsize=8)
    ax.tick_params(axis="x", labelsize=6)
    tick_step = 25 if crop_width <= 120 else 50
    tick_max = int(crop_width // tick_step + 1) * tick_step
    ax.set_xticks(list(range(0, tick_max + 1, tick_step)))
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)

    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".png"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    track_heights = base.parse_track_heights(args.track_height_ratio)
    plot_primates_case(args, track_heights)
    plot_sirenia_case(args, track_heights)
    print(f"Wrote case panels to: {args.output_dir}")
    print("Requested primate species: " + ", ".join(PRIMATE_CASE_SPECIES))
    print("Macaca_mulatta represents Cercopithecinae; more specifically, macaques/Papionini.")


if __name__ == "__main__":
    main()
