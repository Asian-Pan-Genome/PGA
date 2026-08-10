#!/usr/bin/env python3
"""
Plot Fig.6 direction-aware PGA core / duplicon tracks and summarize endpoint TE.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon, Rectangle

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "Arial"
plt.rcParams["mathtext.it"] = "Arial:italic"
plt.rcParams["mathtext.bf"] = "Arial:bold"
plt.rcParams["axes.unicode_minus"] = False


OUTGROUP_PLOT_EXCLUDED = {"Ateles_hybridus", "Lemur_catta"}
STATS_ASSEMBLY_EXCLUDED = {"Colobus_guereza", "Macaca_nigra", "Chlorocebus_sabaeus"}
PLOT_EXCLUDED = OUTGROUP_PLOT_EXCLUDED | STATS_ASSEMBLY_EXCLUDED
STATS_EXTRA_EXCLUDED = {"Rhinopithecus_bieti"}
ALU_GROUP = "SINE/Alu"

TE_COLORS = {
    "SINE/Alu": "#e64b35",
    "LINE/L1": "#4dbbd5",
    "Other SINE": "#f39b7f",
    "Other TE": "#8491b4",
    "DNA": "#7e6148",
    "no endpoint TE": "#bdbdbd",
}


def parse_args() -> argparse.Namespace:
    base_dir = Path(__file__).resolve().parent
    default_input = base_dir / "untangle_duplicon_TE"
    parser = argparse.ArgumentParser(
        description="Draw 25-species PGA core/duplicon tracks and filtered endpoint Alu summaries."
    )
    parser.add_argument("--base-dir", type=Path, default=base_dir)
    parser.add_argument("--input-dir", type=Path, default=default_input)
    parser.add_argument("--species-fa", type=Path, default=base_dir / "apes_owms.fa")
    parser.add_argument("--tree", type=Path, default=base_dir.parent / "295.sp.tree")
    parser.add_argument("--reference", default="Ateles_hybridus")
    parser.add_argument("--m", default="256")
    parser.add_argument("--output-prefix", type=Path, default=None)
    parser.add_argument("--max-join-transition-len", type=int, default=5000)
    parser.add_argument("--flank-kb", type=float, default=10)
    parser.add_argument("--track-height-ratio", default="2:1:1")
    parser.add_argument("--plot-width", type=float, default=8.27)
    parser.add_argument("--plot-height", type=float, default=5.85)
    return parser.parse_args()


def read_graph_species(fasta: Path) -> list[str]:
    names: list[str] = []
    with fasta.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                if name not in names:
                    names.append(name)
    return names


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def prefixed_path(prefix: Path, suffix: str) -> Path:
    return prefix.parent / f"{prefix.name}{suffix}"


def as_int(value: str | int | float | None, default: int = 0) -> int:
    if value in (None, ""):
        return default
    return int(float(value))


def as_float(value: str | int | float | None, default: float = 0.0) -> float:
    if value in (None, ""):
        return default
    return float(value)


def interval_oriented(start: int, end: int, strand: str, seq_len: int) -> tuple[float, float]:
    if strand == "-":
        left = seq_len - end
        right = seq_len - start
    else:
        left = start
        right = end
    if right < left:
        left, right = right, left
    return left / 1000.0, right / 1000.0


def point_oriented(pos: int, strand: str, seq_len: int) -> float:
    return (seq_len - pos if strand == "-" else pos) / 1000.0


def parse_track_heights(spec: str) -> tuple[float, float, float]:
    parts = spec.split(":")
    if len(parts) != 3:
        raise SystemExit("--track-height-ratio must be formatted as gene:duplicon:TE, for example 2:1:1")
    ratios = [float(part) for part in parts]
    if any(value <= 0 for value in ratios):
        raise SystemExit("--track-height-ratio values must be positive")
    max_height = 0.30
    scale = max_height / max(ratios)
    return tuple(value * scale for value in ratios)  # type: ignore[return-value]


def interval_plot(start: int, end: int, strand: str, seq_len: int, crop_start_kb: float) -> tuple[float, float]:
    left, right = interval_oriented(start, end, strand, seq_len)
    return left - crop_start_kb, right - crop_start_kb


def point_plot(pos: int, strand: str, seq_len: int, crop_start_kb: float) -> float:
    return point_oriented(pos, strand, seq_len) - crop_start_kb


def plot_te_group(group: str) -> str:
    if not group:
        return "no endpoint TE"
    if group in TE_COLORS:
        return group
    if group in {"Other LINE", "LINE/L2", "LTR/ERV", "Retroposon/SVA"}:
        return "Other TE"
    return "Other TE"


def load_pruned_tree(tree_path: Path, plot_species: list[str]):
    try:
        from Bio import Phylo
    except ImportError as exc:
        raise SystemExit("Biopython is required to draw the tree.") from exc

    tree = Phylo.read(str(tree_path), "newick")
    keep = set(plot_species)
    for terminal in list(tree.get_terminals()):
        if terminal.name not in keep:
            tree.prune(terminal)
    tree_order = [term.name for term in tree.get_terminals() if term.name in keep]
    missing = [sp for sp in plot_species if sp not in tree_order]
    return tree, tree_order + missing


def draw_tree(ax, tree, species_order: list[str], y_by_species: dict[str, float]) -> None:
    depths = tree.depths()
    if not any(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)

    y_cache = {}

    def clade_y(clade):
        if clade in y_cache:
            return y_cache[clade]
        if clade.is_terminal():
            y_cache[clade] = y_by_species.get(clade.name, 0.0)
        else:
            child_y = [clade_y(child) for child in clade.clades]
            y_cache[clade] = sum(child_y) / len(child_y)
        return y_cache[clade]

    def draw_clade(clade):
        x = depths[clade]
        if clade.clades:
            child_ys = []
            for child in clade.clades:
                cx = depths[child]
                cy = clade_y(child)
                child_ys.append(cy)
                ax.plot([x, cx], [cy, cy], color="black", lw=0.6, solid_capstyle="butt")
                draw_clade(child)
            ax.plot([x, x], [min(child_ys), max(child_ys)], color="black", lw=0.6, solid_capstyle="butt")

    draw_clade(tree.root)
    max_depth = max(depths.values()) if depths else 1
    ax.set_xlim(-max_depth * 0.02, max_depth * 1.05)
    ax.set_ylim(-0.7, len(species_order) - 0.3)
    ax.axis("off")


def grouped_by(rows: list[dict[str, str]], *keys: str) -> dict[tuple, list[dict[str, str]]]:
    out: dict[tuple, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        out[tuple(row.get(key, "") for key in keys)].append(row)
    return out


def endpoint_status(rows: list[dict[str, str]]) -> tuple[bool, str, str]:
    hit_groups = [row.get("te_group", "") for row in rows if row.get("te_hit") == "1" and row.get("te_group")]
    if not hit_groups:
        return False, "no endpoint TE", "noTE"
    if ALU_GROUP in hit_groups:
        return True, ALU_GROUP, ALU_GROUP
    groups = sorted(set(hit_groups))
    simple = "+".join(groups)
    return False, plot_te_group(groups[0]), simple


def build_effective_duplicons(
    duplicon_rows: list[dict[str, str]],
    core_endpoint_rows: list[dict[str, str]],
    duplicon_endpoint_rows: list[dict[str, str]],
    reference: str,
    m_value: str,
    max_join_transition_len: int,
) -> tuple[list[dict], list[dict]]:
    selected_dups = [
        row for row in duplicon_rows if row.get("reference_name") == reference and row.get("m") == m_value
    ]
    core_endpoint_by = grouped_by(
        [row for row in core_endpoint_rows if row.get("reference_name") == reference and row.get("m") == m_value],
        "core_id",
        "endpoint_role",
    )
    duplicon_endpoint_by = grouped_by(
        [row for row in duplicon_endpoint_rows if row.get("reference_name") == reference and row.get("m") == m_value],
        "duplicon_id",
        "endpoint_role",
    )

    effective_rows: list[dict] = []
    effective_endpoint_rows: list[dict] = []
    for row in selected_dups:
        transition_len = as_int(row.get("transition_len"), 0)
        has_transition = row.get("transition_id", "") != ""
        join_transition = has_transition and transition_len <= max_join_transition_len
        if has_transition and not join_transition:
            eff_start = as_int(row["core_raw_start"])
            eff_end = as_int(row["core_raw_end"])
            join_status = "not_joined_long_transition"
            endpoint_source = "core_endpoint"
            endpoint_lookup = core_endpoint_by
            endpoint_key = row["core_id"]
        else:
            eff_start = as_int(row["duplicon_raw_start"])
            eff_end = as_int(row["duplicon_raw_end"])
            join_status = "joined_transition" if has_transition else "terminal_core_only"
            endpoint_source = "duplicon_endpoint"
            endpoint_lookup = duplicon_endpoint_by
            endpoint_key = row["duplicon_id"]
        if eff_end < eff_start:
            eff_start, eff_end = eff_end, eff_start
        out = {
            **row,
            "effective_duplicon_raw_start": eff_start,
            "effective_duplicon_raw_end": eff_end,
            "effective_duplicon_len": eff_end - eff_start,
            "join_status": join_status,
            "endpoint_source": endpoint_source,
            "max_join_transition_len": max_join_transition_len,
        }
        effective_rows.append(out)
        for role in ("biological_start", "biological_end"):
            for endpoint_row in endpoint_lookup.get((endpoint_key, role), []):
                copied = {**endpoint_row}
                copied["duplicon_id"] = row["duplicon_id"]
                copied["core_id"] = row["core_id"]
                copied["endpoint_source"] = endpoint_source
                copied["join_status"] = join_status
                copied["max_join_transition_len"] = max_join_transition_len
                effective_endpoint_rows.append(copied)
    return effective_rows, effective_endpoint_rows


def build_crop_windows(
    species_order: list[str],
    effective_rows: list[dict],
    genes_by_species: dict[str, list[dict[str, str]]],
    length_by_species: dict[str, int],
    pga_strand_by_species: dict[str, str],
    flank_kb: float,
) -> dict[str, tuple[float, float]]:
    windows: dict[str, tuple[float, float]] = {}
    intervals_by_species: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in effective_rows:
        species = row["query_species"]
        seq_len = length_by_species.get(species, 0)
        if not seq_len:
            continue
        left, right = interval_oriented(
            as_int(row["effective_duplicon_raw_start"]),
            as_int(row["effective_duplicon_raw_end"]),
            row.get("pga_strand", "+"),
            seq_len,
        )
        intervals_by_species[species].append((left, right))

    for species in species_order:
        seq_len = length_by_species.get(species, 0)
        strand = pga_strand_by_species.get(species, "+")
        intervals = intervals_by_species.get(species, [])
        if not intervals:
            for gene in genes_by_species.get(species, []):
                intervals.append(
                    interval_oriented(as_int(gene["local_start"]), as_int(gene["local_end"]), strand, seq_len)
                )
        if intervals:
            left = max(0.0, min(start for start, _ in intervals) - flank_kb)
            right = min(seq_len / 1000.0 if seq_len else max(end for _, end in intervals) + flank_kb, max(end for _, end in intervals) + flank_kb)
            if right <= left:
                right = left + 1.0
            windows[species] = (left, right)
    return windows


def build_endpoint_stats(
    endpoint_rows: list[dict[str, str]],
    duplicon_rows: list[dict[str, str]],
    stats_species: set[str],
    reference: str,
    m_value: str,
) -> tuple[list[dict], list[dict], list[dict]]:
    selected_dups = [
        row
        for row in duplicon_rows
        if row.get("reference_name") == reference and row.get("m") == m_value and row.get("query_species") in stats_species
    ]
    selected_ids = {row["duplicon_id"] for row in selected_dups}
    endpoint_groups = grouped_by(
        [
            row
            for row in endpoint_rows
            if row.get("reference_name") == reference
            and row.get("m") == m_value
            and row.get("query_species") in stats_species
            and row.get("duplicon_id") in selected_ids
        ],
        "duplicon_id",
        "endpoint_role",
    )

    per_dup_rows: list[dict] = []
    config_counts: dict[tuple[str, str], int] = defaultdict(int)
    both_alu = 0
    at_least_one_alu = 0
    zero_alu = 0
    for dup in selected_dups:
        did = dup["duplicon_id"]
        start_has, start_plot, start_simple = endpoint_status(endpoint_groups.get((did, "biological_start"), []))
        end_has, end_plot, end_simple = endpoint_status(endpoint_groups.get((did, "biological_end"), []))
        n_alu = int(start_has) + int(end_has)
        both_alu += int(n_alu == 2)
        at_least_one_alu += int(n_alu >= 1)
        zero_alu += int(n_alu == 0)
        config_counts[(start_simple, end_simple)] += 1
        per_dup_rows.append(
            {
                "reference_name": reference,
                "m": m_value,
                "query_species": dup["query_species"],
                "gene": dup["gene"],
                "duplicon_id": did,
                "pga_strand": dup["pga_strand"],
                "join_status": dup.get("join_status", ""),
                "endpoint_source": dup.get("endpoint_source", ""),
                "transition_len": dup.get("transition_len", ""),
                "biological_start_group": start_simple,
                "biological_end_group": end_simple,
                "both_endpoints_alu": int(n_alu == 2),
                "at_least_one_endpoint_alu": int(n_alu >= 1),
                "zero_alu_endpoint": int(n_alu == 0),
            }
        )

    total = len(selected_dups)
    summary_rows = [
        {
            "reference_name": reference,
            "m": m_value,
            "stats_species_n": len(stats_species),
            "duplicon_n": total,
            "both_endpoints_alu_n": both_alu,
            "both_endpoints_alu_frac": f"{both_alu / total:.6f}" if total else "",
            "at_least_one_endpoint_alu_n": at_least_one_alu,
            "at_least_one_endpoint_alu_frac": f"{at_least_one_alu / total:.6f}" if total else "",
            "zero_alu_endpoint_n": zero_alu,
            "zero_alu_endpoint_frac": f"{zero_alu / total:.6f}" if total else "",
        }
    ]
    config_rows = [
        {"biological_start_group": key[0], "biological_end_group": key[1], "n": value}
        for key, value in sorted(config_counts.items(), key=lambda item: (-item[1], item[0]))
    ]
    return summary_rows, config_rows, per_dup_rows


def draw_gene_arrow(ax, start: float, end: float, y: float, strand: str, color: str, height: float) -> None:
    width = end - start
    if width <= 0:
        return
    body_height = height * (2.0 / 3.0)
    head_height = height
    head_length = min(3.0, max(0.8, width * 0.24))
    head_length = min(head_length, width * 0.65)
    if strand == "-":
        points = [
            (end, y - body_height / 2),
            (start + head_length, y - body_height / 2),
            (start + head_length, y - head_height / 2),
            (start, y),
            (start + head_length, y + head_height / 2),
            (start + head_length, y + body_height / 2),
            (end, y + body_height / 2),
        ]
    else:
        points = [
            (start, y - body_height / 2),
            (end - head_length, y - body_height / 2),
            (end - head_length, y - head_height / 2),
            (end, y),
            (end - head_length, y + head_height / 2),
            (end - head_length, y + body_height / 2),
            (start, y + body_height / 2),
        ]
    ax.add_patch(
        Polygon(
            points,
            closed=True,
            facecolor="white",
            edgecolor=color,
            linewidth=0.5,
            joinstyle="miter",
        )
    )


def draw_interval(ax, start: float, end: float, y: float, height: float, color: str, alpha: float, edge: str | None = None) -> None:
    if end <= start:
        return
    ax.add_patch(
        Rectangle(
            (start, y - height / 2),
            end - start,
            height,
            facecolor=color,
            edgecolor=edge if edge else "none",
            linewidth=0.35 if edge else 0,
            alpha=alpha,
        )
    )


def draw_endpoint_ticks(
    ax,
    endpoint_rows: list[dict[str, str]],
    y_by_species: dict[str, float],
    length_by_species: dict[str, int],
    crop_by_species: dict[str, tuple[float, float]],
    y_offset: float,
    height: float,
) -> None:
    for row in endpoint_rows:
        species = row["query_species"]
        if species not in y_by_species or species not in crop_by_species:
            continue
        if row.get("te_hit") != "1" or not row.get("te_group"):
            continue
        seq_len = length_by_species.get(species, 0)
        if not seq_len:
            continue
        x = point_plot(as_int(row["endpoint_raw"]), row.get("pga_strand", "+"), seq_len, crop_by_species[species][0])
        group = plot_te_group(row.get("te_group", ""))
        color = TE_COLORS[group]
        y = y_by_species[species] + y_offset
        ax.plot([x, x], [y - height / 2, y + height / 2], color=color, lw=1.2)


def draw_transition_te(
    ax,
    transition_te_rows: list[dict[str, str]],
    y_by_species: dict[str, float],
    length_by_species: dict[str, int],
    crop_by_species: dict[str, tuple[float, float]],
    y_offset: float,
    height: float,
) -> None:
    for row in transition_te_rows:
        if row.get("te_hit") != "1":
            continue
        species = row["query_species"]
        if species not in y_by_species or species not in crop_by_species:
            continue
        seq_len = length_by_species.get(species, 0)
        if not seq_len:
            continue
        strand = row.get("pga_strand", "+")
        start, end = interval_plot(as_int(row["te_start0"]), as_int(row["te_end0"]), strand, seq_len, crop_by_species[species][0])
        color = TE_COLORS[plot_te_group(row.get("te_group", ""))]
        y = y_by_species[species] + y_offset
        if end - start < 0.08:
            ax.plot([start, start], [y - height / 2, y + height / 2], color=color, lw=0.8)
        else:
            draw_interval(ax, start, end, y, height, color, 0.85)


def plot_tracks(
    mode: str,
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
    transition_rows: list[dict[str, str]],
    transition_te_rows: list[dict[str, str]],
    reference: str,
    m_value: str,
    plot_width: float,
    plot_height: float,
) -> None:
    n = len(species_order)
    y_by_species = {species: n - 1 - idx for idx, species in enumerate(species_order)}
    gene_height, duplicon_height, te_height = track_heights
    gene_y_offset = 0.18
    duplicon_y_offset = -0.10
    fig = plt.figure(figsize=(plot_width, plot_height))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.85, 1.65, 5.9], wspace=0.0)
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_label = fig.add_subplot(gs[0, 1], sharey=ax_tree)
    ax = fig.add_subplot(gs[0, 2], sharey=ax_tree)

    draw_tree(ax_tree, tree, species_order, y_by_species)
    ax_label.set_xlim(0, 1)
    ax_label.set_ylim(-0.7, n - 0.3)
    ax_label.axis("off")
    for species in species_order:
        ax_label.text(0.98, y_by_species[species], species, va="center", ha="right", fontsize=5.2)
    max_x = max((crop_by_species.get(sp, (0.0, 1.0))[1] - crop_by_species.get(sp, (0.0, 1.0))[0] for sp in species_order), default=1.0)

    for species in species_order:
        y = y_by_species[species]
        strand = pga_strand_by_species.get(species, "+")
        seq_len = length_by_species.get(species, 0)
        crop_start = crop_by_species.get(species, (0.0, max_x))[0]
        for gene in genes_by_species.get(species, []):
            start, end = interval_plot(as_int(gene["local_start"]), as_int(gene["local_end"]), strand, seq_len, crop_start)
            draw_gene_arrow(ax, start, end, y + gene_y_offset, "+", "#2b8f7b", gene_height)

    if mode == "core":
        for row in transition_rows:
            species = row["query_species"]
            if species not in y_by_species or species not in crop_by_species:
                continue
            seq_len = length_by_species.get(species, 0)
            if not seq_len:
                continue
            start, end = interval_plot(
                as_int(row["transition_raw_start"]),
                as_int(row["transition_raw_end"]),
                row.get("pga_strand", "+"),
                seq_len,
                crop_by_species[species][0],
            )
            draw_interval(ax, start, end, y_by_species[species] + duplicon_y_offset, duplicon_height, "#d7d7d7", 0.45)

    for row in intervals:
        species = row["query_species"]
        if species not in y_by_species or species not in crop_by_species:
            continue
        seq_len = length_by_species.get(species, 0)
        if not seq_len:
            continue
        strand = row.get("pga_strand", "+")
        if mode == "core":
            start = as_int(row["core_raw_start"])
            end = as_int(row["core_raw_end"])
            color = "#5aa6c8"
            alpha = 0.38
            height = duplicon_height
        else:
            start = as_int(row["effective_duplicon_raw_start"])
            end = as_int(row["effective_duplicon_raw_end"])
            bio_order = as_int(row.get("bio_order", 0))
            color = "#6f6f6f" if bio_order % 2 else "#a8a8a8"
            alpha = 0.38
            height = duplicon_height
        left, right = interval_plot(start, end, strand, seq_len, crop_by_species[species][0])
        draw_interval(ax, left, right, y_by_species[species] + duplicon_y_offset, height, color, alpha)

    if mode == "core":
        draw_transition_te(
            ax,
            transition_te_rows,
            y_by_species,
            length_by_species,
            crop_by_species,
            duplicon_y_offset,
            te_height,
        )
    draw_endpoint_ticks(ax, endpoint_rows, y_by_species, length_by_species, crop_by_species, duplicon_y_offset, te_height)

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

    legend_items = [
        Patch(facecolor="#5aa6c8", alpha=0.38, edgecolor="none", label="duplicon core"),
        Patch(facecolor="#d7d7d7", alpha=0.45, edgecolor="none", label="transition"),
    ]
    if mode == "duplicon":
        legend_items = [
            Patch(facecolor="#6f6f6f", alpha=0.38, edgecolor="none", label="duplicon"),
            Patch(facecolor="#a8a8a8", alpha=0.38, edgecolor="none", label="alternate duplicon"),
        ]
    legend_items.extend(
        [
            Line2D([0], [0], color=TE_COLORS["SINE/Alu"], lw=1.5, label="SINE/Alu"),
            Line2D([0], [0], color=TE_COLORS["LINE/L1"], lw=1.5, label="LINE/L1"),
            Line2D([0], [0], color=TE_COLORS["Other SINE"], lw=1.5, label="Other SINE"),
            Line2D([0], [0], color=TE_COLORS["Other TE"], lw=1.5, label="Other TE"),
            Line2D([0], [0], color=TE_COLORS["DNA"], lw=1.5, label="DNA"),
        ]
    )
    ax.legend(handles=legend_items, loc="upper right", ncol=1, frameon=False, fontsize=6)
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_dir: Path = args.input_dir
    output_prefix = args.output_prefix or (input_dir / f"fig6_{args.reference}.m{args.m}")

    graph_species = read_graph_species(args.species_fa)
    plot_species = [sp for sp in graph_species if sp not in PLOT_EXCLUDED]
    stats_excluded = PLOT_EXCLUDED | STATS_EXTRA_EXCLUDED | STATS_ASSEMBLY_EXCLUDED
    stats_species = [sp for sp in plot_species if sp not in stats_excluded]
    plot_species_set = set(plot_species)
    stats_species_set = set(stats_species)
    track_heights = parse_track_heights(args.track_height_ratio)

    gene_rows = read_tsv(input_dir / "untangle_PGA_genes.local.tsv")
    core_rows = read_tsv(input_dir / "untangle_core_intervals.tsv")
    transition_rows = read_tsv(input_dir / "untangle_transition_regions.tsv")
    duplicon_rows = read_tsv(input_dir / "untangle_duplicon_intervals.tsv")
    core_endpoint_rows = read_tsv(input_dir / "untangle_core_endpoint_TE.tsv")
    transition_te_rows = read_tsv(input_dir / "untangle_transition_region_TE.tsv")
    duplicon_endpoint_rows = read_tsv(input_dir / "untangle_duplicon_endpoint_TE.tsv")
    run_summary_rows = read_tsv(input_dir / "untangle_run_summary.tsv")
    effective_duplicon_rows, effective_endpoint_rows = build_effective_duplicons(
        duplicon_rows,
        core_endpoint_rows,
        duplicon_endpoint_rows,
        args.reference,
        args.m,
        args.max_join_transition_len,
    )

    tree, tree_order = load_pruned_tree(args.tree, plot_species)
    species_order = [sp for sp in tree_order if sp in plot_species_set]
    for sp in plot_species:
        if sp not in species_order:
            species_order.append(sp)

    genes_by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    pga_strand_by_species: dict[str, str] = {}
    length_by_species: dict[str, int] = {}
    pga_cn_by_species: dict[str, int] = {}
    for row in gene_rows:
        species = row["species"]
        if species not in plot_species_set:
            continue
        genes_by_species[species].append(row)
        pga_strand_by_species[species] = row["pga_strand"]
        length_by_species[species] = max(length_by_species.get(species, 0), as_int(row["anchor_end"]) - as_int(row["anchor_start"]))
        pga_cn_by_species[species] = as_int(row["pga_cn"])
    for species in genes_by_species:
        genes_by_species[species].sort(key=lambda row: as_float(row["oriented_start"]))

    selected_core = [
        row
        for row in core_rows
        if row["reference_name"] == args.reference and row["m"] == args.m and row["query_species"] in plot_species_set
    ]
    selected_transitions = [
        row
        for row in transition_rows
        if row["reference_name"] == args.reference and row["m"] == args.m and row["query_species"] in plot_species_set
    ]
    selected_dup = [
        row
        for row in effective_duplicon_rows
        if row["reference_name"] == args.reference and row["m"] == args.m and row["query_species"] in plot_species_set
    ]
    selected_core_endpoint = [
        row
        for row in core_endpoint_rows
        if row["reference_name"] == args.reference and row["m"] == args.m and row["query_species"] in plot_species_set
    ]
    selected_transition_te = [
        row
        for row in transition_te_rows
        if row["reference_name"] == args.reference and row["m"] == args.m and row["query_species"] in plot_species_set
    ]
    selected_dup_endpoint = [
        row
        for row in effective_endpoint_rows
        if row["reference_name"] == args.reference and row["m"] == args.m and row["query_species"] in plot_species_set
    ]
    crop_by_species = build_crop_windows(
        species_order,
        selected_dup,
        genes_by_species,
        length_by_species,
        pga_strand_by_species,
        args.flank_kb,
    )

    stats_rows, config_rows, per_dup_rows = build_endpoint_stats(
        effective_endpoint_rows, effective_duplicon_rows, stats_species_set, args.reference, args.m
    )
    duplicon_n = as_int(stats_rows[0]["duplicon_n"]) if stats_rows else 0
    not_joined_n = sum(
        1
        for row in effective_duplicon_rows
        if row.get("reference_name") == args.reference
        and row.get("m") == args.m
        and row.get("query_species") in plot_species_set
        and row.get("join_status") == "not_joined_long_transition"
    )
    scope_row = {
        "reference_name": args.reference,
        "m": args.m,
        "graph_species_n": len(graph_species),
        "plot_species_n": len(plot_species),
        "stats_species_n": len(stats_species),
        "duplicon_n": duplicon_n,
        "max_join_transition_len": args.max_join_transition_len,
        "not_joined_long_transition_n": not_joined_n,
        "excluded_from_plot": ",".join(sorted(PLOT_EXCLUDED)),
        "excluded_from_stats": ",".join(sorted(stats_excluded)),
        "assembly_or_alignment_suspect_excluded_from_stats": ",".join(sorted(STATS_ASSEMBLY_EXCLUDED)),
        "both_endpoints_alu_n": stats_rows[0]["both_endpoints_alu_n"],
        "both_endpoints_alu_frac": stats_rows[0]["both_endpoints_alu_frac"],
        "at_least_one_endpoint_alu_n": stats_rows[0]["at_least_one_endpoint_alu_n"],
        "at_least_one_endpoint_alu_frac": stats_rows[0]["at_least_one_endpoint_alu_frac"],
        "zero_alu_endpoint_n": stats_rows[0]["zero_alu_endpoint_n"],
        "zero_alu_endpoint_frac": stats_rows[0]["zero_alu_endpoint_frac"],
    }

    summary_by_species = {
        row["query_species"]: row
        for row in run_summary_rows
        if row["reference_name"] == args.reference and row["m"] == args.m
    }
    species_scope_rows = []
    for species in graph_species:
        summary = summary_by_species.get(species, {})
        species_scope_rows.append(
            {
                "species": species,
                "in_plot": int(species in plot_species_set),
                "in_stats": int(species in stats_species_set),
                "exclusion_reason": (
                    "assembly_or_alignment_suspect_incomplete_gene_wrapping"
                    if species in STATS_ASSEMBLY_EXCLUDED
                    else "outgroup_plot_and_stats"
                    if species in OUTGROUP_PLOT_EXCLUDED
                    else "single_copy_OWM_stats_only"
                    if species in STATS_EXTRA_EXCLUDED
                    else ""
                ),
                "pga_cn": pga_cn_by_species.get(species, summary.get("expected_pga_cn", "")),
                "pga_strand": pga_strand_by_species.get(species, summary.get("pga_strand", "")),
                "status": summary.get("status", ""),
            }
        )

    write_tsv(
        prefixed_path(output_prefix, ".scope_summary.tsv"),
        [scope_row],
        [
            "reference_name",
            "m",
            "graph_species_n",
            "plot_species_n",
            "stats_species_n",
            "duplicon_n",
            "max_join_transition_len",
            "not_joined_long_transition_n",
            "excluded_from_plot",
            "excluded_from_stats",
            "assembly_or_alignment_suspect_excluded_from_stats",
            "both_endpoints_alu_n",
            "both_endpoints_alu_frac",
            "at_least_one_endpoint_alu_n",
            "at_least_one_endpoint_alu_frac",
            "zero_alu_endpoint_n",
            "zero_alu_endpoint_frac",
        ],
    )
    write_tsv(
        prefixed_path(output_prefix, ".species_scope.tsv"),
        species_scope_rows,
        ["species", "in_plot", "in_stats", "exclusion_reason", "pga_cn", "pga_strand", "status"],
    )
    write_tsv(
        prefixed_path(output_prefix, ".effective_duplicon_intervals.tsv"),
        selected_dup,
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
            "effective_duplicon_raw_start",
            "effective_duplicon_raw_end",
            "effective_duplicon_len",
            "join_status",
            "endpoint_source",
            "max_join_transition_len",
            "duplicon_raw_start",
            "duplicon_raw_end",
            "duplicon_len",
            "core_raw_start",
            "core_raw_end",
            "transition_raw_start",
            "transition_raw_end",
            "transition_len",
        ],
    )
    write_tsv(
        prefixed_path(output_prefix, ".effective_duplicon_endpoint_TE.tsv"),
        selected_dup_endpoint,
        [
            "untangle_file",
            "reference_name",
            "m",
            "query_species",
            "gene",
            "pga_strand",
            "core_id",
            "duplicon_id",
            "endpoint_role",
            "endpoint_source",
            "join_status",
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
    write_tsv(
        prefixed_path(output_prefix, ".duplicon_endpoint_Alu_stats.tsv"),
        stats_rows,
        [
            "reference_name",
            "m",
            "stats_species_n",
            "duplicon_n",
            "both_endpoints_alu_n",
            "both_endpoints_alu_frac",
            "at_least_one_endpoint_alu_n",
            "at_least_one_endpoint_alu_frac",
            "zero_alu_endpoint_n",
            "zero_alu_endpoint_frac",
        ],
    )
    write_tsv(
        prefixed_path(output_prefix, ".endpoint_configuration_counts.tsv"),
        config_rows,
        ["biological_start_group", "biological_end_group", "n"],
    )
    write_tsv(
        prefixed_path(output_prefix, ".duplicon_endpoint_classification.tsv"),
        per_dup_rows,
        [
            "reference_name",
            "m",
            "query_species",
            "gene",
            "duplicon_id",
            "pga_strand",
            "join_status",
            "endpoint_source",
            "transition_len",
            "biological_start_group",
            "biological_end_group",
            "both_endpoints_alu",
            "at_least_one_endpoint_alu",
            "zero_alu_endpoint",
        ],
    )

    plot_tracks(
        "core",
        prefixed_path(output_prefix, ".core_TE_tracks.pdf"),
        prefixed_path(output_prefix, ".core_TE_tracks.png"),
        tree,
        species_order,
        genes_by_species,
        length_by_species,
        pga_strand_by_species,
        crop_by_species,
        track_heights,
        selected_core,
        selected_core_endpoint,
        selected_transitions,
        selected_transition_te,
        args.reference,
        args.m,
        args.plot_width,
        args.plot_height,
    )
    plot_tracks(
        "duplicon",
        prefixed_path(output_prefix, ".duplicon_endpoint_TE_tracks.pdf"),
        prefixed_path(output_prefix, ".duplicon_endpoint_TE_tracks.png"),
        tree,
        species_order,
        genes_by_species,
        length_by_species,
        pga_strand_by_species,
        crop_by_species,
        track_heights,
        selected_dup,
        selected_dup_endpoint,
        selected_transitions,
        selected_transition_te,
        args.reference,
        args.m,
        args.plot_width,
        args.plot_height,
    )

    print(f"Wrote prefix: {output_prefix}")
    print(f"Plot species: {len(plot_species)}")
    print(f"Stats species: {len(stats_species)}")
    print(f"Duplicons in stats: {duplicon_n}")
    print(f"Long transitions not joined: {not_joined_n}")
    print(f"Both endpoints Alu: {scope_row['both_endpoints_alu_n']} ({float(scope_row['both_endpoints_alu_frac']):.1%})")
    print(
        f"At least one endpoint Alu: {scope_row['at_least_one_endpoint_alu_n']} "
        f"({float(scope_row['at_least_one_endpoint_alu_frac']):.1%})"
    )


if __name__ == "__main__":
    main()
