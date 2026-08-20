#!/usr/bin/env python3
"""Select catalysis-ready water frames from PGA–LSFMAIPP MD trajectories."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import MDAnalysis as mda


ASP_RESIDUES = (32, 215)
MARKER_POSITIONS = (30, 160, 203, 291)
EXPECTED_MARKERS = {
    "PGA34A": "VQAL",
    "PGA34B": "LQAL",
    "PGA5": "VKTV",
}
TARGET_PEPTIDE = ("LEU", "SER", "PHE", "MET", "ALA", "ILE", "PRO", "PRO")
WATER_SELECTION = "resname SOL HOH WAT TIP3 TIP3P and name O OW OH2"


@dataclass(frozen=True)
class Threshold:
    name: str
    asp_min_max_A: float | None
    asp_max_max_A: float
    ow_csc_max_A: float
    angle_min_deg: float
    angle_max_deg: float

    def select(
        self,
        d_asp_min: np.ndarray,
        d_asp_max: np.ndarray,
        d_ow_csc: np.ndarray,
        angle: np.ndarray,
    ) -> np.ndarray:
        mask = (
            (d_asp_max <= self.asp_max_max_A)
            & (d_ow_csc <= self.ow_csc_max_A)
            & (angle >= self.angle_min_deg)
            & (angle <= self.angle_max_deg)
        )
        if self.asp_min_max_A is not None:
            mask &= d_asp_min <= self.asp_min_max_A
        return mask


THRESHOLDS = (
    Threshold("L1", None, 3.5, 3.2, 100.0, 110.0),
    Threshold("L2", None, 3.5, 3.5, 95.0, 115.0),
    Threshold("L3", None, 4.0, 3.5, 95.0, 120.0),
    Threshold("L4", None, 4.0, 4.0, 90.0, 125.0),
    Threshold("L5", 4.0, 5.0, 4.5, 90.0, 130.0),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start-ns", type=float, default=20.0)
    parser.add_argument("--end-ns", type=float, default=100.0)
    return parser.parse_args()


def read_manifest(path: Path) -> pd.DataFrame:
    separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    manifest = pd.read_csv(path, sep=separator, dtype=str)
    required = {"system", "replicate", "topology", "trajectory"}
    missing = required.difference(manifest.columns)
    if missing:
        raise ValueError(f"Manifest is missing columns: {', '.join(sorted(missing))}")
    if manifest.duplicated(["system", "replicate"]).any():
        raise ValueError("Manifest contains duplicate system/replicate rows")

    unknown = sorted(set(manifest["system"]) - set(EXPECTED_MARKERS))
    if unknown:
        raise ValueError(f"Unknown systems: {', '.join(unknown)}")

    for column in ("topology", "trajectory"):
        manifest[column] = manifest[column].map(
            lambda value: str(
                (path.parent / value).resolve()
                if not Path(value).is_absolute()
                else Path(value).resolve()
            )
        )
        absent = [value for value in manifest[column] if not Path(value).is_file()]
        if absent:
            raise FileNotFoundError(f"Missing {column}: {absent[0]}")
    return manifest


def one_letter(resname: str) -> str:
    mapping = {
        "ALA": "A",
        "GLN": "Q",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "VAL": "V",
    }
    try:
        return mapping[resname]
    except KeyError as error:
        raise ValueError(f"Unsupported residue name: {resname}") from error


def validate_markers(universe: mda.Universe, system: str) -> None:
    observed = []
    for residue in MARKER_POSITIONS:
        atoms = universe.select_atoms(f"protein and resid {residue}")
        names = sorted(set(atoms.resnames))
        if len(names) != 1:
            raise ValueError(f"Expected one protein residue at position {residue}")
        observed.append(one_letter(names[0]))
    marker_sequence = "".join(observed)
    expected = EXPECTED_MARKERS[system]
    if marker_sequence != expected:
        raise ValueError(f"{system} markers are {marker_sequence}; expected {expected}")


def find_peptide(universe: mda.Universe):
    residues = list(universe.residues)
    size = len(TARGET_PEPTIDE)
    matches = [
        residues[index : index + size]
        for index in range(len(residues) - size + 1)
        if tuple(residue.resname for residue in residues[index : index + size])
        == TARGET_PEPTIDE
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one LSFMAIPP peptide; found {len(matches)}")
    return matches[0]


def catalytic_oxygens(universe: mda.Universe, residue: int):
    atoms = universe.select_atoms(
        f"protein and resid {residue} and resname ASP ASH and "
        "name OD1 OD2 OT1 OT2 O1 O2"
    )
    if len(atoms) == 0:
        raise ValueError(f"No catalytic carboxyl oxygen found at residue {residue}")
    return atoms


def unique_atom(residue, name: str):
    atoms = residue.atoms.select_atoms(f"name {name}")
    if len(atoms) != 1:
        raise ValueError(
            f"Expected one atom named {name} in {residue.resname}{residue.resid}"
        )
    return atoms[0]


def attack_angle(
    water_positions: np.ndarray,
    carbon_position: np.ndarray,
    oxygen_position: np.ndarray,
) -> np.ndarray:
    water_vectors = water_positions - carbon_position
    oxygen_vector = oxygen_position - carbon_position
    denominator = np.clip(
        np.linalg.norm(water_vectors, axis=1) * np.linalg.norm(oxygen_vector),
        1.0e-12,
        None,
    )
    cosine = np.einsum("ij,j->i", water_vectors, oxygen_vector) / denominator
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def analyse_trajectory(
    record: pd.Series,
    start_ns: float,
    end_ns: float,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    universe = mda.Universe(record.topology, record.trajectory)
    validate_markers(universe, record.system)

    peptide = find_peptide(universe)
    scissile_carbon = unique_atom(peptide[2], "C")
    scissile_oxygen = unique_atom(peptide[2], "O")
    asp32 = catalytic_oxygens(universe, ASP_RESIDUES[0])
    asp215 = catalytic_oxygens(universe, ASP_RESIDUES[1])
    waters = universe.select_atoms(WATER_SELECTION)
    if len(waters) == 0:
        raise ValueError("No water oxygen atoms were found")

    frame_rows: list[dict[str, object]] = []
    l3_water_rows: list[dict[str, object]] = []

    for timestep in universe.trajectory:
        time_ns = float(timestep.time / 1000.0)
        if time_ns < start_ns - 1.0e-6 or time_ns > end_ns + 1.0e-6:
            continue

        positions = waters.positions
        d_asp32 = distance_array(positions, asp32.positions).min(axis=1)
        d_asp215 = distance_array(positions, asp215.positions).min(axis=1)
        d_asp_min = np.minimum(d_asp32, d_asp215)
        d_asp_max = np.maximum(d_asp32, d_asp215)
        d_ow_csc = np.linalg.norm(positions - scissile_carbon.position, axis=1)
        angles = attack_angle(
            positions,
            scissile_carbon.position,
            scissile_oxygen.position,
        )

        selections = {
            threshold.name: threshold.select(
                d_asp_min,
                d_asp_max,
                d_ow_csc,
                angles,
            )
            for threshold in THRESHOLDS
        }
        row: dict[str, object] = {
            "system": record.system,
            "replicate": record.replicate,
            "frame": int(timestep.frame),
            "time_ns": time_ns,
        }
        for name, selected in selections.items():
            count = int(np.count_nonzero(selected))
            row[f"n_{name}_waters"] = count
            row[f"hit_{name}"] = count > 0
        frame_rows.append(row)

        for index in np.flatnonzero(selections["L3"]):
            l3_water_rows.append(
                {
                    "system": record.system,
                    "replicate": record.replicate,
                    "frame": int(timestep.frame),
                    "time_ns": time_ns,
                    "water_resid": int(waters.resids[index]),
                    "water_atomid": int(waters.ids[index]),
                    "d_Asp32_A": float(d_asp32[index]),
                    "d_Asp215_A": float(d_asp215[index]),
                    "d_Ow_Csc_A": float(d_ow_csc[index]),
                    "angle_Ow_Csc_Osc_deg": float(angles[index]),
                }
            )

    return frame_rows, l3_water_rows


def summarise(frame_table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    hit_columns = [f"hit_{threshold.name}" for threshold in THRESHOLDS]
    hit_names = {
        column: column.replace("hit_", "n_") + "_frames" for column in hit_columns
    }
    by_replicate = (
        frame_table.groupby(["system", "replicate"], as_index=False)
        .agg(
            n_frames=("frame", "size"),
            **{column: (column, "sum") for column in hit_columns},
        )
        .rename(columns=hit_names)
    )
    frame_columns = list(hit_names.values())
    by_system = (
        by_replicate.groupby("system", as_index=False)
        .agg(
            n_frames=("n_frames", "sum"),
            **{column: (column, "sum") for column in frame_columns},
        )
    )
    return by_replicate, by_system


def sensitivity_table(by_system: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for threshold in THRESHOLDS:
        for record in by_system.itertuples(index=False):
            rows.append(
                {
                    "level": threshold.name,
                    "system": record.system,
                    "n_hit_frames": getattr(record, f"n_{threshold.name}_frames"),
                    "n_total_frames": record.n_frames,
                    "d_Asp_min_max_A": threshold.asp_min_max_A,
                    "d_Asp_max_max_A": threshold.asp_max_max_A,
                    "d_Ow_Csc_max_A": threshold.ow_csc_max_A,
                    "angle_min_deg": threshold.angle_min_deg,
                    "angle_max_deg": threshold.angle_max_deg,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    manifest = read_manifest(args.manifest.resolve())
    args.output_dir.mkdir(parents=True, exist_ok=True)

    frame_rows: list[dict[str, object]] = []
    water_rows: list[dict[str, object]] = []
    for record in manifest.itertuples(index=False):
        print(f"[analyse] {record.system} {record.replicate}")
        frames, waters = analyse_trajectory(record, args.start_ns, args.end_ns)
        frame_rows.extend(frames)
        water_rows.extend(waters)

    frame_table = pd.DataFrame(frame_rows)
    if frame_table.empty:
        raise ValueError("No trajectory frames fell within the requested time interval")
    frame_table = frame_table.sort_values(["system", "replicate", "time_ns"])

    l3_frames = frame_table.loc[frame_table["hit_L3"]].copy()
    l3_waters = pd.DataFrame(water_rows)
    by_replicate, by_system = summarise(frame_table)
    sensitivity = sensitivity_table(by_system)

    l3_frames.to_csv(args.output_dir / "L3_selected_frames.tsv", sep="\t", index=False)
    l3_waters.to_csv(args.output_dir / "L3_selected_waters.tsv", sep="\t", index=False)
    by_replicate.to_csv(args.output_dir / "counts_by_replicate.tsv", sep="\t", index=False)
    by_system.to_csv(args.output_dir / "counts_by_system.tsv", sep="\t", index=False)
    sensitivity.to_csv(args.output_dir / "sensitivity_L1_L5.tsv", sep="\t", index=False)

    print("\n" + by_system.to_string(index=False))
    print(f"\n[OK] results written to {args.output_dir}")


if __name__ == "__main__":
    main()
