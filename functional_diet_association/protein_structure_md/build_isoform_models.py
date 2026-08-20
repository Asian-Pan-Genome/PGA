#!/usr/bin/env python3
"""Build PGA–LSFMAIPP isoform models for MD."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import pyrosetta
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.toolbox.mutants import mutate_residue


PGA_CHAIN = "A"
PEPTIDE_CHAIN = "B"
PEPTIDE_SEQUENCE = "LSFMAIPP"
PACK_RADIUS_A = 6.0
MARKER_POSITIONS = (30, 160, 203, 291)


@dataclass(frozen=True)
class Mutation:
    residue: int
    source: str
    target: str


@dataclass(frozen=True)
class IsoformSpec:
    name: str
    markers: str
    mutations: tuple[Mutation, ...]
    repack_centres: tuple[int, ...]


ISOFORMS = (
    IsoformSpec("PGA34A", "VQAL", (), (160, 203)),
    IsoformSpec("PGA34B", "LQAL", (Mutation(30, "V", "L"),), (30,)),
    IsoformSpec(
        "PGA5",
        "VKTV",
        (
            Mutation(160, "Q", "K"),
            Mutation(203, "A", "T"),
            Mutation(291, "L", "V"),
        ),
        (160, 203, 291),
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="PGA34A–LSFMAIPP PDB")
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def pose_index(pose, chain: str, pdb_residue: int) -> int:
    index = pose.pdb_info().pdb2pose(chain, pdb_residue)
    if index == 0:
        raise ValueError(f"Residue {chain}{pdb_residue} is absent from the input PDB")
    return index


def chain_sequence(pose, chain: str) -> str:
    info = pose.pdb_info()
    return "".join(
        pose.residue(index).name1()
        for index in range(1, pose.total_residue() + 1)
        if info.chain(index) == chain
    )


def marker_sequence(pose) -> str:
    return "".join(
        pose.residue(pose_index(pose, PGA_CHAIN, residue)).name1()
        for residue in MARKER_POSITIONS
    )


def validate_complex(pose, expected_markers: str) -> None:
    observed_markers = marker_sequence(pose)
    if observed_markers != expected_markers:
        raise ValueError(
            f"Unexpected PGA markers: {observed_markers}; expected {expected_markers}"
        )
    observed_peptide = chain_sequence(pose, PEPTIDE_CHAIN)
    if observed_peptide != PEPTIDE_SEQUENCE:
        raise ValueError(
            f"Unexpected peptide sequence: {observed_peptide}; expected {PEPTIDE_SEQUENCE}"
        )


def apply_mutations(pose, mutations: tuple[Mutation, ...]) -> None:
    for mutation in mutations:
        index = pose_index(pose, PGA_CHAIN, mutation.residue)
        observed = pose.residue(index).name1()
        if observed != mutation.source:
            raise ValueError(
                f"{PGA_CHAIN}{mutation.residue} is {observed}; expected {mutation.source}"
            )
        mutate_residue(pose, index, mutation.target, pack_radius=0.0)


def local_repack(pose, centre_residues: tuple[int, ...]) -> list[int]:
    """Repack side chains within 6 Å of the selected centres."""

    centre_indices = [pose_index(pose, PGA_CHAIN, residue) for residue in centre_residues]
    selected = set(centre_indices)

    for centre_index in centre_indices:
        centre = pose.residue(centre_index).nbr_atom_xyz()
        for index in range(1, pose.total_residue() + 1):
            neighbour = pose.residue(index).nbr_atom_xyz()
            if centre.distance_squared(neighbour) <= PACK_RADIUS_A**2:
                selected.add(index)

    task = pyrosetta.standard_packer_task(pose)
    task.restrict_to_repacking()
    task.or_include_current(True)
    for index in range(1, pose.total_residue() + 1):
        if index not in selected:
            task.nonconst_residue_task(index).prevent_repacking()

    PackRotamersMover(pyrosetta.get_fa_scorefxn(), task).apply(pose)
    return sorted(selected)


def residue_labels(pose, indices: list[int]) -> list[dict[str, object]]:
    info = pose.pdb_info()
    return [
        {
            "chain": info.chain(index),
            "pdb_residue": int(info.number(index)),
            "pose_index": index,
            "residue_name": pose.residue(index).name3(),
        }
        for index in indices
    ]


def build_models(input_pdb: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    template = pose_from_pdb(str(input_pdb))
    validate_complex(template, "VQAL")

    report: dict[str, object] = {
        "input": str(input_pdb.resolve()),
        "pga_chain": PGA_CHAIN,
        "peptide_chain": PEPTIDE_CHAIN,
        "peptide_sequence": PEPTIDE_SEQUENCE,
        "pack_radius_A": PACK_RADIUS_A,
        "isoforms": [],
    }

    for spec in ISOFORMS:
        pose = template.clone()
        apply_mutations(pose, spec.mutations)
        selected = local_repack(pose, spec.repack_centres)
        validate_complex(pose, spec.markers)

        output_pdb = output_dir / f"{spec.name}_LSFMAIPP_repacked.pdb"
        pose.dump_pdb(str(output_pdb))
        report["isoforms"].append(
            {
                "name": spec.name,
                "markers": marker_sequence(pose),
                "mutations": [mutation.__dict__ for mutation in spec.mutations],
                "repack_centres": list(spec.repack_centres),
                "repacked_residues": residue_labels(pose, selected),
                "output": str(output_pdb.resolve()),
            }
        )
        print(f"[OK] {spec.name}: {output_pdb}")

    report_path = output_dir / "model_build_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"[OK] report: {report_path}")


def main() -> None:
    args = parse_args()
    pyrosetta.init("-ignore_unrecognized_res true")
    build_models(args.input, args.output_dir)


if __name__ == "__main__":
    main()
