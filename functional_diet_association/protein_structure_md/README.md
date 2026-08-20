# PGA protein-structure and molecular-dynamics workflow

This directory contains the publication workflow used to construct mature
PGA–LSFMAIPP complexes, run five independent molecular-dynamics replicates,
and identify catalysis-ready water conformations.

## Requirements

- PyRosetta 2026.3
- CHARMM-GUI
- GROMACS 2024
- Python 3 with MDAnalysis, NumPy, and pandas

Residue numbers in the scripts refer to the mature enzyme after removal of the
62-residue signal peptide and activation propeptide.

## 1. Build the three isoform models

AlphaFold 3 was used to predict mature PGA34A bound to the octapeptide
`LSFMAIPP`. Starting from that PDB, run:

```bash
python build_isoform_models.py \
  --input PGA34A_LSFMAIPP.pdb \
  --output-dir models
```

The model-building script reproduces the PyRosetta operations used in the
analysis:

- PGA34A: sequence-preserving local repacking centred on residues 160 and 203;
- PGA34B: V30L followed by local repacking centred on residue 30;
- PGA5: Q160K, A203T, and L291V followed by local repacking centred on the
  three mutated residues.

For every centre, side chains whose residue neighbour atoms are within 6 Å are
repacked with `PackRotamersMover`, `get_fa_scorefxn()`,
`restrict_to_repacking()`, and `or_include_current(True)`. The marker sequences
at mature-enzyme positions 30/160/203/291 are checked as VQAL, LQAL, and VKTV
for PGA34A, PGA34B, and PGA5, respectively. The peptide chain is independently
checked as LSFMAIPP.

## 2. Prepare each system in CHARMM-GUI

Each final complex was prepared under the following conditions:

- protonation states evaluated at pH 1.5;
- D32 manually specified as deprotonated and D215 protonated;
- CHARMM36m force field and TIP3P water;
- approximately cubic 8.40 × 8.40 × 8.40 nm³ periodic box;
- charge neutralization and approximately 0.13 M NaCl.

Place the CHARMM-GUI GROMACS inputs under a system directory:

```text
PGA34A/
└── input/
    ├── step3_input.gro
    ├── topol.top
    ├── index.ndx
    ├── toppar/
    ├── step4.0_minimization.mdp
    ├── step4.1_equilibration.mdp
    └── step5_production.mdp
```

The three MDP files in `md_example/` are the single parameter set used for all
systems and replicates. Copy them into each `input/` directory.

## 3. Run five independent MD replicates

Make GROMACS available in `PATH`, then run:

```bash
bash md_example/run_md_replicates.sh /path/to/PGA34A
```

Replicate numbers can be supplied explicitly:

```bash
bash md_example/run_md_replicates.sh /path/to/PGA34A 1 2 3 4 5
```

The script performs, in order:

1. steepest-descent minimization with positional restraints;
2. 125 ps restrained equilibration at 310 K with independently generated
   initial velocities (`gen-seed = -1`);
3. 100 ns production with a 2 fs time step at 310 K and 1 bar;
4. molecular PBC correction and centring, backbone rotational/translational
   fitting, and extraction of 20–100 ns.

Compressed coordinates are written every 250,000 production steps, equivalent
to 0.5 ns. The final analysis trajectory is therefore expected to contain 161
frames per replicate. By default, the CHARMM-GUI index groups are assumed to be
Protein=1, Backbone=4, and System=0; override `CENTER_GROUP`, `FIT_GROUP`, or
`OUTPUT_GROUP` if an index file uses different group numbers.

Hardware-specific `mdrun` options are deliberately not embedded. They can be
passed without changing the workflow, for example:

```bash
GMX_MDRUN_ARGS="-nb gpu -bonded gpu -pme gpu" \
  bash md_example/run_md_replicates.sh /path/to/PGA34A
```

## 4. Select catalysis-ready water frames

Create a tab-separated manifest with one row per system and replicate, using
`analysis/manifest.example.tsv` as a template. Run:

```bash
python analysis/select_catalysis_ready_frames.py \
  --manifest analysis/manifest.tsv \
  --output-dir analysis_results
```

For each frame, the script identifies the Phe–Met scissile bond in LSFMAIPP.
`Csc` and `Osc` are the carbonyl carbon and oxygen of the Phe residue. For each
catalytic Asp separately, the distance is measured from the same water oxygen
(`Ow`) to the nearer of that Asp residue's two carboxyl oxygens. A frame is
counted once if at least one water satisfies every condition in a threshold
set.

| Level | Asp–Ow criterion | Ow–Csc | Ow–Csc–Osc angle |
|---|---|---:|---:|
| L1 | both ≤ 3.5 Å | ≤ 3.2 Å | 100–110° |
| L2 | both ≤ 3.5 Å | ≤ 3.5 Å | 95–115° |
| L3 (primary) | both ≤ 4.0 Å | ≤ 3.5 Å | 95–120° |
| L4 | both ≤ 4.0 Å | ≤ 4.0 Å | 90–125° |
| L5 | at least one ≤ 4.0 Å and both ≤ 5.0 Å | ≤ 4.5 Å | 90–130° |

The main-text result uses L3. The script writes the selected L3 frames and
waters, counts by replicate and isoform, and the L1–L5 sensitivity table.
