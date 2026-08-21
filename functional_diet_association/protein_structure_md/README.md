# PGA protein-structure and molecular-dynamics workflow

Commands and parameters used for isoform model construction, molecular-dynamics simulation, and catalysis-ready water analysis of PGA–LSFMAIPP complexes.

## Requirements

- AlphaFold 3
- PyRosetta 2026.3
- CHARMM-GUI
- GROMACS 2024
- Python 3 with MDAnalysis, NumPy, and pandas

Residue numbering below refers to the mature enzyme.

## 1. Isoform model construction

The mature PGA34A–LSFMAIPP complex was predicted with AlphaFold 3 after removing the N-terminal 62-residue signal peptide and activation propeptide. PGA34B and PGA5 models were generated from the PGA34A complex with PyRosetta.

Run:

```bash
python build_isoform_models.py \
    --input PGA34A_LSFMAIPP.pdb \
    --output-dir models
```

The model-building script applies the following substitutions and local repacking steps:

```text
PGA34A: no sequence change; repack within 6 Å of residues 160 and 203
PGA34B: V30L; repack within 6 Å of residue 30
PGA5:   Q160K, A203T and L291V; repack within 6 Å of residues 160, 203 and 291
```

Side-chain repacking uses `PackRotamersMover`, `get_fa_scorefxn()`, `restrict_to_repacking()`, and `or_include_current(True)`.

## 2. CHARMM-GUI system preparation

Systems were prepared in CHARMM-GUI with:

```text
pH                              1.5
D32                             deprotonated
D215                            protonated
force field                     CHARMM36m
water model                     TIP3P
periodic box                    ~8.40 × 8.40 × 8.40 nm³
salt                            ~0.13 M NaCl after neutralization
```

For each isoform, retain the CHARMM-GUI topology, coordinates, index file, and `toppar/` directory, and add the MDP files supplied in this repository:

```text
SYSTEM/
└── input/
    ├── step3_input.gro
    ├── topol.top
    ├── index.ndx
    ├── toppar/
    ├── step4.0_minimization.mdp
    ├── step4.1_equilibration.mdp
    └── step5_production.mdp
```

Production simulations use:

```text
dt                      0.002 ps
nsteps                  50,000,000
simulation length       100 ns
temperature             310 K
pressure                1 bar
thermostat              velocity rescaling
barostat                isotropic C-rescale
Coulomb treatment       PME, 1.2-nm cutoff
van der Waals            force switch, 1.0–1.2 nm
constraints             LINCS on bonds involving H
```

## 3. Molecular-dynamics replicates

Five independent replicates were generated for each isoform:

```bash
bash run_md_replicates.sh /path/to/SYSTEM 1 2 3 4 5
```

Each replicate includes restrained energy minimization, restrained equilibration at 310 K, and 100-ns production at 310 K and 1 bar. Independent initial velocities are generated with `gen-seed = -1`.

Trajectories are processed with `-pbc mol -center` and backbone `-fit rot+trans`. Frames from 20–100 ns are retained at 0.5-ns intervals, yielding 161 analysed frames per replicate.

Optional GROMACS hardware arguments can be passed with:

```bash
GMX_MDRUN_ARGS="-nb gpu -bonded gpu -pme gpu" \
    bash run_md_replicates.sh /path/to/SYSTEM 1 2 3 4 5
```

## 4. Catalysis-ready water analysis

Arrange processed trajectories as:

```text
MD_DATA/
├── PGA34A/rep1 ... rep5/
├── PGA34B/rep1 ... rep5/
└── PGA5/rep1 ... rep5/
```

Each replicate directory must contain:

```text
step5_ref.tpr
md_fit_20_100ns.xtc
```

Run:

```bash
python select_catalysis_ready_frames.py \
    --data-dir /path/to/MD_DATA \
    --output-dir analysis_results
```

`Csc` and `Osc` denote the carbonyl carbon and oxygen of Phe at the Phe–Met scissile bond. Each Asp–Ow distance is measured to the nearer carboxyl oxygen of D32 or D215. A frame is counted once if at least one water molecule satisfies all criteria at the selected threshold.

| Level | Asp–Ow criterion | Ow–Csc | Ow–Csc–Osc angle |
| --- | --- | ---: | ---: |
| L1 | both ≤ 3.5 Å | ≤ 3.2 Å | 100–110° |
| L2 | both ≤ 3.5 Å | ≤ 3.5 Å | 95–115° |
| L3 (primary) | both ≤ 4.0 Å | ≤ 3.5 Å | 95–120° |
| L4 | both ≤ 4.0 Å | ≤ 4.0 Å | 90–125° |
| L5 | one ≤ 4.0 Å and both ≤ 5.0 Å | ≤ 4.5 Å | 90–130° |

The primary analysis uses L3. The script also reports the L1–L5 sensitivity results.
