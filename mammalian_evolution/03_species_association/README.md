# Species-level association

Commands and input tables used to select one assembly per species, combine *PGA* copy number with ecological and assembly-quality covariates, and fit phylogenetic generalized least-squares models.

## Inputs

- species-level copy-number tables generated from the classified copy units;
- [`../shared_resources/phylogeny/295_sp.tree`](../shared_resources/phylogeny/295_sp.tree);
- diet and body-mass tables in [`../shared_resources/ecology/`](../shared_resources/ecology/);
- representative-assembly and quality metadata in [`../shared_resources/species_metadata/`](../shared_resources/species_metadata/).

## 1. Prepare the species-level dataset

Representative assemblies are ranked by anchor completeness, contig/scaffold N50, and complete ancestral-gene count.

Prepare species-level copy-number, diet, body-mass, and assembly metadata with:

```bash
python scripts/01_prepare_species_association_input.py \
    --help
```

Diet proportions and higher-order diet groups are generated with:

```bash
python scripts/02_add_diet_group_proportions.py \
    --help
```

Diet records are grouped as carnivore, omnivore, insectivore, or plant-dominant for the higher-order analyses. Body mass is analysed on the `log10` scale. Copy-number and trait records are intersected with the 295-species VertLife subtree before model fitting.

## 2. PGLS analyses

| Script | Procedure |
| --- | --- |
| `scripts/03_run_pgls_body_mass.R` | Fit body-mass PGLS models. |
| `scripts/04_run_pgls_diet_models.R` | Fit higher-order diet PGLS models and export coefficients, AIC, and used-species tables. |
| `scripts/05_run_pgls_models.sh` | Run both PGLS analyses with repository-relative defaults. |

Run both analyses from `mammalian_evolution/` with:

```bash
bash 03_species_association/scripts/05_run_pgls_models.sh
```

Custom input files can be supplied as:

```bash
bash 03_species_association/scripts/05_run_pgls_models.sh \
    bodymass_input.tsv \
    diet_input.tsv \
    tree.nwk \
    output_dir
```

Copy number is modelled as `log1p(CN)`. Body mass and contig N50 are `log10` transformed where included. PGLS uses a Brownian correlation structure with Pagel's lambda estimated by maximum likelihood using `ape` and `nlme`.

## Outputs

The `results/` directory contains:

- prepared species-level copy-number and ecological-trait tables;
- representative-assembly tables;
- PGLS coefficient tables;
- model summaries and AIC tables;
- the species retained in each fitted model.
