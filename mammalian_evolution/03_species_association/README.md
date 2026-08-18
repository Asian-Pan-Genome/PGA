# Stage 03: species-level association

This stage selects one assembly per species, integrates copy number, ecological traits and assembly quality, and fits phylogenetic generalized least-squares (PGLS) models on the VertLife tree.

## Data preparation

Representative assemblies are ranked by anchor completeness, contig/scaffold N50 and complete ancestral-gene count. Deposited tables record all duplicate-species decisions.

Diet proportions and body mass are based on EltonTraits 1.0. Missing species were curated from congeneric records or Animal Diversity Web and retained as authoritative inputs. Detailed diets were grouped as carnivore, omnivore, insectivore or plant-dominant. Body mass was analysed on the `log10` scale.

Copy-number and trait records were intersected with the deposited 295-species VertLife subtree.

## Scripts

| Script | Role |
| --- | --- |
| `scripts/01_prepare_species_association_input.py` | Normalize copy counts, select representative assemblies and join assembly, diet and body-mass metadata. |
| `scripts/02_add_diet_group_proportions.py` | Derive plant/animal food proportions and the higher-order diet group variables used by PGLS. |
| `scripts/03_run_pgls_body_mass.R` | Fit body-mass PGLS models with Pagel's lambda estimated by maximum likelihood. |
| `scripts/04_run_pgls_diet_models.R` | Fit higher-order diet PGLS models and export coefficients, AIC and used-species tables. |
| `scripts/05_run_pgls_models.sh` | Run both deposited PGLS analyses using repository-relative defaults. |

Run both R analyses from the repository root:

```bash
bash 03_species_association/scripts/05_run_pgls_models.sh
```

Custom inputs can be passed as:

```bash
bash 03_species_association/scripts/05_run_pgls_models.sh \
  bodymass_input.tsv diet_input.tsv tree.nwk output_dir
```

## Models

Copy number is modelled as `log1p(CN)`. Body mass and contig N50 use the `log10` scale. PGLS uses a Brownian correlation structure with Pagel's lambda estimated by maximum likelihood in `ape` and `nlme`. Output tables record species retained after tree intersection and complete-case filtering.

Files in `results/` provide the coefficients, model summaries, AIC values and species sets supporting Fig. 6a and related supplementary analyses. No final plotting script is deposited.

## Deposited resources and outputs

- `../shared_resources/ecology/`: 295-species EltonTraits/body-mass data and authoritative diet-classification inputs;
- `../shared_resources/species_metadata/`: representative-assembly decisions, order labels and quality metadata;
- `../shared_resources/phylogeny/295_sp.tree`: the analysis subtree;
- `results/295_*`: prepared copy-number, ecology and representative-assembly tables;
- `results/pga_*_pgls_*`: PGLS coefficients, AIC tables, model summaries and used-species records.

The `*_copy.tsv` files are provenance snapshots from the author analysis and are not exact content duplicates.

## Dependencies

Python 3 with `numpy` and `pandas`; R with `ape` and `nlme`; and Bash for the convenience wrapper.
