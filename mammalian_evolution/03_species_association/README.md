# Stage 03: species-level association

This stage prepares one representative assembly per species, joins copy number to diet, body mass and assembly-quality metadata, intersects the data with the VertLife tree, and fits phylogenetic generalized least-squares (PGLS) models.

## Data preparation

Representative assemblies are selected using, in order, anchor completeness, contig/scaffold N50 and the number of complete ancestral genes. The released tables retain the selected assembly and duplicate-species decision records.

Diet proportions and body mass are based on EltonTraits 1.0. Missing species were curated from congeneric information or Animal Diversity Web and are retained only as explicit authoritative inputs. The detailed diet labels are also grouped into broader carnivore, omnivore, insectivore and plant-dominant categories. Body mass is analysed on the `log10` scale.

The analysis intersects copy-number and trait records with the deposited 295-species VertLife subtree.

## Scripts

| Script | Role |
| --- | --- |
| `scripts/01_prepare_species_association_input.py` | Normalize updated copy-count columns, select representative assemblies and join assembly, diet and body-mass metadata. |
| `scripts/02_add_diet_group_proportions.py` | Derive plant/animal food proportions and the higher-order diet group variables used by PGLS. |
| `scripts/03_run_pgls_body_mass.R` | Fit body-mass PGLS models with Pagel's lambda estimated by maximum likelihood. |
| `scripts/04_run_pgls_diet_models.R` | Fit the higher-order diet PGLS models and export model coefficients, AIC and used-species tables. |
| `scripts/05_run_pgls_models.sh` | Run both deposited PGLS analyses using repository-relative defaults. |

The two R analyses can be reproduced together from the repository root with:

```bash
bash 03_species_association/scripts/05_run_pgls_models.sh
```

Custom inputs can be passed as:

```bash
bash 03_species_association/scripts/05_run_pgls_models.sh \
  bodymass_input.tsv diet_input.tsv tree.nwk output_dir
```

## Models

Copy number is modelled as `log1p(CN)`. Body mass and contig N50 are modelled on the `log10` scale. PGLS models use a Brownian correlation structure with Pagel's lambda estimated by maximum likelihood through `ape` and `nlme`; the deposited outputs record the exact species retained after tree/data intersection and complete-case filtering.

The final coefficient, model-summary, AIC and used-species files in `results/` are the tabular/statistical products supporting Fig. 6a and the associated supplementary analyses. No final figure-generation script is deposited.

## Deposited resources and outputs

- `../shared_resources/ecology/`: the 295-species EltonTraits/body-mass subset and authoritative diet-classification inputs;
- `../shared_resources/species_metadata/`: representative-assembly decisions, order labels and quality metadata;
- `../shared_resources/phylogeny/295_sp.tree`: the analysis subtree;
- `results/295_*`: prepared copy-number, ecology and representative-assembly tables;
- `results/pga_*_pgls_*`: final PGLS coefficients, AIC tables, model summaries and used-species records.

The `*_copy.tsv` files are retained provenance snapshots used during the author analysis; they are not exact duplicate files according to content hashes.

## Dependencies

Python 3 with `numpy` and `pandas`; R with `ape` and `nlme`; and Bash for the convenience wrapper.
