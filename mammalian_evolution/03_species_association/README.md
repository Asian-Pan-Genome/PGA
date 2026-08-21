# Species-level phylogenetic association

Selection of one representative genome per species and phylogeny-informed tests of extant *PGA* copy number against mammalian ecological traits.

Representative assemblies are prioritized by local anchor completeness and assembly contiguity. Species with resolved *PGA* CN and ecological metadata are intersected with the VertLife phylogeny, yielding the 295-species analysis panel.

## Workflow

| Script | Role |
| --- | --- |
| `scripts/01_prepare_species_association_input.py` | Select representative assemblies and integrate CN, ecology, and assembly metadata. |
| `scripts/02_add_diet_group_proportions.py` | Derive dietary proportions and higher-order diet groups. |
| `scripts/03_run_pgls_body_mass.R` | Fit body-mass PGLS models. |
| `scripts/04_run_pgls_diet_models.R` | Fit dietary PGLS models. |
| `scripts/05_run_pgls_models.sh` | Run the deposited PGLS analyses. |

```bash
bash 03_species_association/scripts/05_run_pgls_models.sh
```

*PGA* CN is analysed as `log1p(CN)`. Body mass and contig N50 are log10-transformed where included as covariates. PGLS models estimate Pagel's λ by maximum likelihood.
