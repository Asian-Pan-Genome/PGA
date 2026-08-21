# Functional analyses and diet association

Computational analyses supporting the functional section: protein structural modelling and molecular dynamics, and the regional association between *PGA* copy number and plant-derived protein fraction.

## Protein structure and molecular dynamics

[`protein_structure_md/`](protein_structure_md/) contains construction of the PGA34A/PGA34B/PGA5 protein–substrate models, five independent 100-ns molecular-dynamics replicates per isoform, and analysis of active-site water geometry.

## Association between *PGA* CN and plant-derived protein fraction

Plant-derived protein fraction was calculated from FAOSTAT Food Balances as the proportion of total protein supply derived from vegetal products, averaged over 1961–1980.

```bash
python prepare_faostat.py \
    --input <FAOSTAT_Food_Balances.csv> \
    --output fao_baseline_1961_1980.csv
```

The 1KGP association analysis uses inferred diploid *PGA* CN, population metadata, regional mean PC1–PC3 from [`population_structure/`](../population_structure/), and the FAOSTAT table:

```bash
python pga_diet_association.py \
    --cnv <1KGP_PGA_CN.tsv> \
    --metadata <1KGP_metadata.tsv> \
    --pca <joint_diploid_PCA.tsv> \
    --fao fao_baseline_1961_1980.csv \
    --out-prefix pga_diet
```

The script reports the regional Spearman correlation with a two-sided permutation *P* value and the PC1–PC3-adjusted partial Spearman correlation with a two-sided Freedman–Lane permutation *P* value. Both use 100,000 permutations by default.

FAOSTAT food-supply data are treated as regional dietary proxies rather than individual dietary measurements.
