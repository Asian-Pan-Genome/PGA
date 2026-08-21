# Functional and diet association analyses

Protein-structure modelling, molecular-dynamics simulation, and population-level association between *PGA* copy number and dietary protein composition.

## Analysis structure

| Path | Procedure |
| --- | --- |
| [`protein_structure_md/`](protein_structure_md/) | Construct PGA34A, PGA34B, and PGA5 protein–LSFMAIPP models, run five independent 100-ns molecular-dynamics replicates per isoform, and quantify catalysis-ready water geometries. |
| [`prepare_faostat.py`](prepare_faostat.py) | Calculate country-level plant-derived protein fractions from FAOSTAT Food Balances for 1961–1980. |
| [`pga_diet_association.py`](pga_diet_association.py) | Merge 1KGP *PGA* CN, population metadata, PCA, and FAOSTAT records, and test the association between regional *PGA* CN and plant-derived protein fraction. |

## Protein structure and molecular dynamics

Commands, model-construction steps, molecular-dynamics parameters, and catalysis-ready water criteria are documented in [`protein_structure_md/`](protein_structure_md/).

## Association between *PGA* CN and plant-derived protein fraction

Dietary protein composition was obtained from [FAOSTAT Food Balances](https://www.fao.org/faostat/en/#data/FBSH). The analysis used the old-methodology Food Balances records from 1961–1980.

Prepare the country-level plant-derived protein fraction table with:

```bash
python prepare_faostat.py \
    --input <FAOSTAT_Food_Balances.csv> \
    --output fao_baseline_1961_1980.csv
```

The script calculates the fraction of total protein supply contributed by plant-derived products and averages it across 1961–1980.

The *PGA* CN–diet association uses the 1KGP cohort. Prepare the inferred diploid *PGA* CN table, 1KGP population metadata, the diploid PCA table from [`population_structure/`](../population_structure/), and the FAOSTAT table generated above. Run:

```bash
python pga_diet_association.py \
    --cnv <1KGP_PGA_CN.tsv> \
    --metadata <1KGP_metadata.tsv> \
    --pca <joint_diploid_PCA.tsv> \
    --fao fao_baseline_1961_1980.csv \
    --out-prefix pga_diet
```

The script retains 1KGP samples with defined study superpopulation labels, maps populations to the corresponding FAOSTAT regions, and calculates regional mean diploid *PGA* CN and regional mean PC1–PC3. It reports the unadjusted Spearman correlation with a two-sided permutation *P* value and the PC1–PC3-adjusted partial Spearman correlation with a two-sided Freedman–Lane permutation *P* value. Both tests use 100,000 permutations by default.

Main outputs are:

```text
pga_diet.sample_match.tsv
pga_diet.region_data.tsv
pga_diet.statistics.tsv
pga_diet.raw.pdf
pga_diet.pc_adjusted.pdf
```
