# UK Biobank associations

Association of inferred *PGA34A* copy number with molecular traits and PheCode phenotypes in UK Biobank.

## Association cohort

Starting from *PGA* CN calls generated in [`cn_genotyping/`](../cn_genotyping/), prepare the analysis cohort with:

```bash
python prepare_ukb_association.py \
    --cn UKB.PGA.predicted.tsv \
    --phenotypes UKB.phenotypes.tsv \
    --molecular-fields molecular_fields.txt \
    --out-prefix UKB.PGA34A
```

The script applies the genetic QC, ancestry, sex-concordance, kinship, and missingness filters used in the manuscript and requires non-missing *PGA34A* CN, age, sex, and PC1–PC10. The final study cohort contained 328,897 participants.

## Molecular-trait associations

```bash
python molecular_association.py \
    --cohort UKB.PGA34A.cohort.tsv \
    --proteomics UKB.proteomics.tsv \
    --molecular-fields molecular_fields.txt \
    --threads 32 \
    --out PGA34A.molecular_association.tsv
```

Traits with fewer than 100 non-missing measurements are excluded. Each retained trait is rank-based inverse-normal transformed and tested by OLS with *PGA34A* CN as the predictor, adjusting for age, sex, and PC1–PC10. *P* values are Bonferroni-corrected across tested traits.

## PheWAS

```bash
Rscript PheWAS.R \
    UKB.PGA34A.cohort.tsv \
    UKB.PGA34A.icd10.tsv \
    PGA34A.PheWAS.csv \
    10
```

UKB ICD-10 codes are converted to standard decimal formatting before PheCode mapping. A participant is classified as a case when at least one diagnosis maps to the corresponding PheCode; PheCode exclusions and sex-specific restrictions define eligible controls. Logistic regression tests *PGA34A* CN as a continuous predictor with age, sex, and PC1–PC10 as covariates. Phenotypes with fewer than 20 cases or controls are excluded, and *P* values are Bonferroni-corrected across tested phenotypes.
