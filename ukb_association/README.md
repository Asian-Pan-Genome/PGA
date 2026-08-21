# UK Biobank association analyses

This directory contains the UK Biobank molecular-trait association and PheWAS analyses for *PGA34A* copy number. *PGA* copy-number genotyping is handled separately in [`cn_genotyping`](../cn_genotyping/); the scripts here start from the resulting `Pred_PGA34A` calls.

## Requirements

For molecular association:

- Python 3
- pandas
- NumPy
- SciPy
- statsmodels

For PheWAS:

- R
- [PheWAS](https://github.com/PheWAS/PheWAS)
- dplyr

## 1. Prepare the UKB association cohort

Prepare the *PGA* CN prediction table from [`cn_genotyping`](../cn_genotyping/), the baseline UKB phenotype table, and `molecular_fields.txt`.

`molecular_fields.txt` contains the 314 baseline non-proteomic fields considered in this study. The minimum sample-size filter is applied later by the association script; 313 fields passed this filter in the study data.

Run:

```bash
python prepare_ukb_association.py \
    --cn UKB.PGA.predicted.tsv \
    --phenotypes UKB.phenotypes.tsv \
    --molecular-fields molecular_fields.txt \
    --out-prefix UKB.PGA34A
```

The script applies the UKB genetic QC and ancestry filters used in the study and requires non-missing *PGA34A* CN, age, sex, and PC1-PC10. With the study inputs, the final association cohort contained 328,897 participants.

Main outputs:

```text
UKB.PGA34A.cohort.tsv
UKB.PGA34A.icd10.tsv
```

The first file contains the association cohort, covariates, and baseline molecular fields. The second contains the hospital inpatient ICD-10 diagnoses from UKB Data-Field 41270 in long format for PheWAS.

## 2. Molecular-trait association

Plasma protein measurements are supplied separately as a table with `eid` plus one column per protein. Run:

```bash
python molecular_association.py \
    --cohort UKB.PGA34A.cohort.tsv \
    --proteomics UKB.proteomics.tsv \
    --molecular-fields molecular_fields.txt \
    --threads 32 \
    --out PGA34A.molecular_association.tsv
```

Traits with fewer than 100 non-missing measurements are skipped. Each retained trait is rank-based inverse-normal transformed and tested by OLS with *PGA34A* CN as the predictor, adjusting for age, sex, and PC1-PC10. P values are Bonferroni-corrected across all tested traits.

## 3. PheWAS

Run the PheWAS from the prepared cohort and ICD-10 table:

```bash
Rscript PheWAS.R \
    UKB.PGA34A.cohort.tsv \
    UKB.PGA34A.icd10.tsv \
    PGA34A.PheWAS.csv \
    10
```

`PheWAS.R` restores standard ICD-10 decimal formatting before mapping UKB diagnosis codes to PheCodes. A participant is treated as a case when at least one ICD-10 code maps to the corresponding PheCode. PheCode exclusions and sex-specific restrictions are applied, while eligible participants without inpatient ICD-10 records remain available as controls where appropriate.

Associations are tested by logistic regression with *PGA34A* CN as a continuous predictor, adjusting for age, sex, and PC1-PC10. Phenotypes with fewer than 20 cases or controls are skipped, and P values are Bonferroni-corrected across the successfully tested PheCode phenotypes.
