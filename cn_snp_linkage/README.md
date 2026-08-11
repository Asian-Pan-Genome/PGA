# SNP associations with *PGA34A* copy number

This directory contains the 1KGP analysis used to identify common SNPs associated with *PGA34A* copy number.

## Requirements

- [BCFtools](https://github.com/samtools/bcftools)
- [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
- Python 3 with `pandas`

## Input data

We used the phased 1KGP high-coverage autosomal SNV panel:

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/

A sample metadata table is also required, with sample ID, diploid *PGA34A* copy number, and analysis superpopulation in the first three columns. For example:

```text
Sample	PGA34A_CN	New_superpop
HG00323	0	EUR
HG00268	1	EUR
...
```

In our analysis, AFR combines the AFR-W and AFR-E&S groups; AMR, EAS, EUR, SAS, and CSA are treated separately.

Variant IDs in the input VCF are expected in `CHROM:POS:REF:ALT` format for downstream result merging.

## Prepare population panels, phenotypes, and principal components

Run:

```bash
bash prepare_1kg_population_panels.sh \
    1kGP_high_coverage.autosome.SNV.vcf.gz \
    sample_metadata.tsv \
    population_panels \
    16
```

For each superpopulation, the script generates the sample list and PLINK phenotype file, recalculates allele frequencies after sample subsetting, retains variants with MAF >= 0.05, performs LD pruning with `--indep-pairwise 50 5 0.2`, and calculates the first ten genetic PCs.

## Run SNP-*PGA34A* CN associations

Run:

```bash
bash run_cn_snp_association.sh \
    population_panels \
    association_results \
    16
```

Genome-wide SNP associations with *PGA34A* CN are tested separately in each superpopulation using PLINK linear regression, with PC1-PC10 included as covariates.

## Merge significant associations

Genome-wide significant associations can be combined across superpopulations with:

```bash
python merge_association_results.py \
    --input-dir association_results \
    --p-threshold 5e-8 \
    --output all_populations_aggregated.P5e-8.tsv
```

The output contains each significant variant, the superpopulation(s) in which it was detected, and the corresponding effect sizes and P values.

## Functional annotation

For the manuscript analysis, CN-linked variants in the focal region were further evaluated using GWAS Catalog and GTEx cis-eQTL annotations. Trait-level interpretation and the final selection of reported GWAS associations were manually curated and are not implemented here as an automated workflow.
