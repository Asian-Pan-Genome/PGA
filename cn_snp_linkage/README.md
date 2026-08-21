# SNP associations with *PGA34A* copy number

Genome-wide association of common 1KGP SNPs with diploid *PGA34A* copy number. EUR is the primary analysis population to match the ancestry background of the UK Biobank analyses; the same framework is applied separately to the remaining superpopulations for comparison.

## Population panels and association testing

Prepare population-specific MAF-filtered panels, phenotypes, and the first ten within-population PCs:

```bash
bash prepare_1kg_population_panels.sh \
    1kGP_high_coverage.autosome.SNV.vcf.gz \
    sample_metadata.tsv \
    population_panels \
    16
```

The workflow retains biallelic SNVs with MAF ≥ 0.05, performs LD pruning with `--indep-pairwise 50 5 0.2`, and calculates PC1–PC10 with PLINK v1.9.

Test SNP dosage against *PGA34A* CN within each superpopulation:

```bash
bash run_cn_snp_association.sh \
    population_panels \
    association_results \
    16
```

Models use *PGA34A* CN as the quantitative response and SNP dosage as the predictor, adjusting for PC1–PC10. Genome-wide significant associations can be combined with:

```bash
python merge_association_results.py \
    --input-dir association_results \
    --p-threshold 5e-8 \
    --output all_populations_aggregated.P5e-8.tsv
```

For the manuscript analysis, EUR CN-linked variants within the focal region were further annotated with the GWAS Catalog and GTEx cis-eQTL resources.

1KGP high-coverage phased variants: <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/>
