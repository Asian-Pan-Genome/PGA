# Ancient DNA

This directory contains the temporal analyses of *PGA* copy number and local SNP-panel diversity in ancient Eurasian genomes.

Ancient *PGA* copy-number inference uses the same short-read genotyping framework as the modern NGS cohorts. The genotyping workflow itself is maintained in [`cn_genotyping/`](../cn_genotyping/) and is not duplicated here.

## Requirements

- Python 3
- NumPy
- pandas
- SciPy
- matplotlib
- seaborn
- statsmodels
- [EIGENSOFT](https://github.com/DReichLab/EIG) (`convertf`)

## 1. Temporal analysis of ancient *PGA* copy number

After ancient genomes have been genotyped with [`cn_genotyping/`](../cn_genotyping/), prepare a tab-delimited table containing:

```text
Sample
Region
Date_Mean
Pred_PGA_Total
```

`Region` should contain the labels `East Asia` and/or `West Eurasia`. The script retains the original ancient-CN table convention in which `Date_Mean` is negative for samples in the past (for example, `-9000` for 9 kya).

Run:

```bash
python ancient_copies.PGA_sensitivity.py \
    --input ancient_PGA_CN.tsv \
    --out-dir ancient_PGA_CN_temporal \
    --bin-sizes 1000 2000 3000 \
    --oldest-bin-start 9000
```

The script performs a continuous-age ordinary least-squares regression separately for East and West Eurasia and 1-, 2-, and 3-kyr bin-size sensitivity analyses. Samples older than 9 kya are pooled into the terminal time bin by default.

Main outputs are two multi-page PDFs, one per region, and `ancient_copies.PGA_sensitivity.fit_summary.tsv` containing the slope, two-sided slope-test *P* value and R² for the continuous and binned analyses.

## 2. AADR SNP-panel analysis

The local-diversity analysis uses the AADR v66 2M release in hg19 coordinates. The required source files are:

```text
v66.2M.aadr.PUB.anno
v66.2M.aadr.PUB.geno
v66.2M.aadr.PUB.ind
v66.2M.aadr.PUB.snp
```

### Prepare a chromosome 11 ancient subset

`aadr_prepare_subset.py` retains ancient individuals and, when multiple Genetic IDs represent the same Individual ID, keeps the record with the largest number of SNPs hit on the AADR 2M autosomal targets. The original AADR `.ind` order is preserved so genotype columns remain aligned with the source `.geno` file.

Run:

```bash
python aadr_prepare_subset.py \
    --prefix v66.2M.aadr.PUB \
    --anno v66.2M.aadr.PUB.anno \
    --chrom 11 \
    --out-prefix aadr_v66_2M.hg19.chr11
```

This writes the selected metadata, sample list, a full-order `.ind` file with non-selected samples marked as `Ignore`, and a `convertf` parameter file. The conversion intentionally retains the entire chromosome 11 because the downstream analysis uses chromosome-wide matched-region backgrounds.

Convert the selected samples to chromosome 11 EIGENSTRAT format:

```bash
convertf -p aadr_v66_2M.hg19.chr11.convertf.par
```

The resulting files are:

```text
aadr_v66_2M.hg19.chr11.eigenstrat.geno
aadr_v66_2M.hg19.chr11.eigenstrat.snp
aadr_v66_2M.hg19.chr11.eigenstrat.ind
aadr_v66_2M.hg19.chr11.metadata.tsv
```

### Temporal diversity analysis around *PGA*

Run:

```bash
python aadr_pga_windowed_selection.py \
    --geno aadr_v66_2M.hg19.chr11.eigenstrat.geno \
    --snp aadr_v66_2M.hg19.chr11.eigenstrat.snp \
    --ind aadr_v66_2M.hg19.chr11.eigenstrat.ind \
    --metadata aadr_v66_2M.hg19.chr11.metadata.tsv \
    --anno v66.2M.aadr.PUB.anno \
    --outdir aadr_PGA_selection \
    --bin-sizes 1000 2000 3000 \
    --n-boot 1000 \
    --trend-permutations 9999 \
    --threads 8
```

The analysis uses an 80°E boundary for the East–West comparison. East Eurasian samples are analysed across the DR with the duplicated *PGA* gene cluster masked, whereas West Eurasian samples are analysed across DR-L. Diversity is calculated in 10-kb windows with a 1-kb step and summarized using callable-site-weighted windowed π and mean windowed Tajima's D-like statistics. Temporal summaries use a 1-kyr rolling window and 1-, 2-, and 3-kyr bins, with individual-level bootstrap confidence intervals and matched-length chromosome 11 backgrounds. All coordinates are hg19.

Because AADR genotypes are predominantly pseudo-haploid and the 2M panel is SNP-ascertained, these statistics are interpreted as relative SNP-panel diversity measures rather than conventional whole-sequence estimates.
