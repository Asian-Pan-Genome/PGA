# Ancient DNA

Temporal analyses of *PGA* copy number and local SNP-panel diversity in ancient Eurasian genomes.

## *PGA* copy-number trajectories

Ancient *PGA* CN is inferred with the short-read framework in [`cn_genotyping/`](../cn_genotyping/). From the resulting CN table, run:

```bash
python ancient_copies.PGA_sensitivity.py \
    --input ancient_PGA_CN.tsv \
    --out-dir ancient_PGA_CN_temporal \
    --bin-sizes 1000 2000 3000 \
    --oldest-bin-start 9000
```

The script tests continuous-age trends separately in East and West Eurasia and evaluates alternative temporal bin widths. The analysis uses an 80°E boundary for the East–West comparison.

## AADR regional diversity

The sequence-diversity analysis uses the AADR v66 2M panel in hg19 coordinates. Prepare the chromosome 11 ancient subset with:

```bash
python aadr_prepare_subset.py \
    --prefix v66.2M.aadr.PUB \
    --anno v66.2M.aadr.PUB.anno \
    --chrom 11 \
    --out-prefix aadr_v66_2M.hg19.chr11

convertf -p aadr_v66_2M.hg19.chr11.convertf.par
```

Run the temporal analysis with:

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

East Eurasian analyses focus on the desert region (DR), with the duplicated *PGA* gene cluster masked; West Eurasian analyses focus on DR-L. Statistics are calculated in 10-Kbp windows with a 1-Kbp step and summarized through rolling and fixed temporal bins with individual-level bootstrap uncertainty and chromosome 11 matched-region backgrounds.

Because AADR genotypes are SNP-ascertained and predominantly pseudo-haploid, nucleotide diversity and Tajima's D-like values are interpreted as relative SNP-panel measures rather than conventional whole-sequence estimates.
