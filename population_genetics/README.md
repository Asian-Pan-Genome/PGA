# Population genetics

Population-genetic analyses of sequence variation surrounding the human *PGA* locus using the assembly-derived pangenome panel and the phased 1KGP high-coverage panel.

## Local pangenome variants

A PGGB graph was constructed across the *PGA* structurally variable region from 1,348 complete haplotypes, with GRCh38 retained as the reference path:

```bash
N_SEQS=$(grep -c '^>' input.extend.fa)
pggb -i input.extend.fa -o PGA_extend_PGGB -n "${N_SEQS}" -c 2 -t 32 --skip-viz
```

Variants were decomposed with `vg deconstruct -a`, filtered with `vcfbub`, and normalized with the same `vcfwave`/collapse-bubble procedure used in [`population_structure/`](../population_structure/). Assembly-panel analyses retained individuals with two complete haplotypes and excluded variants overlapping GIAB tandem repeats or the *PGA* gene cluster where specified.

## Diversity and selection statistics

| Analysis | Panel | Implementation |
| --- | --- | --- |
| SNP density | Assembly and 1KGP | 10-Kbp windows, 1-Kbp step |
| Nucleotide diversity (π) | Assembly and 1KGP | VCFtools; MAF ≥ 0.05; 10-Kbp windows, 1-Kbp step |
| Tajima's D | Assembly and 1KGP | VCF-kit; 10-Kbp windows, 1-Kbp step |
| Fst | 1KGP | VCFtools; EAS versus other superpopulations |
| XP-EHH / XP-nSL | 1KGP | selscan; chromosome 11 genetic map |
| Folded Beta1 | Assembly and 1KGP | BetaScan (`-m 0.15 -p 20`) |

For XP-EHH, missing genetic-map positions were interpolated with `predictGMAP` before running selscan. Chromosome-wide distributions were used as empirical backgrounds for regional comparisons.

## Neutral simulation

`EAS_DR_weightedPI_neutral_scan.py` tests whether the EAS diversity depletion can arise under the Relate-inferred EAS demographic history. Each of 10,000 replicates simulates 223 diploid individuals across a 906-Kbp region, retains variants with MAF ≥ 0.05, and records the minimum variant-count-weighted mean π across DR-sized windows.

```bash
python EAS_DR_weightedPI_neutral_scan.py \
    --pi-file PGGB.SMV.noTR.noPGA.dip_only.EAS.vcf.gz.windowed.pi \
    --demography-csv pop_size.Relate.csv \
    --target-region EAS \
    --region chr11:61090000-61506000 \
    --search-length-bp 906000 \
    --n-diploid 223 \
    --maf 0.05 \
    --replicates 10000 \
    --out-prefix EAS.DR.weightedMeanPI.maf5.scanSVR
```

## CN-linked SNPs and local genealogies

Genome-wide SNP association with *PGA34A* CN is implemented in [`cn_snp_linkage/`](../cn_snp_linkage/). Local genealogies at selected CN-linked SNPs were inferred with Relate from phased 1KGP chromosome 11 variants using a mutation rate of `1.25e-8`, initial `N = 20000`, the GRCh38 genetic map, and 28 years per generation for visualization.
