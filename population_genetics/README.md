# Population genetics

This directory contains the population-genetic analyses used to examine genetic diversity, population differentiation, and selection around the human *PGA* locus. Analyses were performed using variants from the local assembly-based pangenome graph and the phased 1KGP high-coverage panel.

## Requirements

- [PGGB](https://github.com/pangenome/pggb)
- [vg](https://github.com/vgteam/vg)
- [vcfbub](https://github.com/pangenome/vcfbub)
- [BCFtools/HTSlib](https://github.com/samtools/bcftools)
- [vcfwave](https://github.com/ekg/vcflib)
- [VCFtools](https://vcftools.github.io/)
- [VCF-kit](https://github.com/AndersenLab/VCF-kit)
- [selscan](https://github.com/szpiech/selscan)
- [predictGMAP](https://github.com/szpiech/predictGMAP)
- [BetaScan](https://github.com/ksiewert/BetaScan)
- [Relate](https://myersgroup.github.io/relate/)
- Python 3 with `msprime`, `numpy`, `pandas`, and `matplotlib`

## 1. Construct the local PGGB graph and decompose variants

Prepare a FASTA containing the sequences spanning the extended *PGA* region. GRCh38 was included as the reference path in our analysis.

```bash
N_SEQS=$(grep -c '^>' input.extend.fa)

pggb \
    -i input.extend.fa \
    -o PGA_extend_PGGB \
    -n "${N_SEQS}" \
    -c 2 \
    -t 32 \
    --skip-viz
```

Decompose graph bubbles against the GRCh38 path and remove nested bubbles:

```bash
vg deconstruct \
    -P GRCh38 \
    -H '?' \
    -a \
    -t 32 \
    <graph.smooth.final.gfa> | \
    bgzip -@ 32 -c > graph.raw.vcf.gz

vcfbub \
    -l 0 \
    -r 100000 \
    -i graph.raw.vcf.gz \
    > graph.unnested.vcf

bcftools norm \
    -m -any \
    graph.unnested.vcf | \
    bcftools view \
        -c 1 \
        -Ov \
        -o graph.biallelic.vcf
```

Subsequent variant-ID annotation, `vcfwave` decomposition, normalization, and duplicate merging followed the same procedure described in [`population_structure/`](../population_structure#2-construct-regional-graphs-and-decompose-variants).

## 2. Prepare population-genetic variant panels

For the assembly-derived panel, analyses were restricted to individuals represented by two complete haplotypes. Small variants were retained with:

```bash
bcftools view \
    -i 'STRLEN(REF) < 50 && STRLEN(ALT) < 50' \
    graph.final.vcf.gz \
    -Ov \
    -o PGGB.SMV.vcf
```

Variants overlapping tandem-repeat regions or the *PGA* gene cluster were excluded. Allele counts and frequencies were recalculated after sample or population subsetting, and invariant sites introduced by subsetting were removed.

The population-specific assembly-panel VCFs are referred to below as:

```text
PGGB.SMV.noTR.noPGA.dip_only.<POP>.vcf.gz
```

For the 1KGP analyses, unrelated individuals and biallelic SNPs from the phased high-coverage chromosome 11 panel were retained and separated by population or superpopulation as required.

## 3. SNP density

SNP density was calculated in 10-kb sliding windows with a 1-kb step across the *PGA* structurally variable region using both the assembly-derived and 1KGP panels. The main-figure analysis used SNPs with MAF >= 0.05 within the corresponding panel or population.

## 4. Nucleotide diversity

Windowed nucleotide diversity was calculated after retaining variants with MAF >= 0.05:

```bash
populations=(
    ALL EAS AFR ARB SAS AMR SEA EUR CSA CAS WAS
    "AFR-E&S" AFR-NE AFR-W
)

for pop in "${populations[@]}"; do
    vcftools \
        --gzvcf "PGGB.SMV.noTR.noPGA.dip_only.${pop}.vcf.gz" \
        --maf 0.05 \
        --window-pi 10000 \
        --window-pi-step 1000 \
        --out "PGGB.SMV.noTR.noPGA.dip_only.${pop}"
done
```

The same analysis was applied to the corresponding 1KGP population panels.

## 5. Tajima's D

Tajima's D was calculated without an additional MAF filter, using 10-kb windows with a 1-kb step:

```bash
for pop in "${populations[@]}"; do
    bcftools view \
        -Ov \
        -o "PGGB.SMV.noTR.noPGA.dip_only.${pop}.vcf" \
        "PGGB.SMV.noTR.noPGA.dip_only.${pop}.vcf.gz"

    vk tajima \
        10000 \
        1000 \
        "PGGB.SMV.noTR.noPGA.dip_only.${pop}.vcf" \
        > "PGGB.SMV.noTR.noPGA.dip_only.${pop}.Tajima.D"
done
```

The same analysis was applied to the corresponding 1KGP population panels.

## 6. Population differentiation

Windowed Weir-Cockerham Fst was calculated from the 1KGP panel between EAS and the other major superpopulations after retaining variants with MAF >= 0.05:

```bash
bcftools view \
    -q 0.05:minor \
    -Oz \
    -o chr11.maf5.vcf.gz \
    chr11.filtered.shapeit2-duohmm-phased.SNP.vcf.gz

for pop in AFR AMR CSA EUR SAS; do
    vcftools \
        --gzvcf chr11.maf5.vcf.gz \
        --weir-fst-pop EAS.pop \
        --weir-fst-pop "${pop}.pop" \
        --fst-window-size 10000 \
        --fst-window-step 1000 \
        --out "EAS_${pop}"
done
```

## 7. XP-EHH and XP-nSL

XP-EHH and XP-nSL were calculated from the MAF >= 0.05 1KGP SNP panel.

For XP-EHH, interpolate missing chromosome 11 genetic-map positions with `predictGMAP`:

```bash
bcftools query \
    -f '%POS\n' \
    chr11.filtered.shapeit2-duohmm-phased.SNP.bi.vcf.gz \
    > query.positions

awk 'NR > 1 {
    OFS="\t"
    print "11", "chr11:"$2, $4, $2
}' hg38.chr11.full.genetic_map.txt \
    > chr11.reference.map

predictGMAP \
    --max-gap 5000000 \
    --query query.positions \
    --ref chr11.reference.map \
    --out chr11.predict.map
```

Prepare population-specific MAF-filtered VCFs:

```bash
for pop in EAS AFR AMR CSA EUR SAS; do
    bcftools view \
        -q 0.05:minor \
        -Ov \
        -o "${pop}.maf5.vcf" \
        "${pop}.vcf"
done
```

XP-EHH:

```bash
for pop in AFR AMR CSA EUR SAS; do
    selscan \
        --xpehh \
        --vcf EAS.maf5.vcf \
        --vcf-ref "${pop}.maf5.vcf" \
        --wagh \
        --map chr11.predict.map \
        --threads 8 \
        --out "EAS_${pop}"

    norm \
        --xpehh \
        --files "EAS_${pop}.xpehh.out"
done
```

XP-nSL:

```bash
for pop in AFR AMR CSA EUR SAS; do
    selscan \
        --xpnsl \
        --vcf EAS.maf5.vcf \
        --vcf-ref "${pop}.maf5.vcf" \
        --threads 8 \
        --out "EAS_${pop}"

    norm \
        --xpnsl \
        --files "EAS_${pop}.xpnsl.out"
done
```

## 8. Folded Beta1

Folded Beta1 scores were calculated without an additional MAF >= 0.05 prefilter:

```bash
populations=(
    EAS AFR ARB SAS AMR SEA EUR CSA CAS WAS
    "AFR-E&S" AFR-NE AFR-W
)

for pop in "${populations[@]}"; do
    python BetaScan.py \
        -i <(
            bcftools query \
                -f '%POS %INFO/AC %INFO/AN %INFO/AF\n' \
                "PGGB.SMV.noTR.noPGA.dip_only.${pop}.vcf.gz" | \
            awk '{
                if ($4 >= 0.5)
                    print $1, $3-$2, $3
                else
                    print $1, $2, $3
            }'
        ) \
        -fold \
        -m 0.15 \
        -p 20 \
        -std \
        -theta 0.001 \
        -o "PGGB.SMV.noTR.noPGA.dip_only.${pop}.beta1_folded"
done
```

The same analysis was applied to the corresponding 1KGP population panels.

## 9. Neutral simulation

The demographic null for the EAS diversity depletion is implemented in [`EAS_DR_weightedPI_neutral_scan.py`](EAS_DR_weightedPI_neutral_scan.py). The script simulates 223 diploid EAS individuals under the Relate-inferred demographic history across a 906-kb search region. Variants with MAF >= 0.05 are retained, and each replicate is scanned with a DR-length interval in 1-kb steps. The minimum variant-count-weighted mean nucleotide diversity is compared with the observed EAS DR statistic.

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

## 10. CN-linked SNPs

Analyses of SNPs and haplotypes linked to *PGA* copy number are described separately in [`cn_snp_linkage/`](../cn_snp_linkage/).

## 11. Local genealogies

Local genealogies for selected CN-linked variants were reconstructed from the phased 1KGP chromosome 11 panel with Relate.

```bash
RelateFileFormats \
    --mode ConvertFromVcf \
    --haps 1KG.ALL.haps \
    --sample 1KG.ALL.sample \
    -i 1KG.ALL

PrepareInputFiles.sh \
    --haps 1KG.ALL.haps \
    --sample 1KG.ALL.sample \
    --ancestor homo_sapiens_ancestor_11.fa \
    --mask StrictMask_chr11.fa.gz \
    -o 1KG.ALL.prepared

RelateParallel.sh \
    --threads 32 \
    --haps 1KG.ALL.prepared.haps.gz \
    --sample 1KG.ALL.prepared.sample.gz \
    --dist 1KG.ALL.prepared.dist.gz \
    -m 1.25e-8 \
    -N 20000 \
    --map hg38.chr11.full.genetic_map.txt \
    -o 1KG.ALL.prepared.relate
```

Extract the genealogy at a selected site:

```bash
TreeView.sh \
    --haps 1KG.ALL.prepared.haps.gz \
    --sample 1KG.ALL.prepared.sample.gz \
    --anc 1KG.ALL.prepared.relate.anc \
    --mut 1KG.ALL.prepared.relate.mut \
    --poplabels 1KG.ALL.poplabels.superpop \
    --bp_of_interest <position> \
    --years_per_gen 28 \
    -o 1KG.ALL.prepared.relate.<position>
```
