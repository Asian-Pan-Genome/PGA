# Population structure

This directory contains the workflow used for joint diploid PCA of the assembly cohort and the 1KGP high-coverage panel.

## Requirements

- [Minigraph-Cactus](https://github.com/ComparativeGenomicsToolkit/cactus)
- [prepare-vcf-MC](https://github.com/eblerjana/genotyping-pipelines/tree/main/prepare-vcf-MC)
- [collapse-bubble](https://github.com/Han-Cao/collapse-bubble)
- [vcfwave](https://github.com/ekg/vcflib)
- [minimap2](https://github.com/lh3/minimap2)
- [BCFtools](https://github.com/samtools/bcftools)
- [SAMtools](https://github.com/samtools/samtools)
- [BEDTools](https://github.com/arq5x/bedtools2)
- [SeqKit](https://github.com/shenwei356/seqkit)
- [PLINK](https://www.cog-genomics.org/plink/)
- [Snakemake](https://snakemake.github.io/)
- Python 3 with `pandas`, `matplotlib`, `numpy`, and `pysam`

## Input

We used 20 fixed 1-Mb regions on chromosome 11 selected from the GRCh38 GIAB easy-region stratification. Prepare these regions as a BED file and extract the corresponding GRCh38 sequences into a FASTA file with headers `GRCh38.region1` to `GRCh38.region20`.

Assembly input is provided as a tab-delimited manifest. For example:

```text
sample	hap	source	fasta
HG005	hap1	HPRC	/path/to/HG005.hap1.fa
HG005	hap2	HPRC	/path/to/HG005.hap2.fa
CHM13v2	hap0	REF	/path/to/CHM13v2.fa
CN1v1	hap0	REF	/path/to/CN1v1.fa
```

The phased 1KGP high-coverage VCF is available from:

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/

## 1. Extract assembly regions

Project the 20 GRCh38 regions onto each assembly and extract the corresponding sequences:

```bash
bash extract_assembly_regions.sh \
    assemblies.tsv \
    GRCh38.regions.fa \
    assembly_regions \
    8
```

Regions that cannot be projected for a haplotype are recorded in `region*/failed.list`.

## 2. Construct regional graphs and decompose variants

For each region, construct a Minigraph-Cactus graph with GRCh38 as the VCF reference and CHM13v2 and CN1v1 as additional references. Variant-processing steps in this workflow borrow scripts and workflow components from [prepare-vcf-MC](https://github.com/eblerjana/genotyping-pipelines/tree/main/prepare-vcf-MC) and [collapse-bubble](https://github.com/Han-Cao/collapse-bubble), together with `vcfwave`.

Set the paths to these external workflows and the sample-sex table at the beginning of `run_minigraph_cactus.sh`, then run:

```bash
bash run_minigraph_cactus.sh \
    assemblies.tsv \
    GRCh38.regions.fa \
    assembly_regions \
    graph \
    32
```

The final assembly-derived SNP file for each region is:

```text
graph/regionN/regionN.SNP.vcf.gz
```

These VCFs retain biallelic SNPs with assembly-panel MAF >= 0.05.

## 3. Joint assembly-1KGP PCA

Run the joint PCA workflow using the same 20 regions from the 1KGP chr11 VCF:

```bash
bash run_joint_pca.sh \
    assemblies.tsv \
    hg38.regions.bed \
    GRCh38.regions.fa \
    graph \
    1KGP.chr11.vcf.gz \
    PGA.copies.tsv \
    assemblies.new_superpop.list \
    1KGP.metadata.txt \
    joint_pca \
    32
```

The workflow standardizes the assembly-derived regional VCFs to a common diploid sample set, converts the corresponding 1KGP variants to the same local coordinates, and retains exact shared SNPs by `CHROM`, `POS`, `REF`, and `ALT`. Sites showing strong panel-specific differences are then removed using individuals represented in both datasets.

The concordance filter uses:

```text
minimum callable overlapping pairs: 10
minimum genotype concordance:        0.98
maximum panel missingness:           0.20
maximum missingness difference:      0.05
```

The filtered joint dataset is further restricted to MAF >= 0.05 and genotype missingness <= 0.05, followed by LD pruning with `--indep-pairwise 50 5 0.2`. The first 10 PCs are then calculated with PLINK.

## 4. Plot PCA

```bash
python plot_joint_pca.py \
    --eigenvec joint_pca/plink/joint.PCA.eigenvec \
    --metadata joint_pca/joint_pca.metadata.tsv \
    --out joint_pca/joint.PCA
```
