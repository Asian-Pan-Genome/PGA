# *PGA* annotation

This directory contains the commands and scripts used to annotate *PGA* genes in human haplotype-resolved assemblies and ape genomes.

## Requirements

- [Liftoff](https://github.com/agshumate/Liftoff)
- [AGAT](https://github.com/NBISweden/AGAT)
- [BEDTools](https://github.com/arq5x/bedtools2)

## Human *PGA* annotation

Human *PGA* copies were annotated directly with Liftoff.

In our analysis, we used the GRCh38.p14 *PGA* locus as the reference, together with a Liftoff annotation database built from the corresponding GENCODE v47 annotation.

```bash
liftoff \
    <target.PGA.fa> \
    <reference.PGA.fa> \
    -sc 0.95 \
    -copies \
    -db <reference.PGA.gff_db> \
    -polish \
    -exclude_partial \
    -p 16 \
    -o <target>.liftoff.gff
```

Here, `<target.PGA.fa>` is the *PGA* locus extracted from a haplotype-resolved assembly.

## Ape *PGA* annotation

Ape *PGA* annotation combines existing TOGA projections with Liftoff-based rescue of additional copies.

The T2T ape genome assemblies used in this study are available from:

https://github.com/marbl/Primates

The corresponding TOGA annotations can be downloaded from:

https://genome.senckenberg.de/download/TOGA2/

Before running the script, prepare a tab-delimited manifest containing the assembly FASTA and TOGA annotation for each genome. For example:

```text
sample_id	species	fasta	toga_gtf
mPanTro3_hap1	Pan_troglodytes	/path/to/mPanTro3.hap1.fa	/path/to/geneAnnotation.gtf.gz
mGorGor1_hap1	Gorilla_gorilla	/path/to/mGorGor1.hap1.fa	/path/to/geneAnnotation.gtf.gz
mSymSyn1_hap1	Symphalangus_syndactylus	/path/to/mSymSyn1.hap1.fa	/path/to/geneAnnotation.gtf.gz
```

Run:

```bash
bash annotate_ape_pga.sh \
    ape_manifest.tsv \
    <reference.PGA.fa> \
    <reference.PGA.gff_db> \
    output \
    16
```

The arguments are:

```text
1. manifest.tsv
2. reference PGA FASTA
3. Liftoff annotation database
4. output directory
5. number of threads (default: 16)
```

For each genome, the script first retains TOGA-annotated *PGA* genes with gene lengths of 8–12 kb. It then extracts the local *PGA* region with 100 kb of flanking sequence on each side, masks the retained TOGA annotations, and runs Liftoff to recover additional *PGA* copies. The TOGA and Liftoff annotations are finally merged into:

```text
<output_dir>/<sample_id>/<sample_id>.TOGA_liftoff.gff
```
