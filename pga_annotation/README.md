# *PGA* annotation

Human *PGA* copies were annotated with Liftoff from GRCh38.p14/GENCODE v47; ape annotations combine TOGA2 projections with Liftoff rescue.

## Human genomes

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

`<target.PGA.fa>` is the *PGA* locus extracted from a haplotype-resolved assembly.

## Ape genomes

Ape annotations use the T2T primate assemblies and hg38-based TOGA2 annotations. Candidate TOGA2 *PGA* models of 8–12 Kbp were retained, the corresponding local locus was masked, and Liftoff was used to recover additional copies. The two annotation sets were then merged.

Prepare a tab-delimited manifest:

```text
sample_id	species	fasta	toga_gtf
mPanTro3_hap1	Pan_troglodytes	/path/to/mPanTro3.hap1.fa	/path/to/geneAnnotation.gtf.gz
mGorGor1_hap1	Gorilla_gorilla	/path/to/mGorGor1.hap1.fa	/path/to/geneAnnotation.gtf.gz
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

The merged annotation for each assembly is written to:

```text
<output_dir>/<sample_id>/<sample_id>.TOGA_liftoff.gff
```

T2T primate assemblies: <https://github.com/marbl/Primates>  
TOGA2 annotations: <https://genome.senckenberg.de/download/TOGA2/>
