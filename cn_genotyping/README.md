# *PGA* copy-number genotyping

This directory contains the scripts used to infer total and paralog-specific *PGA* copy numbers from short-read NGS data. The same core genotyping framework was applied to multiple modern and ancient cohorts, with cohort-specific handling of the input reads and depth baseline.

## Requirements

- [BWA](https://github.com/lh3/bwa)
- [SAMtools](https://github.com/samtools/samtools)
- Python 3
- [pysam](https://github.com/pysam-developers/pysam)
- [Biopython](https://biopython.org/)
- pandas
- NumPy
- scikit-learn
- joblib
- matplotlib and seaborn (benchmark plots only)

## 1. Build the pseudo-*PGA* reference

The paralog-resolved variable-site VCF is generated in [`paralog_evolution/`](../paralog_evolution/). Starting from this VCF and an aligned representative *PGA34A* sequence, build the pseudo-reference with:

```bash
python build_pseudo_reference.py \
    --template PGA34A.template.aligned.fa \
    --vcf PGA.variable_sites.vcf \
    --out-prefix PGA
```

This produces:

```text
PGA.pseudo.fa
PGA.pseudo.vcf
```

The pseudo-reference uses the *PGA34A* sequence as the template. Alignment-gap positions in the template are replaced by the corresponding non-gap allele in the variable-site VCF, with the VCF alleles recoded to the pseudo-reference coordinates.

Generate the paralog-informative site manifest from the pseudo-reference VCF:

```bash
python make_site_manifest.py \
    --vcf PGA.pseudo.vcf \
    --out PGA.manifest.tsv
```

The VCF used here should use the current paralog labels *PGA34A*, *PGA34B* and *PGA5*. Sites with a within-paralog predominant-allele frequency of at least 0.70 are retained by default.

Index the pseudo-reference before read mapping:

```bash
bwa index PGA.pseudo.fa
samtools faidx PGA.pseudo.fa
```

## 2. Prepare reads and depth features

### Pre-aligned modern genomes

For BAM/CRAM datasets with whole-genome alignments, extract reads overlapping the *PGA* genes and calculate the depth features with:

```bash
python extract_pga_reads.py \
    --alignment sample.cram \
    --reference GRCh38.fa \
    --pga-bed GRCh38.PGA.bed \
    --control-bed GRCh38.control.bed \
    --out-prefix sample \
    --threads 16
```

The *PGA* BED file contains the reference *PGA3*, *PGA4* and *PGA5* intervals in the fourth column. For the GRCh38-based analyses in this study, the control BED was:

```text
chr11	59501328	60532588
chr11	61506000	62093042
```

APG alignments were duplicate-filtered before this step. 1KGP and UK Biobank alignments were processed directly with the corresponding reference-specific *PGA* and control intervals.

The extracted read pairs are remapped to the pseudo-*PGA* reference, for example:

```bash
bwa mem -t 16 PGA.pseudo.fa sample.R1.fq.gz sample.R2.fq.gz | \
    samtools sort -@ 16 -o sample.pseudo.bam
samtools index sample.pseudo.bam
```

### Raw-read cohorts

HGDP paired-end reads and ancient single-end reads were first mapped to the pseudo-*PGA* reference to retain candidate *PGA* reads. The candidate reads were then remapped to the appropriate whole-genome reference, followed by final *PGA* read extraction and pseudo-reference remapping.

Because this enrichment step does not retain the flanking control regions, the depth baseline for these samples was estimated from total input sequence yield divided by 3.1 Gbp. `set_sequence_yield_baseline.py` can update the `Baseline_depth` field either from the BWA log of the initial mapping:

```bash
python set_sequence_yield_baseline.py \
    --depth sample.depth.tsv \
    --bwa-log sample.initial_pseudo.bwa.log \
    --out sample.depth.normalized.tsv
```

or from a metadata table containing `run_accession` and `base_count`:

```bash
python set_sequence_yield_baseline.py \
    --depth sample.depth.tsv \
    --metadata ancient_samples.tsv \
    --sample SAMPLE \
    --out sample.depth.normalized.tsv
```

Ancient single-end reads were mapped with BWA aln/samse using:

```text
-n 0.01 -l 16500 -o 2
```

For ancient samples available only as pre-aligned single-end BAM/CRAM files, use:

```bash
python extract_pga_reads_single_end.py \
    --alignment sample.bam \
    --reference reference.fa \
    --pga-bed reference.PGA.bed \
    --control-bed reference.control.bed \
    --out-prefix sample
```

Archaic genomes were handled from their pre-aligned files using reference-matched *PGA* and control intervals; the T2T-CN1-aligned samples in this study were remapped to the pseudo-reference with BWA aln/samse.

## 3. Extract continuous CN features

Prepare a tab-delimited sample list containing the pseudo-reference alignment and depth file for each sample:

```text
Sample	Alignment	Depth
sample1	/path/to/sample1.pseudo.bam	/path/to/sample1.depth.tsv
sample2	/path/to/sample2.pseudo.bam	/path/to/sample2.depth.tsv
```

Then run:

```bash
python pga_cn_gmm.py extract \
    --sample-list samples.tsv \
    --reference PGA.pseudo.fa \
    --manifest PGA.manifest.tsv \
    --out cohort.features.tsv \
    --threads 16
```

The resulting table contains continuous depth-derived features for total *PGA*, total *PGA34* and *PGA5*, together with the informative-site feature used for *PGA34B*.

## 4. Fit the CN models and genotype samples

Assembly-resolved CN labels are used to calibrate the constrained one-dimensional GMMs. The truth table should contain `Sample`, `PGA34A`, `PGA34B` and `PGA5`. It may contain either one diploid row per sample or two haplotype rows per sample.

Fit the models:

```bash
python pga_cn_gmm.py fit \
    --features train.features.tsv \
    --truth assembly_truth.tsv \
    --model-out PGA_CN.model.pkl
```

For cross-cohort calibration, the target-cohort feature distribution can also be supplied during fitting:

```bash
python pga_cn_gmm.py fit \
    --features train.features.tsv \
    --truth assembly_truth.tsv \
    --unlabeled-features target.features.tsv \
    --model-out PGA_CN.model.pkl
```

Predict integer CN states with:

```bash
python pga_cn_gmm.py predict \
    --features target.features.tsv \
    --model PGA_CN.model.pkl \
    --out target.predicted.tsv
```

Separate models are fitted for total *PGA*, total *PGA34*, *PGA5* and *PGA34B*. *PGA34A* CN is calculated as total *PGA34* CN minus *PGA34B* CN.

## 5. Benchmarking

Compare predicted CNs with assembly-resolved truth labels using exact integer-copy accuracy and RMSE:

```bash
python benchmark_cn.py \
    --predictions target.predicted.tsv \
    --truth assembly_truth.tsv \
    --out benchmark.tsv
```

True-versus-predicted CN heatmaps can be generated with:

```bash
python plot_cn_benchmark.py \
    --predictions target.predicted.tsv \
    --truth assembly_truth.tsv \
    --out-prefix copies_true_vs_pred
```
