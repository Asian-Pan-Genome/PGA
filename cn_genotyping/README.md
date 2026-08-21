# *PGA* copy-number genotyping

Depth-based, assembly-calibrated inference of total and paralog-specific *PGA* copy number from short-read data.

## Pseudo-reference and informative sites

The paralog-resolved variable-site set is generated in [`paralog_evolution/`](../paralog_evolution/). Build a pseudo-*PGA* reference from an aligned *PGA34A* template:

```bash
python build_pseudo_reference.py \
    --template PGA34A.template.aligned.fa \
    --vcf PGA.variable_sites.vcf \
    --out-prefix PGA

python make_site_manifest.py \
    --vcf PGA.pseudo.vcf \
    --out PGA.manifest.tsv
```

The informative-site manifest uses the publication-level classes *PGA34A*, *PGA34B*, and *PGA5* and retains sites with a predominant within-class allele frequency ≥ 0.70 by default.

## Read extraction and depth normalization

For pre-aligned whole-genome BAM/CRAM data:

```bash
python extract_pga_reads.py \
    --alignment sample.cram \
    --reference GRCh38.fa \
    --pga-bed GRCh38.PGA.bed \
    --control-bed GRCh38.control.bed \
    --out-prefix sample \
    --threads 16
```

The GRCh38 analyses used the following local diploid-depth controls:

```text
chr11	59501328	60532588
chr11	61506000	62093042
```

Extracted reads are remapped to the pseudo-*PGA* reference before feature calculation. Raw-read cohorts, including HGDP and ancient genomes, use an initial *PGA*-enrichment mapping step; their baseline depth is estimated from total sequence yield divided by 3.1 Gbp. `set_sequence_yield_baseline.py` implements this normalization. Pre-aligned single-end ancient genomes can be processed with `extract_pga_reads_single_end.py`.

## CN inference

Prepare a sample table containing the pseudo-reference alignment and depth file for each sample, then extract continuous CN features:

```bash
python pga_cn_gmm.py extract \
    --sample-list samples.tsv \
    --reference PGA.pseudo.fa \
    --manifest PGA.manifest.tsv \
    --out cohort.features.tsv \
    --threads 16
```

Fit constrained one-dimensional GMMs with assembly-resolved truth labels:

```bash
python pga_cn_gmm.py fit \
    --features train.features.tsv \
    --truth assembly_truth.tsv \
    --model-out PGA_CN.model.pkl
```

Target-cohort features can be supplied through `--unlabeled-features` during fitting. Predict integer CN states with:

```bash
python pga_cn_gmm.py predict \
    --features target.features.tsv \
    --model PGA_CN.model.pkl \
    --out target.predicted.tsv
```

Separate models are fitted for total *PGA*, total *PGA34*, *PGA5*, and *PGA34B*. *PGA34A* CN is calculated as total *PGA34* CN minus *PGA34B* CN.

## Benchmarking

```bash
python benchmark_cn.py \
    --predictions target.predicted.tsv \
    --truth assembly_truth.tsv \
    --out benchmark.tsv

python plot_cn_benchmark.py \
    --predictions target.predicted.tsv \
    --truth assembly_truth.tsv \
    --out-prefix copies_true_vs_pred
```

Benchmarking reports exact integer-copy accuracy and RMSE against assembly-resolved CN labels.
