# Genome-wide copy-number differentiation

This directory contains the analysis used to compare *PGA* copy-number differentiation with genome-wide tandemly duplicated gene families.

## Requirements

- [OrthoFinder](https://github.com/davidemms/OrthoFinder)
- [AGAT](https://github.com/NBISweden/AGAT)
- [SeqKit](https://github.com/shenwei356/seqkit)
- [Liftoff](https://github.com/agshumate/Liftoff)
- Python 3 with `pandas` and `numpy`

## 1. Define tandemly duplicated gene families

We used Ensembl-canonical protein-coding genes from GENCODE v47 on chromosomes 1-22 and X. Protein sequences were extracted from GRCh38.p14 and renamed with their GENCODE gene symbols.

```bash
rg 'gene_type=protein_coding' gencode.v47.annotation.gff3 \
  | awk '$3=="transcript"' \
  | rg -v 'Ensembl_canonical' \
  | cut -f 9 \
  | awk -F '[=;]' '{print $2}' \
  > noncanonical_transcripts.txt

agat_sp_filter_feature_from_kill_list.pl \
  --gff gencode.v47.annotation.gff3 \
  --kill_list noncanonical_transcripts.txt \
  -o gencode.v47.Ensembl_canonical.gff3

rg '^chr([1-9]|1[0-9]|2[0-2]|X)\\t' gencode.v47.Ensembl_canonical.gff3 \
  > gencode.v47.Ensembl_canonical.1-22X.gff3

agat_sp_extract_sequences.pl \
  -g gencode.v47.Ensembl_canonical.1-22X.gff3 \
  -f GRCh38.p14.genome.fa \
  -t cds -p --keep_attributes \
  -o GRCh38.Ensembl_canonical.transcript.pep.1-22X.fa
```

Rename FASTA headers to gene symbols:

```bash
rg '>' GRCh38.Ensembl_canonical.transcript.pep.1-22X.fa \
  | rg -o 'gene_name=[^ ]+' \
  | cut -d '=' -f 2 \
  | paste <(rg '>' GRCh38.Ensembl_canonical.transcript.pep.1-22X.fa | sed 's/^>//') - \
  > GRCh38.Ensembl_canonical.transcript.pep.1-22X.fa.kv.txt

seqkit replace \
  -p '(.*)' \
  -r '{kv}' \
  -k GRCh38.Ensembl_canonical.transcript.pep.1-22X.fa.kv.txt \
  GRCh38.Ensembl_canonical.transcript.pep.1-22X.fa \
  > GRCh38.Ensembl_canonical.transcript.pep.1-22X.renamed.fa

mv GRCh38.Ensembl_canonical.transcript.pep.1-22X.renamed.fa \
   GRCh38.Ensembl_canonical.transcript.pep.1-22X.fa
```

We supplied two identical copies of this protein set to OrthoFinder to group duplicated genes into orthogroups:

```bash
cp GRCh38.Ensembl_canonical.transcript.pep.1-22X.fa \
   GRCh38.Ensembl_canonical.transcript.pep.1-22X.cp.fa

orthofinder -f orthofinder_input -t 16 -a 16
```

The OrthoFinder orthogroups were then separated by chromosome and genomic proximity to define local tandemly duplicated gene families. Genes in the same orthogroup were grouped while successive members remained within 5 Mb; unassigned genes were retained as singleton families.

```bash
python define_local_gene_families.py \
  --orthogroups OrthoFinder/Results/Orthogroups/Orthogroups.tsv \
  --unassigned OrthoFinder/Results/Orthogroups/Orthogroups_UnassignedGenes.tsv \
  --gff3 gencode.v47.Ensembl_canonical.1-22X.gff3 \
  --reference-column GRCh38.Ensembl_canonical.transcript.pep.1-22X \
  --duplicate-column GRCh38.Ensembl_canonical.transcript.pep.1-22X.cp \
  --max-extension 5000000 \
  --output local_gene_families.tsv
```

This produced 18,575 gene families in the analysis used for the manuscript.

## 2. Build the genome-wide gene-family CN matrix

Gene copy numbers were obtained from Liftoff annotations of each haplotype assembly. For each haplotype, prepare a two-column table containing gene symbol and copy number:

```text
gene    CN
PGA3    1
PGA4    2
PGA5    1
```

Prepare a manifest listing the per-haplotype CN tables:

```text
sample_hap    gene_cn_file
sample1.hap1  /path/to/sample1.hap1.gene_cn.tsv
sample1.hap2  /path/to/sample1.hap2.gene_cn.tsv
```

Then sum the CNs of all member genes within each family:

```bash
python build_gene_family_cn_matrix.py \
  --families local_gene_families.tsv \
  --manifest gene_cn_manifest.tsv \
  --output all_samples.merge.counts
```

## 3. Calculate genome-wide Vst ranks

For each population pair, `calculate_vst.py` calculates Vst and the absolute difference in mean CN for every retained gene family, and ranks families by:

```text
Vst x |Delta CN|
```

The script uses the full exclusion pattern from the original analysis to remove selected repetitive or highly copy-variable gene families before ranking. After filtering, 17,289 gene families remained.

### Total *PGA*

Use the *PGA* family CN already present in the genome-wide matrix:

```bash
python calculate_vst.py \
  --cn-matrix all_samples.merge.counts \
  --sample-table PGA.copies.tsv \
  --population-table assemblies.new_superpop.list \
  --target-gene PGA3 \
  --target-label PGA \
  --full-pair EAS EUR \
  --processes 8 \
  --out-prefix PGA
```

### *PGA34A*

For the paralog-specific analysis, replace the *PGA* family row with haplotype-level *PGA34A* CN (`PGA34A1 + PGA34A2`) and use the same genome-wide background:

```bash
python calculate_vst.py \
  --cn-matrix all_samples.merge.counts \
  --sample-table PGA34A1_A2_B.tsv \
  --population-table assemblies.new_superpop.list \
  --target-gene PGA3 \
  --target-label PGA34A \
  --target-cn-table PGA34A1_A2_B.tsv \
  --target-cn-columns PGA34A1 PGA34A2 \
  --full-pair EAS EUR \
  --processes 8 \
  --out-prefix PGA34A
```

The main outputs are:

```text
<prefix>.pairwise_rank.tsv
<prefix>.EAS_vs_EUR.genomewide.tsv
```
