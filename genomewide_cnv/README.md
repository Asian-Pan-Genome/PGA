# Genome-wide copy-number differentiation

Genome-wide comparison of *PGA* copy-number differentiation with local protein-coding gene families.

## Gene-family definition

Ensembl-canonical protein-coding genes from GENCODE v47 were clustered with OrthoFinder. Orthogroups were then split by chromosome and genomic proximity to define local families; successive members separated by more than 5 Mbp were assigned to different families, and unassigned genes were retained as singletons.

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

This procedure generated 18,575 local gene families before exclusion of highly repetitive gene classes.

## Gene-family CN matrix

For each haplotype, provide a table of gene symbol and annotated copy number, then sum member-gene CN within each family:

```bash
python build_gene_family_cn_matrix.py \
    --families local_gene_families.tsv \
    --manifest gene_cn_manifest.tsv \
    --output all_samples.merge.counts
```

## Population differentiation

`calculate_vst.py` jointly calculates Vst, the absolute difference in mean CN, and the `Vst × |ΔCN|` ranking for each population contrast. After exclusion of repetitive gene classes, 17,289 families were retained in the manuscript analysis.

Total *PGA*:

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

For the paralog-specific analysis, the target family row is replaced by haplotype-level *PGA34A* CN while retaining the same genome-wide background:

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
