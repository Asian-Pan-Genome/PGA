# Paralog evolution

This directory contains the sequence-based analyses used to classify human *PGA* paralogs, summarize paralog-variable sites, and identify coding variants with predicted functional effects.

## Requirements

- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [trimAl](https://github.com/inab/trimal)
- [IQ-TREE](https://iqtree.github.io/)
- [iTOL](https://itol.embl.de/)
- Python 3
- [Biopython](https://biopython.org/)
- [SeqKit](https://bioinf.shenwei.me/seqkit/)
- [minimap2](https://github.com/lh3/minimap2), including `paftools.js`
- `bgzip`
- [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), including `filter_vep`

## 1. Phylogenetic classification of *PGA* copies

Two full-length sequence sets are used:

- all *PGA* copies, to separate *PGA5* from *PGA34*;
- *PGA34* copies, to resolve the *PGA34A1*, *PGA34A2* and *PGA34B* lineages.

For the all-*PGA* tree:

```bash
mafft --auto --thread -1 all.PGA.fa > all.PGA.mafft.fa
trimal -in all.PGA.mafft.fa -out all.PGA.mafft.trim.fa -automated1 -keepheader
iqtree -s all.PGA.mafft.trim.fa --prefix PGA -T auto -B 1000 -bnni -safe
```

For the *PGA34* tree:

```bash
mafft --auto --thread -1 all.PGA34.fa > all.PGA34.mafft.fa
trimal -in all.PGA34.mafft.fa -out all.PGA34.mafft.trim.fa -automated1 -keepheader
iqtree -s all.PGA34.mafft.trim.fa --prefix PGA34 -T auto -B 1000 -bnni -safe
```

The resulting trees were visualized in iTOL and lineage assignments were made manually from the phylogenetic clustering. *PGA34A1* and *PGA34A2* are retained as sublineages for sequence-divergence analyses but are combined as *PGA34A* for downstream copy-number analyses.

## 2. Variable and lineage-diagnostic sites

`find_var_site_in_msa.py` summarizes variable columns from a classified FASTA alignment. Sequence IDs must begin with the lineage label followed by `#`, for example:

```text
PGA34A1#sample1.hap1.PGA
PGA34A2#sample2.hap2.PGA
PGA34B#sample3.hap1.PGA
```

For each variable position, the script reports the predominant REF/ALT state (`GT`) and its within-lineage frequency (`GF`) for each requested group.

For example, to compare the three *PGA34* lineages:

```bash
python find_var_site_in_msa.py \
    --alignment all.PGA34.mafft.trim.classified.fa \
    --reference <reference_sequence_id> \
    --groups PGA34A1 PGA34A2 PGA34B \
    --output PGA34.variable_sites.vcf
```

The script itself does not impose a fixed-site threshold. Lineage-diagnostic sites were identified manually from the output by inspecting positions with opposing predominant allele states and high within-lineage genotype frequencies.

For the variable-site VCF used by [`cn_genotyping/`](../cn_genotyping/), relabel *PGA34A1* and *PGA34A2* copies with the shared `PGA34A#` prefix and run:

```bash
python find_var_site_in_msa.py \
    --alignment all.PGA.classified.fa \
    --reference <PGA34A_reference_sequence_id> \
    --groups PGA34A PGA34B PGA5 \
    --output PGA.variable_sites.vcf
```

This keeps the publication-level paralog nomenclature consistent between `paralog_evolution/` and `cn_genotyping/`.

## 3. High-impact variants

Variants were called for individual *PGA* copies relative to the corresponding paralog reference sequence and annotated with VEP. Prepare a text file containing one copy ID per line and a FASTA containing those copies, then run:

```bash
bash annotate_pga_variants.sh \
    <copy_list> \
    <copies_fasta> \
    <reference_fasta> \
    <reference_gff> \
    <output_dir> \
    16
```

The script extracts each copy with SeqKit, aligns it to the selected reference with minimap2, calls variants with `paftools.js`, annotates them with VEP, and filters out `MODIFIER` annotations using:

```text
IMPACT != MODIFIER
```

In the analysis reported in the manuscript, the variants retained after this filtering step were all annotated as `HIGH` impact.
