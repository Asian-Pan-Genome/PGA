# Paralog evolution

Phylogenetic classification of human *PGA* copies, paralog-divergent sites, predicted high-impact variants, and sequence context of the *PGA34* 904-bp deletion.

## Phylogenetic classification

Full-length human *PGA* sequences were analysed with ape homologs as outgroups. The first tree separates *PGA5* and *PGA34*; a second tree resolves *PGA34A1*, *PGA34A2*, and *PGA34B*.

```bash
mafft --auto --thread -1 all.PGA.fa > all.PGA.mafft.fa
trimal -in all.PGA.mafft.fa -out all.PGA.mafft.trim.fa -automated1 -keepheader
iqtree -s all.PGA.mafft.trim.fa --prefix PGA -T auto -B 1000 -bnni -safe

mafft --auto --thread -1 all.PGA34.fa > all.PGA34.mafft.fa
trimal -in all.PGA34.mafft.fa -out all.PGA34.mafft.trim.fa -automated1 -keepheader
iqtree -s all.PGA34.mafft.trim.fa --prefix PGA34 -T auto -B 1000 -bnni -safe
```

*PGA34A1* and *PGA34A2* are retained for sequence-divergence analyses and combined as *PGA34A* for downstream CN analyses.

## Paralog-divergent sites

`find_var_site_in_msa.py` summarizes variable alignment columns after paralog assignment. Sequence identifiers should begin with the paralog label followed by `#`.

```bash
python find_var_site_in_msa.py \
    --alignment all.PGA.classified.fa \
    --reference <PGA34A_reference_sequence_id> \
    --groups PGA34A PGA34B PGA5 \
    --output PGA.variable_sites.vcf
```

This paralog-resolved variable-site set is also used to construct the pseudo-*PGA* reference in [`cn_genotyping/`](../cn_genotyping/).

## Predicted high-impact variants

```bash
bash annotate_pga_variants.sh \
    <copy_list> \
    <copies_fasta> \
    <reference_fasta> \
    <reference_gff> \
    <output_dir> \
    16
```

Each copy is aligned to the corresponding paralog reference with minimap2/`paftools.js` and annotated with Ensembl VEP. The retained variants in the manuscript analysis were all annotated as `HIGH` impact.

The *PGA34* 904-bp deletion was examined separately in a local Minigraph–Cactus graph containing human and ape *PGA* sequences. The deletion junction and flanking sequences were extracted from the graph and compared for local sequence features, including microhomology.
