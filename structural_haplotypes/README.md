# Structural haplotypes

Code used to resolve structural haplotypes (SHs) of the human *PGA* locus and to examine recurrent NAHR.

## Requirements

- [PGR-TK](https://github.com/Sema4-Research/pgr-tk)
- [ODGI](https://github.com/pangenome/odgi)
- Python 3 with `numpy`, `pandas`, `scipy`, `ete3`, `scikit-learn`, `hdbscan` and `matplotlib`
- R with `tidyverse`, `ggtree`, `gggenes`, `aplot`, `RColorBrewer` and `ape`
- `bcftools`, `tabix`, `bgzip` and `samtools`
- [vcflib](https://github.com/vcflib/vcflib), including `vcfcreatemulti`
- MAFFT
- trimAl
- IQ-TREE 2
- minimap2
- [TreeCluster](https://github.com/niemasd/TreeCluster) for the flanking-tree clustering benchmark

## 1. Principal bundles and SH assignment

Run PGR-TK on the 1,348 gap-free *PGA* haplotypes:

```bash
pgr-pbundle-decomp \
    -r 4 \
    --min-span 12 \
    --bundle-length-cutoff 10 \
    --min-branch-size 16 \
    PGA_region.000.nogap.fa \
    PGA_48_56_4_12_10_16.
```

The raw PGR-TK BED contains a number of short or fragmented bundle calls. Merge them into the principal-bundle representation with:

```bash
python curate_principal_bundles.py \
    PGA_48_56_4_12_10_16.bed \
    PGA.principal_bundles.bed
```

Then combine the bundle paths with the paralog annotation table:

```bash
python assign_structural_haplotypes.py \
    PGA.principal_bundles.bed \
    PGA34A_B.tsv \
    PGA.structural_haplotypes.tsv
```

`PGA34A_B.tsv` should contain `Sample`, `Hap`, `Source` and `PGAs`. The final analysis contains 36 SHs. In the SH labels, A, B and C denote *PGA34A*, *PGA34B* and *PGA5*, respectively; X and Y denote *PGA34/5* and *PGA5/34* hybrid genes.

## 2. PGGB Jaccard-distance tree

The SH plot is ordered by UPGMA clustering of pairwise path Jaccard distances from the PGGB graph. Generate the pairwise distances with `odgi similarity -d`, then run:

```bash
python pggb_jaccard_upgma.py \
    --distance pggb.jaccard.distance.tsv \
    --haplotypes PGA.structural_haplotypes.tsv \
    --output PGA.jaccard.upgma.newick
```

The ODGI table should contain `group.a`, `group.b` and `jaccard.distance`.

## 3. Gene track and SH plot

Generate the *PGA* gene track with:

```bash
python generate_gene_track.py \
    PGA_region.000.hit \
    id.list \
    PGA34A_B.tsv \
    PGA.gene_track.bed
```

`id.list` should contain `Sample`, `Hap` and `GFF File`.

After selecting one representative haplotype for each SH, subset the principal-bundle and gene-track BED files to those representatives and plot them with:

```bash
Rscript plot_structural_haplotypes.R \
    representative.principal_bundles.bed \
    representative.gene_track.bed \
    PGA.jaccard.upgma.newick \
    PGA.structural_haplotypes.pdf
```

Representative selection is only used for visualization and is not part of SH assignment.

## 4. Upstream and downstream flanking trees

The upstream tree (UT) uses GRCh38 `chr11:61191000-61203514`, and the downstream tree (DT) uses `chr11:61251444-61263958`.

`build_flanking_phylogeny.sh` takes a graph-derived VCF, GRCh38 FASTA, one interval and the orthologous chimpanzee sequence. The chimpanzee FASTA should use `chimpanzee` as the sequence ID.

For example:

```bash
bash build_flanking_phylogeny.sh \
    graph.vcf.gz \
    GRCh38.fa \
    chr11:61191000-61203514 \
    chimpanzee.upstream.fa \
    PGA.upstream \
    16

bash build_flanking_phylogeny.sh \
    graph.vcf.gz \
    GRCh38.fa \
    chr11:61251444-61263958 \
    chimpanzee.downstream.fa \
    PGA.downstream \
    16
```

The script keeps variants up to 1 kb, merges overlapping indel records with `vcfcreatemulti`, reconstructs one consensus sequence per haplotype, and then runs MAFFT, trimAl and IQ-TREE. The resulting trees are `PGA.upstream.treefile` and `PGA.downstream.treefile`.

## 5. Flanking-tree clustering and candidate ancestral NAHR

Raw UT and DT clusters are defined from patristic-distance matrices. The benchmark script compares the tree-only clustering procedures and HDBSCAN parameter settings:

```bash
python PGA_flanking_tree_only_clusters.py \
    --ut-tree PGA.upstream.treefile \
    --dt-tree PGA.downstream.treefile \
    --annotation PGA.structural_haplotypes.tsv \
    --outdir flanking_tree_only_clusters
```

The final analysis uses HDBSCAN with `min_cluster_size=5`. Refine the raw clusters into SH-enriched monophyletic clades with:

```bash
python PGA_flanking_tree_SH_refined_clusters.py \
    --raw-dir flanking_tree_only_clusters \
    --ut-tree PGA.upstream.treefile \
    --dt-tree PGA.downstream.treefile \
    --outdir flanking_tree_SH_refined_clusters
```

The final refinement uses exact-SH purity >= 0.80 and a minimum clade size of 10. `ancestral_NAHR_candidates.hdbscan.tsv` contains the candidate groups together with left- and right-context comparator haplotypes.

For each candidate, extract the focal *PGA* sequence and the two comparator sequences, align them with MAFFT, and convert the alignment to a focal-coordinate variant table:

```bash
mafft --auto --thread -1 candidate.fa > candidate.mafft.fa

python alignment_to_focal_variants.py \
    --alignment candidate.mafft.fa \
    --focal <focal_haplotype> \
    --left <left_context_haplotype> \
    --right <right_context_haplotype> \
    --output candidate.focal_variants.tsv
```

Candidate events and crossover intervals were then checked from the informative allele switches in this table.

## 6. Trio-based de novo NAHR

For the phased trios, infer the parental origin of each child haplotype from 10-kb bins across the +/-10-Mb *PGA* flanks:

```bash
python infer_trio_haplotype_origin.py \
    --ped ped.list \
    --root trio_assemblies \
    --out-prefix PGA.trios \
    --processes 16
```

The pedigree table should contain `Sample`, `Father` and `Mother`. For each sample, the script expects phased FASTA files named `<sample>.hap1.fa` and `<sample>.hap2.fa`; the child directory should also contain `<child>.hap1.liftoff.gff_polished` and `<child>.hap2.liftoff.gff_polished`.

The script first compares each child bin against all four parental haplotypes, assigns the transmitting parent, and then repeats the comparison between the two haplotypes of that parent. It writes a transmission summary and the bin-level parental-origin calls.

Candidate de novo events were selected by comparing these parental-origin assignments with the child and parental SHs. This step was inspected manually.

For a candidate event, prepare a tab-delimited file containing:

```text
child_hap	parent_hap
<child_haplotype>	<parental_haplotype>
```

Then generate child-reference informative-site tables for the principal-bundle gene and intergenic units with:

```bash
python nahr_mafft_to_childref_vcf.py \
    --bed PGA.principal_bundles.bed \
    --pairs candidate_pairs.tsv \
    --fasta-root trio_assemblies \
    --outdir trio_NAHR_informative_sites
```

Final de novo NAHR classification and crossover-interval assignment were made by inspecting the SH inheritance pattern together with the parent-informative allele switches.
