# Structural haplotypes

Code used to resolve structural haplotypes (SHs) of the human *PGA* locus, characterize its duplicated architecture, and examine recurrent NAHR.

## Requirements

- [PGR-TK](https://github.com/Sema4-Research/pgr-tk)
- [ODGI](https://github.com/pangenome/odgi)
- [MUMmer4](https://github.com/mummer4/mummer) and `paftools.js`
- [RepeatMasker](https://www.repeatmasker.org/)
- Python 3 with `numpy`, `pandas`, `scipy`, `ete3`, `scikit-learn` and `hdbscan`
- R with `optparse`, `tidyverse`, `ggtree`, `gggenes`, `aplot`, `RColorBrewer` and `ape`
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

Merge the raw bundle calls into the principal-bundle representation:

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

`PGA34A_B.tsv` should contain `Sample`, `Hap`, `Source` and `PGAs`. The final analysis contains 36 SHs. A, B and C denote *PGA34A*, *PGA34B* and *PGA5*, respectively; X and Y denote hybrid *PGA34*/*PGA5* and *PGA5*/*PGA34* genes.

## 2. PGGB Jaccard-distance tree

The SH plot is ordered by UPGMA clustering of pairwise path Jaccard distances from the PGGB graph. Generate the distances with `odgi similarity -d`, then run:

```bash
python pggb_jaccard_upgma.py \
    --distance pggb.jaccard.distance.tsv \
    --haplotypes PGA.structural_haplotypes.tsv \
    --output PGA.jaccard.upgma.newick
```

The ODGI table should contain `group.a`, `group.b` and `jaccard.distance`.

## 3. Gene track

Generate the *PGA* gene track with:

```bash
python generate_gene_track.py \
    PGA_region.000.hit \
    id.list \
    PGA34A_B.tsv \
    PGA.gene_track.bed
```

`id.list` should contain `Sample`, `Hap` and `GFF File`.

## 4. Duplicons and repeat annotation

Duplicons were defined from self-alignments of one representative haplotype from each SH. For each representative haplotype, run:

```bash
nucmer --maxmatch --nosimplify \
    <haplotype>.fa <haplotype>.fa \
    -p <haplotype>.self_aln

paftools.js delta2paf \
    <haplotype>.self_aln.delta \
    > <haplotype>.self_aln.paf

RepeatMasker -pa 16 -dir . -species human <haplotype>.fa
```

Prepare a tab-delimited manifest:

```text
haplotype	paf	repeatmasker
AFG05.hap1	/path/to/AFG05.hap1.self_aln.paf	/path/to/AFG05.hap1.fa.out
BAN02.hap2	/path/to/BAN02.hap2.self_aln.paf	/path/to/BAN02.hap2.fa.out
```

Call duplicated cores, final duplicons and repeats directly overlapping their endpoints with:

```bash
Rscript call_human_duplicons.R \
    --manifest representative_haplotypes.tsv \
    --gene-track representative.gene_track.bed \
    --output-prefix PGA
```

The caller merges locally collinear self-alignment blocks separated by no more than 1 kb, removes alignments shorter than the 2-kb threshold and alignments within 5 kb of the self-alignment diagonal, and intersects the reference- and query-axis projections. For each *PGA* gene, the shortest qualifying core spanning the complete gene body is retained; candidates extending into an adjacent *PGA* gene are excluded.

Each duplicated core is extended to the beginning of the next downstream core to define the final duplicon. The terminal sequence downstream of the distal *PGA5* core is not included. RepeatMasker annotations are reported only when they directly overlap a duplicon boundary.

The main outputs are:

```text
PGA.duplicons.tsv
PGA.duplicon_endpoint_repeats.tsv
```

## 5. Representative SH plot

After selecting one representative haplotype for each SH, subset the principal-bundle and gene tracks to those representatives and plot them together with the duplicon tracks:

```bash
Rscript plot_structural_haplotypes.R \
    representative.principal_bundles.bed \
    representative.gene_track.bed \
    PGA.duplicons.tsv \
    PGA.duplicon_endpoint_repeats.tsv \
    PGA.jaccard.upgma.newick \
    PGA.structural_haplotypes.pdf
```

Representative selection is used only for visualization and is not part of SH assignment.

## 6. Upstream and downstream flanking trees

The upstream tree (UT) uses GRCh38 `chr11:61191000-61203514`, and the downstream tree (DT) uses `chr11:61251444-61263958`.

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

## 7. Flanking-tree clustering and candidate ancestral NAHR

Benchmark tree-only clustering procedures and HDBSCAN parameter settings with:

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

The final refinement uses exact-SH purity >= 0.80 and a minimum clade size of 10.

For each candidate ancestral NAHR event, align the focal haplotype with the two flanking-context comparators and convert the alignment to a focal-coordinate variant table:

```bash
mafft --auto --thread -1 candidate.fa > candidate.mafft.fa

python alignment_to_focal_variants.py \
    --alignment candidate.mafft.fa \
    --focal <focal_haplotype> \
    --left <left_context_haplotype> \
    --right <right_context_haplotype> \
    --output candidate.focal_variants.tsv
```

Candidate events and crossover intervals were checked from the informative allele switches.

## 8. Trio-based de novo NAHR

Infer the parental origin of each child haplotype from 10-kb bins across the +/-10-Mb *PGA* flanks:

```bash
python infer_trio_haplotype_origin.py \
    --ped ped.list \
    --root trio_assemblies \
    --out-prefix PGA.trios \
    --processes 16
```

The pedigree table should contain `Sample`, `Father` and `Mother`. The script first compares each child bin against all four parental haplotypes, assigns the transmitting parent, and then repeats the comparison between the two haplotypes of that parent.

Candidate de novo events were selected by comparing the parental-origin assignments with the child and parental SHs. This step was inspected manually.

For each candidate event, prepare:

```text
child_hap	parent_hap
<child_haplotype>	<parental_haplotype>
```

Then generate child-reference informative-site tables for principal-bundle gene and intergenic units:

```bash
python nahr_mafft_to_childref_vcf.py \
    --bed PGA.principal_bundles.bed \
    --pairs candidate_pairs.tsv \
    --fasta-root trio_assemblies \
    --outdir trio_NAHR_informative_sites
```

Final de novo NAHR classification and crossover-interval assignment were made by inspecting the SH inheritance pattern together with the parent-informative allele switches.
