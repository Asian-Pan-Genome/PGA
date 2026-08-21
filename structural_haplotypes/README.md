# Structural haplotypes and NAHR

Structural decomposition of the human *PGA* locus, structural-haplotype (SH) assignment, duplicon definition, and inference of recurrent NAHR.

## Principal bundles and SHs

Run PGR-TK on the 1,348 gap-free *PGA* haplotypes and curate the resulting bundles:

```bash
pgr-pbundle-decomp \
    -r 4 \
    --min-span 12 \
    --bundle-length-cutoff 10 \
    --min-branch-size 16 \
    PGA_region.000.nogap.fa \
    PGA_48_56_4_12_10_16.

python curate_principal_bundles.py \
    PGA_48_56_4_12_10_16.bed \
    PGA.principal_bundles.bed

python assign_structural_haplotypes.py \
    PGA.principal_bundles.bed \
    PGA34A_B.tsv \
    PGA.structural_haplotypes.tsv
```

The final analysis resolved 36 SHs. A, B, and C denote *PGA34A*, *PGA34B*, and *PGA5*; X and Y denote hybrid *PGA34*/*PGA5* and *PGA5*/*PGA34* structures.

SHs in the main structural display are ordered by UPGMA clustering of PGGB path-Jaccard distances:

```bash
python pggb_jaccard_upgma.py \
    --distance pggb.jaccard.distance.tsv \
    --haplotypes PGA.structural_haplotypes.tsv \
    --output PGA.jaccard.upgma.newick
```

## Duplicons and structural display

Duplicated cores were inferred from self-alignments of one representative haplotype per SH. For each representative:

```bash
nucmer --maxmatch --nosimplify <haplotype>.fa <haplotype>.fa -p <haplotype>.self_aln
paftools.js delta2paf <haplotype>.self_aln.delta > <haplotype>.self_aln.paf
RepeatMasker -pa 16 -dir . -species human <haplotype>.fa
```

Call gene-spanning duplicated cores, extend them through the following internal spacer to define duplicons, and record repeats directly overlapping duplicon endpoints:

```bash
Rscript call_human_duplicons.R \
    --manifest representative_haplotypes.tsv \
    --gene-track representative.gene_track.bed \
    --output-prefix PGA
```

The structural-haplotype plot combines principal bundles, gene annotations, duplicons, and the Jaccard-distance tree:

```bash
Rscript plot_structural_haplotypes.R \
    representative.principal_bundles.bed \
    representative.gene_track.bed \
    PGA.duplicons.tsv \
    PGA.duplicon_endpoint_repeats.tsv \
    PGA.jaccard.upgma.newick \
    PGA.structural_haplotypes.pdf
```

## Flanking phylogenies and recent NAHR candidates

The upstream tree (UT) uses GRCh38 `chr11:61191000-61203514`; the downstream tree (DT) uses `chr11:61251444-61263958`.

```bash
bash build_flanking_phylogeny.sh \
    graph.vcf.gz GRCh38.fa chr11:61191000-61203514 \
    chimpanzee.upstream.fa PGA.upstream 16

bash build_flanking_phylogeny.sh \
    graph.vcf.gz GRCh38.fa chr11:61251444-61263958 \
    chimpanzee.downstream.fa PGA.downstream 16
```

Raw flanking clusters are inferred from patristic distances with HDBSCAN (`min_cluster_size = 5`) and refined to monophyletic clades with exact-SH purity ≥ 0.80 and at least 10 haplotypes:

```bash
python PGA_flanking_tree_only_clusters.py \
    --ut-tree PGA.upstream.treefile \
    --dt-tree PGA.downstream.treefile \
    --annotation PGA.structural_haplotypes.tsv \
    --outdir flanking_tree_only_clusters

python PGA_flanking_tree_SH_refined_clusters.py \
    --raw-dir flanking_tree_only_clusters \
    --ut-tree PGA.upstream.treefile \
    --dt-tree PGA.downstream.treefile \
    --outdir flanking_tree_SH_refined_clusters
```

Candidates nominated from UT/DT discordance were evaluated by local multiple-sequence alignment. `alignment_to_focal_variants.py` converts each alignment to a focal-coordinate informative-site table; crossover locations are reported as intervals bounded by allele switches rather than single-base breakpoints.

## Trio-based de novo NAHR

Parental origin of child haplotypes was inferred from 10-Kbp bins across ±10 Mbp around the locus:

```bash
python infer_trio_haplotype_origin.py \
    --ped ped.list \
    --root trio_assemblies \
    --out-prefix PGA.trios \
    --processes 16
```

Candidate de novo structural transitions were evaluated against parental SHs. Informative sites in homologous genic and intergenic units were generated with:

```bash
python nahr_mafft_to_childref_vcf.py \
    --bed PGA.principal_bundles.bed \
    --pairs candidate_pairs.tsv \
    --fasta-root trio_assemblies \
    --outdir trio_NAHR_informative_sites
```

Final event classification and crossover intervals were based on SH inheritance together with parent-informative allele switches.
