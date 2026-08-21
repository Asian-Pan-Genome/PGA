# Copy classification

Commands used to classify candidate mammalian pepsinogen copy units by alignment-free k-mer similarity and to reconstruct the intact-protein phylogeny.

## Inputs

- merged candidate-copy FASTA from [`../01_anchor_annotation/`](../01_anchor_annotation/);
- assembly-level candidate-ID lists and locus BED files;
- the primary cluster-assignment table in `results/`;
- candidate protein FASTA for intact-ORF filtering and phylogeny reconstruction.

## 1. K-mer classification

Classification was repeated for `k = 13, 15, 17, 19`; `k = 15` was used for the primary analysis.

```bash
python scripts/01_run_kmer_classification.py \
    --fasta /path/to/all_candidate_units.fasta \
    --outdir /path/to/kmer_output \
    --k_values 13,15,17,19
```

For each value of `k`, the script:

1. generates canonical presence/absence k-mer sets and excludes ambiguous sequence windows;
2. removes k-mers present in more than 95% of candidate sequences;
3. converts exact pairwise Jaccard similarities to Mash distances;
4. performs average-linkage UPGMA clustering;
5. calculates principal-coordinate analysis coordinates and cross-k summaries.

## 2. Assign copy classes

| Script | Procedure |
| --- | --- |
| `scripts/02_update_cluster_assignments.py` | Convert the selected numeric clusters to the copy-class labels used downstream. |
| `scripts/03_assign_cluster_labels_to_units.py` | Join selected class labels to local units and generate numbered assembly-level copy labels. |
| `scripts/03_run_cluster_assignment_by_assembly.sh` | Apply the deposited assignment table to each assembly. |

Apply copy assignments at the assembly level with:

```bash
bash scripts/03_run_cluster_assignment_by_assembly.sh \
    /path/to/per_assembly \
    /path/to/cluster_assignments.tsv \
    /path/to/output
```

## 3. Intact-protein phylogeny

Remove protein sequences containing internal stop codons:

```bash
python scripts/04_filter_intact_orf.py \
    --input candidate_proteins.fasta \
    --output intact_orf_proteins.fasta
```

Build the protein tree with:

```bash
bash scripts/05_build_protein_phylogeny.sh \
    intact_orf_proteins.fasta \
    /path/to/tree_output \
    16
```

The phylogeny used MAFFT v7.505 with `--auto`, trimAl v1.4.rev15, and IQ-TREE v2.1.4 with automatic model selection, 1,000 ultrafast bootstrap replicates, and 1,000 SH-aLRT replicates.

## Outputs

The deposited outputs include:

- cluster assignments and PCoA coordinates for `k = 13, 15, 17, 19`;
- the selected copy-class assignment table;
- assembly-level copy labels and copy-count summaries;
- candidate and intact-ORF protein FASTAs;
- aligned and trimmed protein FASTAs;
- IQ-TREE output files.
