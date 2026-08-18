# Stage 02: copy classification

This stage classifies local candidate units by alignment-free k-mer content, evaluates sensitivity across `k = 13, 15, 17, 19`, propagates the selected labels back to assembly-level units, filters intact protein sequences and reconstructs the protein phylogeny.

## Inputs

- the merged candidate FASTA from Stage 01;
- assembly-level candidate-ID lists and locus BED files produced by Stage 01;
- the deposited primary cluster assignment table in `results/`;
- the final protein FASTA for intact-ORF filtering and phylogenetic reconstruction.

## Alignment-free classification

`scripts/01_run_kmer_classification.py` applies the following workflow independently for each requested k:

1. generate canonical presence/absence k-mer sets while skipping windows containing ambiguous bases;
2. remove k-mers present in more than 95% of sequences;
3. compute exact pairwise Jaccard similarities and transform them to Mash distances;
4. perform average-linkage/UPGMA clustering;
5. generate principal-coordinate (PCoA) coordinates and clustering diagnostics.

The primary analysis uses `k = 15`. The deposited sensitivity tables cover 2,902 classified units: 647 *PGA*, 1,688 *PAG*, 566 *Pepf*-like units and one divergent *PGA*-like unit. Relative to `k = 15`, one assignment changed at `k = 17` and nine changed at `k = 19` (adjusted Rand indices 0.999 and 0.995); `k = 13` retained the four broad classes but split them into many small subclusters. These comparisons support Supplementary Table S12.

Example:

```bash
python scripts/01_run_kmer_classification.py \
  --fasta /path/to/all_candidate_units.fasta \
  --outdir /path/to/kmer_output \
  --k_values 13,15,17,19
```

The classifier writes diagnostic plots together with its matrices and tables. Only the compact assignment and PCoA tables used for audit and downstream analysis are deposited; no final manuscript plotting script is included.

## Scripts

| Script | Role |
| --- | --- |
| `scripts/01_run_kmer_classification.py` | Exact canonical k-mer, Mash-distance, UPGMA and PCoA workflow with cross-k sensitivity summaries. |
| `scripts/02_update_cluster_assignments.py` | Convert numeric primary clusters to the biological labels used downstream. |
| `scripts/03_assign_cluster_labels_to_units.py` | Join selected labels to candidate local units and create numbered assembly-level copy labels. |
| `scripts/03_run_cluster_assignment_by_assembly.sh` | Assembly-level wrapper for applying the deposited assignment table. |
| `scripts/04_filter_intact_orf.py` | Remove protein sequences with internal stop codons while allowing a terminal stop. |
| `scripts/05_build_protein_phylogeny.sh` | Run MAFFT, trimAl and IQ-TREE with model selection and 1,000 ultrafast bootstrap replicates. |

Example protein-tree run:

```bash
python scripts/04_filter_intact_orf.py \
  --input candidate_proteins.fasta \
  --output intact_orf_proteins.fasta

bash scripts/05_build_protein_phylogeny.sh \
  intact_orf_proteins.fasta /path/to/tree_output 16
```

The reported phylogeny used MAFFT 7.505 (`--auto`), trimAl 1.4.rev15 and IQ-TREE 2.1.4 with automatic model selection, 1,000 ultrafast bootstrap replicates and 1,000 SH-aLRT replicates.

## Deposited outputs

- cluster assignments and PCoA coordinates for `k = 13, 15, 17, 19`;
- the selected biological-label assignment table and copy-count summary;
- merged candidate, intact-ORF, aligned and trimmed FASTAs;
- final tree files and the retained Newick trees.

Large distance matrices, diagnostic graphics and ordinary IQ-TREE runtime/report files are not part of the GitHub release.

## Expression support

Human *PGA5*, mouse *Pepf* and cattle *PAG1* expression support was evaluated from the GTEx, mouse ENCODE and CattleGTEx web resources cited in the top-level README. The repository does not mirror or reprocess those external expression matrices.

## Dependencies

Python 3 with `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib` and Biopython; Bash; MAFFT; trimAl; and IQ-TREE 2 (`iqtree2`).
