# Mammalian *PGA* copy classification

Classification of local pepsinogen copy units into canonical *PGA* and divergent *PGA*-like classes using alignment-free k-mer structure, with protein phylogeny as an independent sequence-based check.

## K-mer classification

`scripts/01_run_kmer_classification.py` constructs strand-invariant k-mer presence/absence profiles, removes k-mers present in >95% of sequences, converts pairwise Jaccard similarities to Mash distances, and generates UPGMA and PCoA representations.

```bash
python scripts/01_run_kmer_classification.py \
    --fasta /path/to/all_candidate_units.fasta \
    --outdir /path/to/kmer_output \
    --k_values 13,15,17,19
```

The primary analysis uses `k = 15`; `k = 13, 17, 19` provide sensitivity analyses. Across 2,902 refined units, the primary classification contains 647 canonical *PGA*, 1,688 *PAG*, 566 *Pepf*-like, and one additional divergent *PGA*-like sequence.

## Copy labels and protein phylogeny

| Script | Role |
| --- | --- |
| `scripts/02_update_cluster_assignments.py` | Convert primary clusters to biological copy classes. |
| `scripts/03_assign_cluster_labels_to_units.py` | Apply selected classes to assembly-level copy units. |
| `scripts/03_run_cluster_assignment_by_assembly.sh` | Assembly-level wrapper for copy assignment. |
| `scripts/04_filter_intact_orf.py` | Retain intact protein sequences for phylogeny. |
| `scripts/05_build_protein_phylogeny.sh` | Run MAFFT, trimAl, and IQ-TREE. |

```bash
python scripts/04_filter_intact_orf.py \
    --input candidate_proteins.fasta \
    --output intact_orf_proteins.fasta

bash scripts/05_build_protein_phylogeny.sh \
    intact_orf_proteins.fasta /path/to/tree_output 16
```

The resulting class assignments provide species-level *PGA* CN for downstream comparative analyses.
