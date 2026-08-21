# Mammalian anchor-locus annotation

Synteny-aware extraction of the `VPS37C–PGA–VWCE` interval and resolution of local pepsinogen gene-copy units from TOGA2 projections.

TOGA2 models from human, mouse, cattle, and African elephant references are integrated through the whitelists in [`../shared_resources/assembly_and_toga/whitelists/`](../shared_resources/assembly_and_toga/whitelists/). Assemblies are retained only when the two anchor genes occur on the same contig and the local interval is gap-free; subsequent copy-unit QC yields the 479-genome, 377-species comparative panel.

## Scripts

| Script | Role |
| --- | --- |
| `scripts/00a_extract_reference_intervals.sh` | Extract reference `VPS37C–VWCE` intervals. |
| `scripts/00b_build_canonical_interval_whitelist.py` | Build reference interval and target-family whitelists. |
| `scripts/00c_assemble_four_reference_whitelists.sh` | Merge the four reference whitelists. |
| `scripts/01_extract_anchor_locus.sh` | Extract the anchor interval and local TOGA2 models and record locus QC. |
| `scripts/02_refine_toga_local_units.py` | Resolve complete and partial local gene-copy units. |
| `scripts/02_run_refine_toga_local_units.sh` | Apply local-unit refinement across assemblies. |
| `scripts/03_extract_copy_sequences.sh` | Extract candidate copy sequences for classification. |
| `scripts/04_annotate_anchor_repeats.sh` | Annotate repeats in the anchor interval. |
| `scripts/05_merge_copy_sequences.sh` | Merge assembly-level copy sequences and tables. |

## Assembly-level run

```bash
bash scripts/01_extract_anchor_locus.sh \
    ASSEMBLY_ID /path/to/TOGA2/ASSEMBLY_ID /path/to/ASSEMBLY_ID.2bit /path/to/per_assembly

bash scripts/02_run_refine_toga_local_units.sh \
    ASSEMBLY_ID /path/to/TOGA2/ASSEMBLY_ID /path/to/per_assembly

bash scripts/03_extract_copy_sequences.sh \
    ASSEMBLY_ID /path/to/TOGA2/ASSEMBLY_ID /path/to/per_assembly
```

The merged candidate FASTA and copy-unit tables provide the input to [`../02_copy_classification/`](../02_copy_classification/).
