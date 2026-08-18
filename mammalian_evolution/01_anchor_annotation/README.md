# 01: anchor-locus annotation

This folder contains scripts used to extract the gap-free `VPS37C–VWCE` interval and annotate local *PGA*, *PGA*-like and *PAG*-like copy units. It integrates TOGA2 models from four references, resolves independent units and extracts candidate sequences for classification.

## Inputs

- `../shared_resources/assembly_and_toga/assemblies_and_species.tsv`: assembly accessions, sources and quality metadata.
- `../shared_resources/assembly_and_toga/whitelists/`: deposited four-reference target and canonical-interval whitelists.
- Per assembly: TOGA2 `query_annotation.bed`, `nucleotide.fa.gz` and `protein.fa.gz`.
- Per assembly: the corresponding genome as a UCSC 2bit file.
- Merged transcript BEDs for human GRCh38, mouse GRCm38, cattle ARS-UCD2.0 and elephant mEleMax1, required only to rebuild the whitelists.


## Scripts

| Script | Role |
| --- | --- |
| `scripts/00a_extract_reference_intervals.sh` | Build a canonical `VPS37C–VWCE` interval whitelist for each reference. |
| `scripts/00b_build_canonical_interval_whitelist.py` | Parse BED12 transcript names and retain the canonical anchor interval and target-family models. |
| `scripts/00c_assemble_four_reference_whitelists.sh` | Merge reference whitelists and generate transcript, source-name and gene-label target lists. |
| `scripts/01_extract_anchor_locus.sh` | Extract the anchor interval and local TOGA2 models, detect gaps and optionally generate Gepard QC. |
| `scripts/02_refine_toga_local_units.py` | Resolve complete and partial local units and record fused, stretched, fragmentary and non-target models. |
| `scripts/02_run_refine_toga_local_units.sh` | Assembly-level wrapper for local-unit refinement with the deposited four-reference whitelist. |
| `scripts/03_extract_copy_sequences.sh` | Extract candidate proteins and add assembly identifiers to FASTA headers. |
| `scripts/04_annotate_anchor_repeats.sh` | Run RepeatMasker and TRF on the extracted anchor locus. The optional TRF converter is intentionally not public. |
| `scripts/05_merge_copy_sequences.sh` | Merge assembly-level FASTAs and copy tables for Stage 02. |

## Typical assembly-level run

```bash
bash scripts/01_extract_anchor_locus.sh \
  ASSEMBLY_ID /path/to/TOGA2/ASSEMBLY_ID /path/to/ASSEMBLY_ID.2bit /path/to/per_assembly

bash scripts/02_run_refine_toga_local_units.sh \
  ASSEMBLY_ID /path/to/TOGA2/ASSEMBLY_ID /path/to/per_assembly

bash scripts/03_extract_copy_sequences.sh \
  ASSEMBLY_ID /path/to/TOGA2/ASSEMBLY_ID /path/to/per_assembly
```

Optional repeat annotation:

```bash
bash scripts/04_annotate_anchor_repeats.sh ASSEMBLY_ID /path/to/per_assembly 16
```


## Define copy unit and QC

Primary seeds are compact nine-exon models spanning 4–15 kbp. Nine-exon models up to 18 kbp and models with at least five exons in the accepted span can support complete or partial units. Multi-unit, strongly stretched or attached fragment models are recorded but not counted as additional copies. Only four-reference whitelist matches can seed a target-family unit.

QC tables in `results/` record missing anchors, anchors on different contigs, sequence gaps and other extraction failures.

## Outputs

- final assembly lists and anchor-locus QC summaries;
- gap-free anchor-interval and species metadata;
- anchor-locus length summaries;
- the merged candidate protein FASTA (for Stage 02).

