# Anchor-locus annotation

Commands and scripts used to extract the mammalian `VPS37C–VWCE` anchor interval, integrate TOGA2 gene models, and define local pepsinogen copy units.

## Inputs

- [`../shared_resources/assembly_and_toga/assemblies_and_species.tsv`](../shared_resources/assembly_and_toga/assemblies_and_species.tsv): assembly accessions, sources, and quality metadata.
- [`../shared_resources/assembly_and_toga/whitelists/`](../shared_resources/assembly_and_toga/whitelists/): four-reference target and canonical-interval whitelists.
- Per assembly: TOGA2 `query_annotation.bed`, `nucleotide.fa.gz`, and `protein.fa.gz`.
- Per assembly: the corresponding genome in UCSC 2bit format.
- Merged transcript BEDs for human GRCh38, mouse GRCm38, cattle ARS-UCD2.0, and elephant mEleMax1 when rebuilding the whitelists.

## Scripts

| Script | Procedure |
| --- | --- |
| `scripts/00a_extract_reference_intervals.sh` | Extract reference `VPS37C–VWCE` intervals. |
| `scripts/00b_build_canonical_interval_whitelist.py` | Parse BED12 transcript names and retain canonical anchor intervals and target-family models. |
| `scripts/00c_assemble_four_reference_whitelists.sh` | Merge the four reference whitelists and generate transcript, source-name, and gene-label target lists. |
| `scripts/01_extract_anchor_locus.sh` | Extract the anchor interval and local TOGA2 models, record sequence gaps, and optionally generate Gepard QC files. |
| `scripts/02_refine_toga_local_units.py` | Resolve complete and partial local copy units and record fused, stretched, fragmentary, and non-target models. |
| `scripts/02_run_refine_toga_local_units.sh` | Assembly-level wrapper for local-unit refinement. |
| `scripts/03_extract_copy_sequences.sh` | Extract candidate copy sequences and add assembly identifiers to FASTA headers. |
| `scripts/04_annotate_anchor_repeats.sh` | Run RepeatMasker and TRF on the extracted anchor interval. |
| `scripts/05_merge_copy_sequences.sh` | Merge assembly-level FASTAs and copy tables for copy classification. |

## 1. Extract the anchor interval

For each assembly:

```bash
bash scripts/01_extract_anchor_locus.sh \
    ASSEMBLY_ID \
    /path/to/TOGA2/ASSEMBLY_ID \
    /path/to/ASSEMBLY_ID.2bit \
    /path/to/per_assembly
```

The script locates `VPS37C` and `VWCE`, extracts the intervening locus, records anchor or sequence-gap failures, and subsets local TOGA2 models.

## 2. Define local copy units

Run:

```bash
bash scripts/02_run_refine_toga_local_units.sh \
    ASSEMBLY_ID \
    /path/to/TOGA2/ASSEMBLY_ID \
    /path/to/per_assembly
```

Primary copy-unit seeds are compact nine-exon models spanning 4–15 Kbp. Nine-exon models up to 18 Kbp and models with at least five exons within the accepted span can be retained as supporting complete or partial units. Multi-unit, strongly stretched, or attached-fragment models are recorded but are not counted as additional copies. Only models matching the four-reference whitelist can seed a target-family unit.

## 3. Extract candidate sequences

```bash
bash scripts/03_extract_copy_sequences.sh \
    ASSEMBLY_ID \
    /path/to/TOGA2/ASSEMBLY_ID \
    /path/to/per_assembly
```

Merge assembly-level outputs for downstream classification with:

```bash
bash scripts/05_merge_copy_sequences.sh \
    /path/to/per_assembly \
    /path/to/output
```

## 4. Repeat annotation

Repeat annotation of the extracted anchor interval can be generated with:

```bash
bash scripts/04_annotate_anchor_repeats.sh \
    ASSEMBLY_ID \
    /path/to/per_assembly \
    16
```

## Outputs

The deposited outputs include:

- assembly-level anchor-locus QC tables;
- the retained gap-free assembly list;
- anchor-locus length summaries;
- candidate copy tables;
- the merged candidate protein FASTA used for copy classification.
