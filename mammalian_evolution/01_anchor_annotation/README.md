# Stage 01: anchor-locus annotation

This stage implements the Supplementary Note 11 workflow from genome-assembly metadata and syntenic-anchor definition through local *PGA*/*PGA*-like/*PAG*-like copy-unit annotation. It extracts the gap-free interval between `VPS37C` and `VWCE`, integrates TOGA2 models from four source references, separates independent local units from fused, stretched or fragmentary projections, extracts candidate sequences and optionally annotates repeats across the anchor locus.

## Inputs

- `../shared_resources/assembly_and_toga/assemblies_and_species.tsv`: assembly accessions, sources and quality metadata.
- `../shared_resources/assembly_and_toga/whitelists/`: deposited four-reference target and canonical-interval whitelists.
- Per assembly: TOGA2 `query_annotation.bed`, `nucleotide.fa.gz` and `protein.fa.gz`.
- Per assembly: the corresponding genome as a UCSC 2bit file.
- For rebuilding the deposited whitelists only: the author-curated transcript BEDs for human GRCh38, mouse GRCm38, cattle ARS-UCD2.0 and elephant mEleMax1.

Whole assemblies, bulk TOGA2 directories and the complete four-reference transcript BED collections are not stored in GitHub. See the top-level README for download locations.

## Scripts

| Script | Role |
| --- | --- |
| `scripts/00a_extract_reference_intervals.sh` | Build one canonical `VPS37C–VWCE` interval whitelist per source reference from transcript BEDs. |
| `scripts/00b_build_canonical_interval_whitelist.py` | Parse BED12 transcript names and retain the canonical anchor interval and target-family models. |
| `scripts/00c_assemble_four_reference_whitelists.sh` | Merge the four reference-specific whitelists and produce source-transcript, source-name and gene-label target lists. |
| `scripts/01_extract_anchor_locus.sh` | Require both anchors on the same contig, extract the interval and local TOGA2 models, detect gaps, and optionally generate Gepard QC. |
| `scripts/02_refine_toga_local_units.py` | Refine local TOGA2 models into independent complete or partial units while recording fused/bridging, stretched, fragmentary and non-target rearrangement evidence. |
| `scripts/02_run_refine_toga_local_units.sh` | Assembly-level wrapper for local-unit refinement with the deposited four-reference whitelist. |
| `scripts/03_extract_copy_sequences.sh` | Extract candidate protein sequences and add assembly identifiers to FASTA headers. |
| `scripts/04_annotate_anchor_repeats.sh` | Run RepeatMasker and TRF on the extracted anchor locus. The optional TRF converter is intentionally not public. |
| `scripts/05_merge_copy_sequences.sh` | Merge assembly-level copy FASTAs and copy tables for downstream classification. |

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

The optional Gepard and legacy clustering arguments accepted by `01_extract_anchor_locus.sh` reproduce historical QC/export steps; they are not required for the released local-unit workflow.

## Unit logic and QC

The manuscript-defined logic uses compact nine-exon models spanning 4–15 kbp as primary seeds; nine-exon models up to 18 kbp and models with at least five exons in the accepted span can support complete or partial local units. Models spanning multiple local units, strongly stretched projections and fragments attached to an existing unit are recorded but are not counted as additional copies. Only models matching the four-reference target whitelist can seed a target-family unit.

The stage records assemblies where the anchors are absent, occur on different contigs, contain sequence gaps, or otherwise fail locus extraction. These QC lists are retained in `results/` so that the final assembly set is auditable.

## Deposited outputs

- the final assembly lists and anchor-locus QC summaries;
- the gap-free anchor-interval/species metadata table;
- anchor-locus length summaries;
- the merged candidate protein FASTA used by Stage 02.

Large per-assembly locus FASTAs, raw TOGA2 model directories and full RepeatMasker/TRF working outputs are not deposited.

## Dependencies

Bash, GNU core utilities, `awk`, `bc`, Python 3, `pandas`, UCSC `twoBitToFa` and `seqkit` are required. RepeatMasker plus TRF are required only for `04_annotate_anchor_repeats.sh`; Java/Gepard is optional.
