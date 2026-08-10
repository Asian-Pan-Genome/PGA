# Mammalian evolution

This directory contains comparative analyses of *PGA* evolution across mammals, including mammalian gene annotation, phylogenetic analyses and repeat-context analyses of lineage-specific *PGA* expansions.

## Requirements

- PGGB
- ODGI
- RepeatMasker
- MUMmer4
- `paftools.js`
- `samtools`
- `bedtools`
- Python 3
- R

## Mammalian *PGA* annotation and comparative analyses

Mammalian *PGA* loci were annotated from TOGA2 gene models within the syntenic `VPS37C–VWCE` interval. TOGA2 annotations were downloaded from:

<https://genome.senckenberg.de/download/TOGA2/>

The current analysis includes integration of TOGA2 models from multiple reference genomes, classification of canonical *PGA* and *PGA*-like gene copies, representative-genome selection, and phylogenetic analyses of *PGA* copy-number evolution and dietary ecology. The time-calibrated mammalian phylogeny was obtained from VertLife:

<https://vertlife.org/phylosubsets/>

The scripts and commands for this part will be added after the collaborator workflow is finalized.

## Repeat context

Repeat-context analyses are provided in [`repeat_context/`](repeat_context/).

### Apes and Old World monkeys

Build the 27-sequence PGGB graph and generate `odgi untangle` projections with:

```bash
THREADS=16 bash repeat_context/run_primate_pggb_untangle.sh \
    apes_owms.fa \
    apes_owms
```

The wrapper generates projections onto *Ateles hybridus* and *Rhinopithecus bieti* using `-m` values of 128, 256, 500, 1000 and 2000 bp. The primary analysis uses *Ateles hybridus* with `-m 256`.

Call duplicated cores, candidate duplicons and endpoint repeat annotations with:

```bash
python repeat_context/call_untangle_duplicon_TE.py \
    --base-dir /path/to/data \
    --species-fa /path/to/data/apes_owms.fa \
    --untangle-dir /path/to/data/apes_owms \
    --output-dir /path/to/data/untangle_duplicon_TE \
    --te-window 0
```

Generate the final duplicon intervals and repeat summary with:

```bash
python repeat_context/plot_untangle_duplicon_TE.py \
    --base-dir /path/to/data \
    --input-dir /path/to/data/untangle_duplicon_TE \
    --species-fa /path/to/data/apes_owms.fa \
    --tree /path/to/mammalian_tree.nwk \
    --reference Ateles_hybridus \
    --m 256 \
    --max-join-transition-len 5000 \
    --output-prefix /path/to/output/Ateles_hybridus.m256
```

Representative primate and Sirenia panels can be generated with `plot_fig6b_case_panels.py`.

### RepeatMasker annotation

RepeatMasker annotations were generated using the mammalian Dfam library:

```bash
RepeatMasker \
    -engine rmblast \
    -species mammalia \
    -gff \
    -xsmall \
    anchor_locus.fa
```

### Self- and pairwise-alignment analyses

For Sirenia, Perissodactyla and Lagomorpha, local homology was examined with MUMmer4 self-alignment:

```bash
nucmer -t 16 --maxmatch --nosimplify \
    target.pga.anchor.locus.fa \
    target.pga.anchor.locus.fa \
    -p target.self_aln

paftools.js delta2paf target.self_aln.delta > target.self_aln.paf
```

Plot and summarize the self-alignment with:

```bash
Rscript repeat_context/plot_mammal_self_alignment.R \
    target.self_aln.paf \
    -o target.pga.self_aln \
    --rm target.pga.anchor.locus.fa.out \
    --gene-bed toga.PGA_like.local.v4.assign_candidate_ids.bed \
    --anchor-bed pga.anchor.locus.bed \
    --merge-dist 2000 \
    --te-window 0 \
    -q 4000 -m 4000 -r 4000
```

Where a single-copy outgroup is available, pairwise alignment can be analysed with `plot_mammal_pairwise_alignment.R`. For Sirenia, *Dugong dugon* was used as the single-copy outgroup for *Trichechus inunguis*.

The directory also contains dataset-specific plotting scripts for representative primate, Sirenia, Perissodactyla and Lagomorpha loci.
