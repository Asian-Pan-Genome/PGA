# Mammalian evolution

This directory contains comparative analyses of *PGA* evolution across mammals, including mammalian gene annotation, phylogenetic analyses and repeat-context analyses of lineage-specific *PGA* expansions.

## Mammalian *PGA* annotation and comparative analyses

Mammalian *PGA* loci were annotated from TOGA2 gene models within the syntenic `VPS37C–VWCE` interval. TOGA2 annotations were downloaded from:

<https://genome.senckenberg.de/download/TOGA2/>

The analysis includes integration of TOGA2 models from multiple reference genomes, classification of canonical *PGA* and *PGA*-like genes, representative-genome selection, and phylogenetic analyses of *PGA* copy-number evolution and dietary ecology. The mammalian phylogeny was obtained from VertLife:

<https://vertlife.org/phylosubsets/>

Scripts and commands for this part will be added after the collaborator workflow is finalized.

## Repeat context

### Requirements

- PGGB
- ODGI
- RepeatMasker
- MUMmer4 and `paftools.js`
- `samtools` and `bedtools`
- Python 3 with `matplotlib` and `biopython`

### Apes and Old World monkeys

Build the 27-sequence PGGB graph and generate `odgi untangle` projections with:

```bash
THREADS=16 bash repeat_context/run_primate_pggb_untangle.sh apes_owms.fa apes_owms
```

The primary analysis uses *Ateles hybridus* with `-m 256`; the wrapper also generates the other tested `-m` values and *Rhinopithecus bieti* projections.

Call duplicated cores and endpoint repeat annotations with:

```bash
python repeat_context/call_untangle_duplicon_TE.py \
    --base-dir /path/to/data \
    --species-fa /path/to/data/apes_owms.fa \
    --untangle-dir /path/to/data/apes_owms \
    --output-dir /path/to/data/untangle_duplicon_TE \
    --te-window 0
```

Generate the final duplicon intervals, repeat summary and plots with:

```bash
python repeat_context/plot_untangle_duplicon_TE.py \
    --base-dir /path/to/data \
    --input-dir /path/to/data/untangle_duplicon_TE \
    --species-fa /path/to/data/apes_owms.fa \
    --tree /path/to/mammalian_tree.nwk \
    --reference Ateles_hybridus \
    --m 256 \
    --max-join-transition-len 5000
```

Representative primate and Sirenia panels are generated with `plot_fig6b_case_panels.py`.

### Repeat annotation and other mammalian lineages

RepeatMasker annotations were generated using the mammalian Dfam library. Sirenia, Perissodactyla and Lagomorpha loci were examined with MUMmer4 self-alignment; *Dugong dugon* was additionally used as a single-copy outgroup for pairwise comparison with *Trichechus inunguis*.

```bash
RepeatMasker -engine rmblast -species mammalia -gff -xsmall anchor_locus.fa

nucmer -t 16 --maxmatch --nosimplify \
    target.pga.anchor.locus.fa target.pga.anchor.locus.fa \
    -p target.self_aln

paftools.js delta2paf target.self_aln.delta > target.self_aln.paf
```
