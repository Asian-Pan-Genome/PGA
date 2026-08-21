# Repeat context of mammalian *PGA* expansions

Commands used to define duplicated blocks and repeat-associated boundaries in selected mammalian *PGA* loci.

## Requirements

- PGGB
- ODGI
- RepeatMasker
- MUMmer4 and `paftools.js`
- SAMtools and BEDTools
- Python 3 with Biopython and Matplotlib

## 1. Apes and Old World monkeys

Build the 27-sequence PGGB graph and generate `odgi untangle` projections with:

```bash
THREADS=16 bash run_primate_pggb_untangle.sh \
    apes_owms.fa \
    apes_owms
```

The primary projection uses *Ateles hybridus* with `-m 256`. The wrapper also generates the tested alternative `-m` settings and projections using *Rhinopithecus bieti*.

Call duplicated cores and endpoint repeat annotations with:

```bash
python call_untangle_duplicon_TE.py \
    --base-dir /path/to/data \
    --species-fa /path/to/data/apes_owms.fa \
    --untangle-dir /path/to/data/apes_owms \
    --output-dir /path/to/data/untangle_duplicon_TE \
    --te-window 0
```

Generate duplicon intervals, endpoint summaries, and local plots with:

```bash
python plot_untangle_duplicon_TE.py \
    --base-dir /path/to/data \
    --input-dir /path/to/data/untangle_duplicon_TE \
    --species-fa /path/to/data/apes_owms.fa \
    --tree /path/to/mammalian_tree.nwk \
    --reference Ateles_hybridus \
    --m 256 \
    --max-join-transition-len 5000
```

Representative primate and Sirenia locus panels are generated with:

```bash
python plot_fig6b_case_panels.py --help
```

## 2. Repeat annotation

RepeatMasker annotations were generated with the mammalian Dfam library:

```bash
RepeatMasker \
    -engine rmblast \
    -species mammalia \
    -gff \
    -xsmall \
    anchor_locus.fa
```

## 3. Self-alignment of expanded loci

Sirenia, Perissodactyla, and Lagomorpha anchor loci were analysed with MUMmer4 self-alignment:

```bash
nucmer \
    -t 16 \
    --maxmatch \
    --nosimplify \
    target.pga.anchor.locus.fa \
    target.pga.anchor.locus.fa \
    -p target.self_aln

paftools.js delta2paf \
    target.self_aln.delta \
    > target.self_aln.paf
```

For the Sirenia comparison, *Dugong dugon* is used as the single-copy reference for pairwise comparison with *Trichechus inunguis*.

## Outputs

The workflow generates:

- PGGB/ODGI projection files for the primate analysis;
- duplicated-core and duplicon interval tables;
- repeat annotations intersecting duplicon endpoints;
- self-alignment PAF files for selected mammalian loci;
- local structural plots and representative case panels.
