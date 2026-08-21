# Repeat context of mammalian *PGA* expansions

Structural analysis of duplicated units and repeat-associated boundaries in mammalian lineages with expanded *PGA* loci.

## Apes and Old World monkeys

A PGGB graph is projected onto a single-copy outgroup with `odgi untangle`:

```bash
THREADS=16 bash run_primate_pggb_untangle.sh apes_owms.fa apes_owms
```

The primary projection uses *Ateles hybridus* with `-m 256`; neighbouring merge-distance settings and the single-copy Old World monkey *Rhinopithecus bieti* provide sensitivity comparisons.

Call candidate duplicated cores and endpoint repeat annotations with:

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

`plot_fig6b_case_panels.py` generates representative primate and Sirenia locus panels.

## Other mammalian lineages

Sirenia, Perissodactyla, and Lagomorpha loci are examined with MUMmer4 self-alignment and RepeatMasker. Where a suitable single-copy outgroup is available, pairwise alignment is used to compare duplicated-block boundaries; *Dugong dugon* provides the outgroup comparison for *Trichechus inunguis*.

```bash
RepeatMasker -engine rmblast -species mammalia -gff -xsmall anchor_locus.fa
nucmer -t 16 --maxmatch --nosimplify target.pga.anchor.locus.fa target.pga.anchor.locus.fa -p target.self_aln
paftools.js delta2paf target.self_aln.delta > target.self_aln.paf
```

These analyses describe structural and repeat context; repeat overlap alone is not treated as proof of a specific duplication mechanism.
