# Mammalian evolution

This directory contains the code and compact analysis products used for the mammalian *PGA* annotation and comparative analyses described in Supplementary Text 11. The released workflow covers the analysis from mammalian genome assembly selection and syntenic-anchor definition through independent local expansion episodes and leave-one-order-out analyses, including the k-mer sensitivity analysis.

## Directory structure

```text
mammalian_evolution/
├── shared_resources/          # Files shared across different analyses.
├── 01_anchor_annotation/      # VPS37C–VWCE locus extraction and local copy units
├── 02_copy_classification/    # k-mer classes, sensitivity analysis and protein tree
├── 03_species_association/    # species/ecology table preparation and PGLS
├── 04_ancestral_dynamics/     # ancestral CN, branch rates, episodes, LOO and null tests
└── repeat_context/            # repeat-context workflow
```

Each subfolder has its own README with inputs, scripts and outputs:

1. [Anchor annotation](01_anchor_annotation/README.md)
2. [Copy classification](02_copy_classification/README.md)
3. [Species association](03_species_association/README.md)
4. [Ancestral dynamics](04_ancestral_dynamics/README.md)


## External data acquisition

### Genome assemblies and TOGA2 annotations

Assembly accessions, source repositories and quality metadata are listed in [`shared_resources/assembly_and_toga/assemblies_and_species.tsv`](shared_resources/assembly_and_toga/assemblies_and_species.tsv). NCBI accessions can be obtained through [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/genome/); rows with a non-NCBI source retain their original source URL in the same table.

TOGA2 annotations were downloaded from:

<https://genome.senckenberg.de/download/TOGA2/>

For each assembly, the first stage expects the TOGA2 `query_annotation.bed`, `nucleotide.fa.gz` and `protein.fa.gz` files plus the corresponding genome in UCSC 2bit format. The deposited four-reference whitelists were generated from author-curated transcript BEDs for human GRCh38, mouse GRCm38, cattle ARS-UCD2.0 and elephant mEleMax1; those full BED collections are not mirrored.

### Phylogeny and ecological traits

The time-calibrated mammalian phylogeny was obtained from [VertLife](https://vertlife.org/phylosubsets/). Diet and body-mass information was derived from [EltonTraits 1.0](https://figshare.com/collections/EltonTraits_1_0_Species-level_foraging_attributes_of_the_world_s_birds_and_mammals/3306933), with documented author curation for missing species using congeneric information or [Animal Diversity Web](https://animaldiversity.org/). Only the analysis subset and authoritative curation inputs required by the released scripts are deposited.

### Expression support

Expression support was inspected directly from the [GTEx Portal](https://gtexportal.org/home/), [ENCODE](https://www.encodeproject.org/) and the [CattleGTEx portal](https://cattlegtex.farmgtex.org/). These sources support the human *PGA5*, mouse *Pepf* and cattle *PAG1* expression comparisons described in the manuscript.

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
