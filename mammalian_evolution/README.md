# Mammalian evolution

Comparative analysis of mammalian *PGA* evolution, spanning synteny-aware annotation and copy classification, phylogenetic association, ancestral copy-number dynamics, and repeat context of lineage-specific expansions.

The filtered comparative panel contains 479 genomes from 377 eutherian species; phylogeny-informed analyses use 295 species matched to the VertLife tree.

## Analysis structure

| Directory | Analysis |
| --- | --- |
| [`01_anchor_annotation/`](01_anchor_annotation/) | Extract the `VPS37C–PGA–VWCE` syntenic interval and resolve local pepsinogen copy units. |
| [`02_copy_classification/`](02_copy_classification/) | Classify canonical *PGA* and divergent *PGA*-like copies using k-mers and protein phylogeny. |
| [`03_species_association/`](03_species_association/) | Select one representative genome per species and test extant CN–ecology associations with PGLS. |
| [`04_ancestral_dynamics/`](04_ancestral_dynamics/) | Reconstruct ancestral *PGA* CN and analyse branch-level gains, losses, and expansion episodes. |
| [`repeat_context/`](repeat_context/) | Resolve duplicated blocks and repeat-associated boundaries in expanded lineages. |
| [`shared_resources/`](shared_resources/) | Assembly metadata, ecological traits, species metadata, and the 295-species phylogeny. |

## Data sources

Genome accessions and assembly metadata are listed in [`shared_resources/assembly_and_toga/assemblies_and_species.tsv`](shared_resources/assembly_and_toga/assemblies_and_species.tsv). Local gene models were derived from [TOGA2](https://genome.senckenberg.de/download/TOGA2/). Phylogenetic analyses use [VertLife](https://vertlife.org/phylosubsets/); dietary traits and body mass are based primarily on [EltonTraits 1.0](https://figshare.com/collections/EltonTraits_1_0_Species-level_foraging_attributes_of_the_world_s_birds_and_mammals/3306933), supplemented where necessary from [Animal Diversity Web](https://animaldiversity.org/).
