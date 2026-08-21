# Mammalian *PGA* analysis

Code and processed resources used for mammalian *PGA* annotation and comparative analyses.

The filtered dataset contains 479 genome assemblies from 377 eutherian species. Species-level phylogenetic analyses use the 295 species intersecting the VertLife phylogeny.

## Analysis structure

| Directory | Procedure |
| --- | --- |
| [`01_anchor_annotation/`](01_anchor_annotation/) | Extract the `VPS37C–VWCE` anchor interval, integrate TOGA2 annotations, and define local pepsinogen copy units. |
| [`02_copy_classification/`](02_copy_classification/) | Classify candidate copy units by k-mer similarity and reconstruct the protein phylogeny. |
| [`03_species_association/`](03_species_association/) | Select one assembly per species, merge ecological covariates, and run PGLS models. |
| [`04_ancestral_dynamics/`](04_ancestral_dynamics/) | Reconstruct ancestral CN states and analyse branch-level gains, losses, and expansion episodes. |
| [`repeat_context/`](repeat_context/) | Define duplicated blocks and repeat context in selected expanded lineages. |
| [`shared_resources/`](shared_resources/) | Assembly metadata, whitelists, ecological traits, species metadata, and phylogeny. |

## External data

Genome accessions, assembly sources, and quality metadata are listed in [`shared_resources/assembly_and_toga/assemblies_and_species.tsv`](shared_resources/assembly_and_toga/assemblies_and_species.tsv).

TOGA2 annotations were obtained from:

<https://genome.senckenberg.de/download/TOGA2/>

The time-calibrated mammalian phylogeny was obtained from [VertLife](https://vertlife.org/phylosubsets/). Diet and body-mass records were obtained primarily from [EltonTraits 1.0](https://figshare.com/collections/EltonTraits_1_0_Species-level_foraging_attributes_of_the_world_s_birds_and_mammals/3306933); missing diet annotations were supplemented from congeneric records or [Animal Diversity Web](https://animaldiversity.org/).

Expression support for selected reference genes was inspected in the [GTEx Portal](https://gtexportal.org/home/), [ENCODE](https://www.encodeproject.org/), and [CattleGTEx](https://cattlegtex.farmgtex.org/).

Commands, input formats, parameters, and output files are documented in the corresponding analysis directories.
