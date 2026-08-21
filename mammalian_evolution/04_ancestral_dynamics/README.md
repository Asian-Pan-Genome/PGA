# Ancestral copy-number dynamics

Commands and input tables used to reconstruct ancestral *PGA* copy number, derive branch-specific gains and losses, and run branch-level and sensitivity analyses.

## Inputs

- [`../shared_resources/phylogeny/295_sp.tree`](../shared_resources/phylogeny/295_sp.tree);
- [`../shared_resources/species_metadata/295_species_order.tsv`](../shared_resources/species_metadata/295_species_order.tsv);
- [`../shared_resources/species_metadata/295_quality.tsv`](../shared_resources/species_metadata/295_quality.tsv);
- [`../03_species_association/results/295_pga_cn_primary.tsv`](../03_species_association/results/295_pga_cn_primary.tsv);
- [`../03_species_association/results/295_diet_binary.tsv`](../03_species_association/results/295_diet_binary.tsv) and the detailed diet-state tables.

## Scripts

| Script | Procedure |
| --- | --- |
| `scripts/01_extract_species_order.py` | Generate mammalian order labels from taxonomic lineages when rebuilding the species-order table. |
| `scripts/02_reconstruct_branch_cn_and_test_rates.R` | Reconstruct ancestral CN, derive branch gains and losses, and fit count and event models. |
| `scripts/03_test_binary_branch_rates.R` | Test gain/loss counts and events for plant-dominant versus non-plant diet states on all branches or diet-stable branches. |
| `scripts/04_collapse_independent_expansion_episodes.R` | Collapse adjacent gain branches into independent expansion episodes and summarize order membership. |
| `scripts/05_leave_one_order_out.R` | Repeat branch analyses after excluding one focal mammalian order at a time. |
| `scripts/06_run_phylogenetic_null.R` | Run the phylogenetically structured null analysis with a fixed random seed. |
| `scripts/07_run_simmap_null.R` | Run conditional stochastic diet-state mapping and summarize the null analyses. |
| `scripts/08_run_high_quality_subset.R` | Repeat ancestral reconstruction and event analyses after assembly-quality filtering. |
| `scripts/09_run_high_quality_subset.sh` | Wrapper for the 10-Mbp, gap-free high-quality subset. |

## 1. Ancestral CN and branch events

Ancestral copy number is reconstructed by maximum likelihood as a continuous trait and then rounded to non-negative integer states. Branch gain and loss counts are calculated from parent and child node states.

Count models include branch duration as an offset and starting copy number as a covariate. Negative-binomial models are used when selected by the dispersion check implemented in the analysis scripts. Event-based analyses are retained in parallel.

## 2. Independent expansion episodes

Adjacent gain branches are collapsed into independent expansion episodes. The primary threshold is `CN >= 4`; `CN >= 5` is evaluated as a sensitivity threshold.

Leave-one-order-out analyses are run for:

```text
Lagomorpha
Perissodactyla
Primates
Sirenia
```

## 3. Phylogenetic null analyses

Run the phylogenetic and stochastic-mapping null procedures with:

```bash
Rscript scripts/06_run_phylogenetic_null.R --help
Rscript scripts/07_run_simmap_null.R --help
```

The scripts use the deposited species tree, copy-number table, and diet-state tables and write compact observed and null summary tables to `results/`.

## 4. High-quality assembly subset

The high-quality subset requires a gap-free anchor interval and contig N50 of at least 10 Mbp.

```bash
bash scripts/09_run_high_quality_subset.sh
```

## Outputs

The `results/` directory contains:

- reconstructed branch-event tables;
- gain/loss count-model coefficients and summaries;
- binary diet-state branch tests;
- independent expansion-episode tables;
- leave-one-order-out summaries;
- phylogenetic-null and stochastic-map summaries;
- high-quality-subset species lists, branch events, and rate-test outputs.
