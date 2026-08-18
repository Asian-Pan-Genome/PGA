# Stage 04: ancestral copy-number dynamics

This stage reconstructs ancestral *PGA* copy number, derives branch-specific gain and loss events, tests diet-associated event rates, collapses adjacent gain branches into independent local expansion episodes, performs leave-one-order-out analyses, and evaluates phylogenetic null and high-quality assembly subsets.

## Inputs

- `../shared_resources/phylogeny/295_sp.tree`;
- `../shared_resources/species_metadata/295_species_order.tsv`;
- `../shared_resources/species_metadata/295_quality.tsv`;
- `../03_species_association/results/295_pga_cn_primary.tsv`;
- `../03_species_association/results/295_diet_binary.tsv` and the detailed diet-state tables.

## Scripts

| Script | Role |
| --- | --- |
| `scripts/01_extract_species_order.py` | Extract mammalian order labels from taxonomic lineages when rebuilding the deposited order table. |
| `scripts/02_reconstruct_branch_cn_and_test_rates.R` | Reconstruct continuous ancestral CN, round to non-negative copy states, derive branch gains/losses and fit count/event-rate models. |
| `scripts/03_test_binary_branch_rates.R` | Recode diet to plant-dominant versus non-plant and test gain/loss counts and events on all or diet-stable branches. |
| `scripts/04_collapse_independent_expansion_episodes.R` | Collapse adjacent gain branches into independent expansion episodes and summarize their order membership. |
| `scripts/05_leave_one_order_out.R` | Refit the branch models after excluding one focal mammalian order at a time. |
| `scripts/06_run_phylogenetic_null.R` | Generate the phylogenetically structured null analysis with a fixed random seed. |
| `scripts/07_run_simmap_null.R` | Run conditional stochastic diet-state mapping and summarize the corresponding null tests. |
| `scripts/08_run_high_quality_subset.R` | Repeat ancestral reconstruction and event tests after assembly-quality filtering. |
| `scripts/09_run_high_quality_subset.sh` | Convenience wrapper for the deposited 10-Mbp, gap-free high-quality subset. |

## Analysis conventions

Ancestral copy number is reconstructed by maximum likelihood as a continuous trait and rounded to non-negative states before branch gains and losses are derived. Count models include branch duration as an offset and starting copy number as a covariate; negative-binomial models are used when supported by dispersion, with complementary event tests retained in the outputs.

Independent local expansions are defined by collapsing adjacent gain-bearing branches and applying the manuscript copy-number thresholds (`CN >= 4`, with `CN >= 5` as a sensitivity threshold). Leave-one-order-out analyses evaluate Lagomorpha, Perissodactyla, Primates and Sirenia.

The high-quality sensitivity analysis requires a gap-free anchor interval and contig N50 of at least 10 Mbp:

```bash
bash 04_ancestral_dynamics/scripts/09_run_high_quality_subset.sh
```

## Deposited outputs

The `results/` directory contains:

- branch-event tables and count-model coefficients/summaries;
- binary-diet results for all branches and diet-stable branches;
- independent-expansion episode membership and order summaries;
- leave-one-order-out statistics;
- compact phylogenetic-null and stochastic-map summaries plus observed test tables;
- high-quality species lists, branch events and binary-rate tests.

Full per-replicate null simulations and stochastic maps are intentionally not deposited. The compact summaries contain the observed statistics and empirical comparisons required to audit the manuscript results. No final manuscript plotting script is included.

## Dependencies

Python 3 with `pandas`; R with `ape`, `MASS` and `phytools`; and Bash for the high-quality subset wrapper.
