# Ancestral copy-number dynamics

Reconstruction of ancestral *PGA* copy number and analysis of branch-level gains and losses in relation to dietary state. Additional analyses evaluate independent expansion episodes and robustness to order composition, trait mapping, and assembly quality.

## Workflow

| Script | Role |
| --- | --- |
| `scripts/01_extract_species_order.py` | Prepare mammalian order labels. |
| `scripts/02_reconstruct_branch_cn_and_test_rates.R` | Reconstruct ancestral CN and derive branch-level gains and losses. |
| `scripts/03_test_binary_branch_rates.R` | Compare gain/loss dynamics between plant-dominant and other branches. |
| `scripts/04_collapse_independent_expansion_episodes.R` | Collapse adjacent gain branches into local expansion episodes. |
| `scripts/05_leave_one_order_out.R` | Repeat branch models after excluding focal mammalian orders. |
| `scripts/06_run_phylogenetic_null.R` | Evaluate a phylogenetically structured null. |
| `scripts/07_run_simmap_null.R` | Evaluate stochastic dietary-state mappings. |
| `scripts/08_run_high_quality_subset.R` | Repeat analyses in a higher-quality assembly subset. |
| `scripts/09_run_high_quality_subset.sh` | Wrapper for the gap-free, contig-N50 ≥10-Mbp subset. |

Ancestral CN is reconstructed as a continuous trait and rounded to non-negative integer states for branch-event summaries. Count models use branch duration as an offset and starting CN as a covariate; negative-binomial models are used when required by overdispersion. Expansion onset is evaluated at CN ≥4, with CN ≥5 as a sensitivity threshold.

Leave-one-order-out analyses evaluate Lagomorpha, Perissodactyla, Primates, and Sirenia. The high-quality subset can be rerun with:

```bash
bash 04_ancestral_dynamics/scripts/09_run_high_quality_subset.sh
```
