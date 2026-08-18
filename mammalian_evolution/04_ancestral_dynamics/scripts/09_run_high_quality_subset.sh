#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
module_dir="$(cd -- "${script_dir}/.." && pwd)"
reorganized_root="$(cd -- "${script_dir}/../.." && pwd)"

tree_file=${1:-"${reorganized_root}/shared_resources/phylogeny/295_sp.tree"}
cn_file=${2:-"${reorganized_root}/03_species_association/results/295_pga_cn_primary.tsv"}
diet_file=${3:-"${reorganized_root}/03_species_association/results/295_diet_binary.tsv"}
quality_file=${4:-"${reorganized_root}/shared_resources/species_metadata/295_quality.tsv"}
out_prefix=${5:-"${module_dir}/results/295_pga_analysis5_n50_10mb_gapfree"}
min_contig_n50=${6:-10000000}
require_gap_free=${7:-TRUE}

command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript is not on PATH" >&2; exit 1; }
for input in "${tree_file}" "${cn_file}" "${diet_file}" "${quality_file}"; do
  [[ -s "${input}" ]] || { echo "ERROR: missing input: ${input}" >&2; exit 1; }
done
mkdir -p "$(dirname "${out_prefix}")"

Rscript "${script_dir}/08_run_high_quality_subset.R" \
  "${tree_file}" "${cn_file}" "${diet_file}" "${quality_file}" \
  "${out_prefix}" "${min_contig_n50}" "${require_gap_free}"

echo "High-quality subset analysis completed: ${out_prefix}"
