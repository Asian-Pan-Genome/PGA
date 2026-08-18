#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
module_dir="$(cd -- "${script_dir}/.." && pwd)"
reorganized_root="$(cd -- "${script_dir}/../.." && pwd)"

bodymass_input=${1:-"${module_dir}/results/295_pga_diet_bodymass_best_assembly.tsv"}
diet_input=${2:-"${module_dir}/results/295_pga_diet_bodymass_high_order_ratio.tsv"}
tree_file=${3:-"${reorganized_root}/shared_resources/phylogeny/295_sp.tree"}
output_dir=${4:-"${module_dir}/results"}

command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript is not on PATH" >&2; exit 1; }
for input in "${bodymass_input}" "${diet_input}" "${tree_file}"; do
  [[ -s "${input}" ]] || { echo "ERROR: missing input: ${input}" >&2; exit 1; }
done
mkdir -p "${output_dir}"

Rscript "${script_dir}/03_run_pgls_body_mass.R"   "${bodymass_input}" "${tree_file}"   "${output_dir}/pga_diet_bodymass_pgls"

Rscript "${script_dir}/04_run_pgls_diet_models.R"   "${diet_input}" "${tree_file}"   "${output_dir}/pga_high_order_diet_pgls"

echo "PGLS models completed: ${output_dir}"
