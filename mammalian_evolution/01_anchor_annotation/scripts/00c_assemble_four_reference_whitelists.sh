#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
reorganized_root="$(cd -- "${script_dir}/../.." && pwd)"
whitelist_dir=${1:-"${reorganized_root}/shared_resources/assembly_and_toga/whitelists"}

combined="${whitelist_dir}/four_ref_vps37c_vwce_canonical_interval_whitelist.tsv"
target_source="${whitelist_dir}/four_ref_pga_like_target_source_whitelist.tsv"

inputs=(
  "${whitelist_dir}/hg38_vps37c_vwce_canonical_interval_whitelist.tsv"
  "${whitelist_dir}/mm10_vps37c_vwce_canonical_interval_whitelist.tsv"
  "${whitelist_dir}/hlbostau10_vps37c_vwce_canonical_interval_whitelist.tsv"
  "${whitelist_dir}/hlelemaxind3a_vps37c_vwce_canonical_interval_whitelist.tsv"
)

for input in "${inputs[@]}"; do
  [[ -s "${input}" ]] || { echo "ERROR: missing whitelist: ${input}" >&2; exit 1; }
done

awk 'FNR==1 && NR!=1 {next} {print}' "${inputs[@]}" > "${combined}"

awk -F'\t' '
NR==1 {print; next}
$9!="VPS37C" && $9!="VWCE" && $9!="Vps37c" && $9!="Vwce" && $9!="CD5" {print}
' "${combined}" > "${target_source}"

awk -F'\t' 'NR>1{print $8}' "${target_source}" |
  sort -u > "${whitelist_dir}/four_ref_pga_like_source_transcript_id_whitelist.txt"
awk -F'\t' 'NR>1{print $10}' "${target_source}" |
  sort -u > "${whitelist_dir}/four_ref_pga_like_source_name_whitelist.txt"
awk -F'\t' 'NR>1{print $9}' "${target_source}" |
  sort -u > "${whitelist_dir}/four_ref_pga_like_gene_label_whitelist.txt"

echo "Combined whitelists written to: ${whitelist_dir}"
