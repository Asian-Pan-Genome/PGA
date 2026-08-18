#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 <assembly> <per_assembly_output_root> [cluster_assignments.tsv]" >&2
  exit 2
fi

assembly=$1
output_root=$2
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
module_dir="$(cd -- "${script_dir}/.." && pwd)"
assignments=${3:-"${module_dir}/results/cluster_assignments_updated_v2.tsv"}
outdir="${output_root}/${assembly}"
subset="${outdir}/${assembly}.cluster_assignments.updated.v2.tsv"

[[ -s "${assignments}" ]] || { echo "ERROR: missing assignment table: ${assignments}" >&2; exit 1; }
[[ -s "${outdir}/toga.PGA_like.local.v3.assign_candidate_ids.txt" ]] || {
  echo "ERROR: missing candidate IDs for ${assembly}" >&2
  exit 1
}
[[ -s "${outdir}/all.VPS37C-VWCE.locus_annotations.bed" ]] || {
  echo "ERROR: missing locus BED for ${assembly}" >&2
  exit 1
}
mkdir -p "${outdir}"

awk -v assembly="${assembly}" 'NR==1 || index($0, assembly)>0'   "${assignments}" > "${subset}"

python3 "${script_dir}/03_assign_cluster_labels_to_units.py" \
  --assign "${subset}" \
  --ids "${outdir}/toga.PGA_like.local.v3.assign_candidate_ids.txt" \
  --bed "${outdir}/all.VPS37C-VWCE.locus_annotations.bed" \
  -o "${outdir}/toga.PGA_like.local.v4.assign_candidate_ids.bed"

echo "Cluster labels assigned for: ${assembly}"
