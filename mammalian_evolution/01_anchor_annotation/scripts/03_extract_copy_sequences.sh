#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 <assembly> <toga_assembly_dir> <output_root> [output_prefix]" >&2
  exit 2
fi

species=$1
toga_dir=$2
output_root=$3
prefix=${4:-toga.PGA_like.local.v3.assign_candidate_ids}
outdir="${output_root}/${species}"
protein="${toga_dir}/protein.fa.gz"
id_list="${outdir}/toga.PGA_like.local.v3.assign_candidate_ids.txt"

[[ -s "${protein}" ]] || { echo "ERROR: missing protein FASTA: ${protein}" >&2; exit 1; }
[[ -s "${id_list}" ]] || { echo "ERROR: missing candidate ID list: ${id_list}" >&2; exit 1; }
mkdir -p "${outdir}"

seqkit grep -w 0 -f "${id_list}" "${protein}" |
  seqkit replace -p '(.+)' -r "\$1#${species}"   > "${outdir}/${prefix}.protein.fa"

echo "Protein sequences: $(grep -c '^>' "${outdir}/${prefix}.protein.fa" || true)"
echo "Copy-sequence extraction completed: ${outdir}/${prefix}.protein.fa"
