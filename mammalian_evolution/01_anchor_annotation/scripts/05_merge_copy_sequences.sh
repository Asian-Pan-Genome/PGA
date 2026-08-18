#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <per_assembly_output_root> [merged_output_dir]" >&2
  exit 2
fi

base=$1
output_dir=${2:-"$PWD"}
mkdir -p "${output_dir}"

mapfile -t nucleotide_files < <(find "${base}" -type f -name "PGA_like.copy.nucleotide.fa" | sort)
mapfile -t protein_files < <(find "${base}" -type f -name "PGA_like.copy.protein.fa" | sort)
mapfile -t cluster_files < <(find "${base}" -type f -name "toga.PGA_like.copy.clusters.tsv" | sort)

[[ ${#nucleotide_files[@]} -gt 0 ]] || { echo "ERROR: no nucleotide FASTA files found" >&2; exit 1; }
[[ ${#protein_files[@]} -gt 0 ]] || { echo "ERROR: no protein FASTA files found" >&2; exit 1; }
[[ ${#cluster_files[@]} -gt 0 ]] || { echo "ERROR: no cluster tables found" >&2; exit 1; }

cat "${nucleotide_files[@]}" > "${output_dir}/all.PGA_like.copy.nucleotide.fa"
cat "${protein_files[@]}" > "${output_dir}/all.PGA_like.copy.protein.fa"

head -n 1 "${cluster_files[0]}" > "${output_dir}/all.PGA_like.copy.clusters.tsv"
for file in "${cluster_files[@]}"; do
  tail -n +2 "${file}"
done >> "${output_dir}/all.PGA_like.copy.clusters.tsv"

echo "Merged outputs written to: ${output_dir}"
