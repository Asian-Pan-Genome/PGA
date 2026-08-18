#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
  echo "Usage: $0 <protein.fa> [output_dir] [threads]" >&2
  exit 2
fi

seq=$1
output_dir=${2:-"$PWD"}
threads=${3:-16}
name=$(basename "${seq}")
stem=${name%.*}
alignment="${output_dir}/${stem}.aln.fa"
trimmed="${output_dir}/${stem}.trim.fa"

[[ -s "${seq}" ]] || { echo "ERROR: missing protein FASTA: ${seq}" >&2; exit 1; }
for command_name in mafft trimal iqtree2; do
  command -v "${command_name}" >/dev/null 2>&1 || {
    echo "ERROR: required command is not on PATH: ${command_name}" >&2
    exit 1
  }
done
mkdir -p "${output_dir}"

mafft --auto --reorder --thread "${threads}" "${seq}" > "${alignment}"
trimal -in "${alignment}" -out "${trimmed}" -automated1
iqtree2 -s "${trimmed}" -m MFP -bb 1000 -alrt 1000 -T AUTO

echo "Protein alignment and IQ-TREE analysis completed: ${trimmed}"
