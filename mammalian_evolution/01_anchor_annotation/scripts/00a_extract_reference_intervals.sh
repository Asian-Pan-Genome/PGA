#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <reference_bed_dir> [output_dir]" >&2
  exit 2
fi

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
reorganized_root="$(cd -- "${script_dir}/../.." && pwd)"
bed_dir=$1
output_dir=${2:-"${reorganized_root}/shared_resources/assembly_and_toga/whitelists"}
builder="${script_dir}/00b_build_canonical_interval_whitelist.py"

mkdir -p "${output_dir}"

for ref_prefix in hg38 mm10 HLbosTau10 HLeleMaxInd3A; do
  bed="${bed_dir}/${ref_prefix}.toga.transcripts.bed"
  if [[ ! -s "${bed}" ]]; then
    echo "ERROR: missing reference BED: ${bed}" >&2
    exit 1
  fi

  normalized_prefix=$(printf '%s' "${ref_prefix}" | tr '[:upper:]' '[:lower:]')
  python3 "${builder}" \
    --bed "${bed}" \
    --ref-prefix "${ref_prefix}" \
    -o "${output_dir}/${normalized_prefix}_vps37c_vwce_canonical_interval_whitelist.tsv"
done

echo "Reference interval whitelists written to: ${output_dir}"
