#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 <assembly> <toga_assembly_dir> <output_root> [whitelist.tsv]" >&2
  exit 2
fi

species=$1
toga_dir=$2
output_root=$3
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
reorganized_root="$(cd -- "${script_dir}/../.." && pwd)"
whitelist=${4:-"${reorganized_root}/shared_resources/assembly_and_toga/whitelists/four_ref_pga_like_target_source_whitelist.tsv"}
outdir="${output_root}/${species}"

[[ -s "${outdir}/all.VPS37C-VWCE.locus_annotations.bed" ]] || {
  echo "ERROR: run 01_extract_anchor_locus.sh first for ${species}" >&2
  exit 1
}
[[ -s "${whitelist}" ]] || { echo "ERROR: missing whitelist: ${whitelist}" >&2; exit 1; }
mkdir -p "${outdir}"

python3 "${script_dir}/02_refine_toga_local_units.py" \
  --bed "${outdir}/all.VPS37C-VWCE.locus_annotations.bed" \
  --nucleotide "${toga_dir}/nucleotide.fa.gz" \
  --protein "${toga_dir}/protein.fa.gz" \
  --assembly "${species}" \
  --out-prefix "${outdir}/toga.PGA_like.local.v3" \
  --whitelist-tsv "${whitelist}" \
  --target-exons 9 \
  --min-seed-span 4000 \
  --max-compact-seed-span 15000 \
  --max-seed-span 18000 \
  --partial-min-exons 5 \
  --weak-min-exons 2 \
  --seed-cluster-overlap 0.35 \
  --seed-exonic-overlap 0.50 \
  --seed-shared-exons 3 \
  --attach-overlap 0.35 \
  --attach-exonic-overlap 0.35 \
  --attach-shared-exons 2 \
  --attach-min-bp 500 \
  --fused-hard-span 35000 \
  --fused-fold 2.0 \
  --stretch-fold 1.35 \
  --stretch-exon-ov 0.30 \
  --stretch-exon-ov-unit 0.50 \
  --stretch-shared-exons 3 \
  --prefer-id-prefix hg38.ENST \
  --rearrange-min-non-target-genes 3 \
  --rearrange-min-non-target-frac 0.30

grep -w -f "${outdir}/toga.PGA_like.local.v3.complete_unit_representative_ids.txt" \
  "${outdir}/all.VPS37C-VWCE.locus_annotations.bed" \
  > "${outdir}/toga.PGA_like.local.v3.complete_unit_representative_ids.bed"

echo "Local-unit refinement completed: ${outdir}"
