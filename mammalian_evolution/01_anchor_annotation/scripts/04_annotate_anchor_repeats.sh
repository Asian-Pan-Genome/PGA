#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 5 ]]; then
  cat >&2 <<'USAGE'
Usage:
  04_annotate_anchor_repeats.sh <assembly> <output_root> [threads] \
    [trf_dat2gff3.py] [RepeatMasker_command]

The TRF converter is optional because it is retained locally under review_hold
and is not a GitHub candidate. Without it, the TRF .dat output is still kept.
USAGE
  exit 2
fi

species=$1
output_root=$2
threads=${3:-8}
trf_converter=${4:-}
repeatmasker_cmd=${5:-RepeatMasker}
outdir="${output_root}/${species}"
seq="${outdir}/pga.anchor.locus.fa"
renamed_seq="${outdir}/${species}.pga.anchor.locus.fa"

[[ -s "${seq}" ]] || { echo "ERROR: missing anchor sequence: ${seq}" >&2; exit 1; }
command -v "${repeatmasker_cmd}" >/dev/null 2>&1 || {
  echo "ERROR: RepeatMasker command not found: ${repeatmasker_cmd}" >&2
  exit 1
}
command -v trf >/dev/null 2>&1 || { echo "ERROR: trf is not on PATH" >&2; exit 1; }

sed "s/${species}/pga_anchor_locus/g" "${seq}" > "${renamed_seq}"

"${repeatmasker_cmd}"   -engine rmblast   -species mammalia   -pa "${threads}"   -gff   -xsmall   -dir "${outdir}"   "${renamed_seq}"

(
  cd "${outdir}"
  trf "$(basename "${renamed_seq}")" 2 7 7 80 10 50 500 -d -h
)

raw_trf="${renamed_seq}.2.7.7.80.10.50.500.dat"
final_trf="${outdir}/${species}.pga.anchor.locus.trf.dat"
[[ -s "${raw_trf}" ]] || { echo "ERROR: expected TRF output not found: ${raw_trf}" >&2; exit 1; }
mv "${raw_trf}" "${final_trf}"

if [[ -n "${trf_converter}" ]]; then
  [[ -s "${trf_converter}" ]] || {
    echo "ERROR: missing TRF converter: ${trf_converter}" >&2
    exit 1
  }
  python3 "${trf_converter}" "${final_trf}"     "${outdir}/${species}.pga.anchor.locus.trf.gff"
else
  echo "INFO: no TRF converter supplied; retaining the TRF .dat file only"
fi

echo "Repeat annotation completed: ${outdir}"
