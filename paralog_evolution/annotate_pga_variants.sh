#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 5 || $# -gt 6 ]]; then
    echo "Usage: $0 <copy_list> <copies_fasta> <reference_fasta> <reference_gff> <output_dir> [threads]" >&2
    exit 1
fi

COPY_LIST=$1
COPIES_FASTA=$2
REFERENCE_FASTA=$3
REFERENCE_GFF=$4
OUTDIR=$5
THREADS=${6:-16}

mkdir -p "${OUTDIR}"

tmpdir=$(mktemp -d)
trap 'rm -rf "${tmpdir}"' EXIT

while IFS= read -r copy_id || [[ -n "${copy_id}" ]]; do
    [[ -z "${copy_id}" ]] && continue

    safe_id=${copy_id//\//_}

    seqkit grep -p "${copy_id}" "${COPIES_FASTA}" > "${tmpdir}/query.fa"

    nseq=$(grep -c '^>' "${tmpdir}/query.fa" || true)
    if [[ "${nseq}" -ne 1 ]]; then
        echo "ERROR: expected one FASTA record for '${copy_id}', found ${nseq}" >&2
        exit 1
    fi

    minimap2 \
        -x asm5 \
        --cs \
        -c \
        -t "${THREADS}" \
        "${REFERENCE_FASTA}" \
        "${tmpdir}/query.fa" \
        > "${tmpdir}/query.paf"

    paftools.js call \
        -l 0 \
        -L 0 \
        -q 0 \
        -f "${REFERENCE_FASTA}" \
        "${tmpdir}/query.paf" | \
        bgzip -c > "${tmpdir}/query.vcf.gz"

    vep \
        --fasta "${REFERENCE_FASTA}" \
        --gff "${REFERENCE_GFF}" \
        --force_overwrite \
        --input_file "${tmpdir}/query.vcf.gz" \
        --output_file "${tmpdir}/query.vep.txt"

    filter_vep \
        --input_file "${tmpdir}/query.vep.txt" \
        --output_file "${OUTDIR}/${safe_id}.vep.filter.txt" \
        --force_overwrite \
        --filter 'IMPACT != MODIFIER'

done < "${COPY_LIST}"
