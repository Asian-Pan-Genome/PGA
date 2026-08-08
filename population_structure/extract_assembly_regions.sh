#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: $0 <assemblies.tsv> <GRCh38.regions.fa> <output_dir> [threads]" >&2
    exit 1
fi

MANIFEST=$1
REF_REGIONS_FA=$2
OUTDIR=$3
THREADS=${4:-8}
N_REGIONS=20

mkdir -p "${OUTDIR}"
for i in $(seq 1 "${N_REGIONS}"); do
    mkdir -p "${OUTDIR}/region${i}"
    : > "${OUTDIR}/region${i}/failed.list"
done

sample_to_graph_id() {
    local sample=$1
    local hap=$2
    case "${hap}" in
        hap1) printf '%s.1' "${sample}" ;;
        hap2) printf '%s.2' "${sample}" ;;
        hap0) printf '%s' "${sample}" ;;
        *)
            echo "Unsupported hap value: ${hap}" >&2
            return 1
            ;;
    esac
}

while IFS=$'\t' read -r sample hap source fasta; do
    [[ "${sample}" == "sample" ]] && continue
    [[ -n "${sample}" && -n "${hap}" && -n "${fasta}" ]] || continue
    [[ -s "${fasta}" ]] || { echo "Missing FASTA: ${fasta}" >&2; exit 1; }

    graph_id=$(sample_to_graph_id "${sample}" "${hap}")
    paf="${OUTDIR}/${sample}.${hap}.paf"

    minimap2 \
        -x asm5 \
        -c \
        --secondary=no \
        -t "${THREADS}" \
        <(seqkit grep -r -p 'chr11' "${fasta}") \
        "${REF_REGIONS_FA}" \
        > "${paf}"

    for i in $(seq 1 "${N_REGIONS}"); do
        region="region${i}"
        query_name="GRCh38.${region}"
        tmp_bed=$(mktemp)

        paftools.js liftover -q 0 "${paf}" \
            <(printf '%s\t0\t1000000\n' "${query_name}") \
            | bedtools sort -i - \
            > "${tmp_bed}"

        if [[ ! -s "${tmp_bed}" ]]; then
            echo "${graph_id}" >> "${OUTDIR}/${region}/failed.list"
            rm -f "${tmp_bed}"
            continue
        fi

        chrom=$(awk 'NR==1 {print $1}' "${tmp_bed}")
        start=$(awk 'NR==1 {print $2}' "${tmp_bed}")
        end=$(awk 'END {print $3}' "${tmp_bed}")
        span_bed=$(mktemp)
        printf '%s\t%s\t%s\n' "${chrom}" "${start}" "${end}" > "${span_bed}"

        out_fa="${OUTDIR}/${region}/${graph_id}.fa"
        bedtools getfasta -fi "${fasta}" -bed "${span_bed}" \
            | sed "1s/^>.*/>${graph_id}.${region}/" \
            > "${out_fa}"

        if [[ ! -s "${out_fa}" ]]; then
            echo "${graph_id}" >> "${OUTDIR}/${region}/failed.list"
        fi

        rm -f "${tmp_bed}" "${span_bed}"
    done
done < "${MANIFEST}"
