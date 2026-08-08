#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
    echo "Usage: $0 <manifest.tsv> <reference.PGA.fa> <reference.PGA.gff_db> <output_dir> [threads]" >&2
    exit 1
fi

MANIFEST=$(readlink -f "$1")
REF_PGA_FASTA=$(readlink -f "$2")
REF_PGA_GFF_DB=$(readlink -f "$3")
OUTDIR=$4
THREADS=${5:-16}

FLANK_BP=100000
MIN_GENE_BP=8000
MAX_GENE_BP=12000

mkdir -p "${OUTDIR}"
OUTDIR=$(cd "${OUTDIR}" && pwd)

while IFS=$'\t' read -r sample_id species fasta toga_gtf; do
    [[ -z "${sample_id}" ]] && continue
    [[ "${sample_id}" == "sample_id" ]] && continue

    fasta=$(readlink -f "${fasta}")
    toga_gtf=$(readlink -f "${toga_gtf}")

    echo "[INFO] Processing ${sample_id} (${species})"

    (
        workdir="${OUTDIR}/${sample_id}"
        mkdir -p "${workdir}"
        cd "${workdir}"

        # Extract candidate TOGA PGA annotations.
        gzip -cd "${toga_gtf}" \
            | grep -E 'PGA.\.' \
            > PGA.TOGA.gtf

        agat_convert_sp_gxf2gxf.pl \
            -g PGA.TOGA.gtf \
            -o PGA.TOGA.gff

        awk -F'\t' '$3 != "transcript"' PGA.TOGA.gff \
            > PGA.TOGA.no_transcript.gff

        # Retain full-length TOGA PGA gene models.
        awk -F'\t' \
            -v min_bp="${MIN_GENE_BP}" \
            -v max_bp="${MAX_GENE_BP}" \
            '
            $3 == "gene" {
                gene_len = $5 - $4 + 1
                if (gene_len >= min_bp && gene_len <= max_bp) {
                    split($9, a, ";")
                    sub(/^ID=/, "", a[1])
                    print a[1]
                }
            }
            ' PGA.TOGA.no_transcript.gff \
            > PGA.TOGA.keep_ids.txt

        if [[ ! -s PGA.TOGA.keep_ids.txt ]]; then
            echo "[ERROR] No 8-12 kb TOGA PGA gene models retained for ${sample_id}." >&2
            exit 1
        fi

        agat_sp_filter_feature_from_keep_list.pl \
            -f PGA.TOGA.no_transcript.gff \
            --kl PGA.TOGA.keep_ids.txt \
            -o PGA.TOGA.filtered.gff

        agat_sp_fix_features_locations_duplicated.pl \
            -g PGA.TOGA.filtered.gff \
            -o PGA.TOGA.curated.gff

        # Define and extract the local PGA interval (+/- 100 kb).
        n_contigs=$(awk -F'\t' '$3 == "gene" {print $1}' PGA.TOGA.curated.gff | sort -u | wc -l)
        if [[ "${n_contigs}" -ne 1 ]]; then
            echo "[ERROR] Expected retained PGA genes on one contig for ${sample_id}; found ${n_contigs}." >&2
            exit 1
        fi

        chrom=$(awk -F'\t' '$3 == "gene" {print $1; exit}' PGA.TOGA.curated.gff)
        gene_start=$(awk -F'\t' '$3 == "gene" {if (min == "" || $4 < min) min = $4} END {print min}' PGA.TOGA.curated.gff)
        gene_end=$(awk -F'\t' '$3 == "gene" {if (max == "" || $5 > max) max = $5} END {print max}' PGA.TOGA.curated.gff)

        interval_start=$((gene_start - 1 - FLANK_BP))
        interval_end=$((gene_end + FLANK_BP))
        (( interval_start < 0 )) && interval_start=0

        printf '%s\t%s\t%s\n' "${chrom}" "${interval_start}" "${interval_end}" > PGA.interval.bed

        bedtools getfasta \
            -fi "${fasta}" \
            -bed PGA.interval.bed \
            > "${sample_id}.PGA.fa"

        local_header="${chrom}:${interval_start}-${interval_end}"

        # Convert genomic GFF coordinates to coordinates on the extracted local sequence.
        awk -F'\t' \
            -v OFS='\t' \
            -v header="${local_header}" \
            -v start="${interval_start}" \
            '
            !/^#/ {
                $1 = header
                $4 -= start
                $5 -= start
                if ($2 == "stdin") $2 = "TOGA"
                print
            }
            ' PGA.TOGA.curated.gff \
            > PGA.TOGA.local.gff

        # Mask retained TOGA copies and rescue additional PGA copies with Liftoff.
        bedtools maskfasta \
            -fi "${sample_id}.PGA.fa" \
            -bed PGA.TOGA.local.gff \
            -fo "${sample_id}.PGA.hard_mask.fa"

        liftoff \
            "${sample_id}.PGA.hard_mask.fa" \
            "${REF_PGA_FASTA}" \
            -sc 0.90 \
            -copies \
            -db "${REF_PGA_GFF_DB}" \
            -polish \
            -exclude_partial \
            -p "${THREADS}" \
            -dir liftoff_intermediate \
            -o "${sample_id}.liftoff_rescue.gff"

        # Standardize rescued annotations and merge them with TOGA annotations.
        agat_sp_manage_IDs.pl \
            --gff "${sample_id}.liftoff_rescue.gff_polished" \
            -o "${sample_id}.liftoff_rescue.fixed.gff"

        cat \
            PGA.TOGA.local.gff \
            <(grep -v '^#' "${sample_id}.liftoff_rescue.fixed.gff") \
            > "${sample_id}.TOGA_liftoff.unsorted.gff"

        agat_convert_sp_gxf2gxf.pl \
            -g "${sample_id}.TOGA_liftoff.unsorted.gff" \
            -o "${sample_id}.TOGA_liftoff.gff"

        echo "[INFO] Final annotation: ${workdir}/${sample_id}.TOGA_liftoff.gff"
    )
done < "${MANIFEST}"
