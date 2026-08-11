#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 <population_panel_dir> <output_dir> [threads]"
    exit 1
fi

PANEL_DIR=$1
OUTDIR=$2
THREADS=${3:-16}

POPULATIONS=(AFR AMR EAS EUR SAS CSA)

mkdir -p "${OUTDIR}"

for chrom in $(seq 1 22); do
    CHROM="chr${chrom}"
    CHROM_DIR="${OUTDIR}/${CHROM}"
    mkdir -p "${CHROM_DIR}"

    for pop in "${POPULATIONS[@]}"; do
        echo "[INFO] Testing ${pop} ${CHROM}"

        POP_VCF="${PANEL_DIR}/1KG.${pop}.vcf.gz"
        PCA_FILE="${PANEL_DIR}/${pop}.Global_PCA.eigenvec"
        PHENO_FILE="${PANEL_DIR}/pheno.${pop}.txt"

        CHR_VCF="${CHROM_DIR}/1KG.${pop}.${CHROM}.vcf.gz"
        RAW_PREFIX="${CHROM_DIR}/${pop}.${CHROM}.assoc_raw"
        ASSOC_PREFIX="${CHROM_DIR}/${pop}.${CHROM}.assoc_Final_Result"

        bcftools view \
            -r "${CHROM}" \
            "${POP_VCF}" \
            --threads "${THREADS}" \
            -Oz \
            -o "${CHR_VCF}"

        plink \
            --vcf "${CHR_VCF}" \
            --make-bed \
            --keep-allele-order \
            --double-id \
            --threads "${THREADS}" \
            --out "${RAW_PREFIX}"

        plink \
            --bfile "${RAW_PREFIX}" \
            --pheno "${PHENO_FILE}" \
            --linear \
            --covar "${PCA_FILE}" \
            --hide-covar \
            --allow-no-sex \
            --threads "${THREADS}" \
            --out "${ASSOC_PREFIX}"
    done
done
