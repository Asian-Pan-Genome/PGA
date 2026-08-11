#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: $0 <1KGP.vcf.gz> <sample_metadata.tsv> <output_dir> [threads]"
    exit 1
fi

VCF=$1
METADATA=$2
OUTDIR=$3
THREADS=${4:-16}

POPULATIONS=(AFR AMR EAS EUR SAS CSA)

mkdir -p "${OUTDIR}"

for pop in "${POPULATIONS[@]}"; do
    echo "[INFO] Preparing ${pop}"

    SAMPLE_LIST="${OUTDIR}/${pop}.samples.txt"
    PHENO_FILE="${OUTDIR}/pheno.${pop}.txt"
    POP_VCF="${OUTDIR}/1KG.${pop}.vcf.gz"

    # AFR combines AFR-W and AFR-E&S; other superpopulations are matched exactly.
    awk -F '\t' -v pop="${pop}" '
        NR > 1 && $2 != "" && $3 != "" &&
        ($3 == pop || (pop == "AFR" && $3 ~ /^AFR-/)) {
            print $1
        }
    ' "${METADATA}" > "${SAMPLE_LIST}"

    # PLINK phenotype format: FID, IID, phenotype.
    awk -F '\t' -v pop="${pop}" '
        BEGIN {OFS="\t"}
        NR > 1 && $2 != "" && $3 != "" &&
        ($3 == pop || (pop == "AFR" && $3 ~ /^AFR-/)) {
            print $1, $1, $2
        }
    ' "${METADATA}" > "${PHENO_FILE}"

    if [[ ! -s "${SAMPLE_LIST}" ]]; then
        echo "[ERROR] No samples found for ${pop}" >&2
        exit 1
    fi

    bcftools view \
        -S "${SAMPLE_LIST}" \
        "${VCF}" \
        --threads "${THREADS}" \
        -Ou | \
    bcftools +fill-tags \
        - \
        --threads "${THREADS}" \
        -Ou \
        -- -t AC,AF,AN | \
    bcftools view \
        -a \
        -c 1 \
        --threads "${THREADS}" \
        -Ou | \
    bcftools view \
        -q 0.05:minor \
        --threads "${THREADS}" \
        -Oz \
        -o "${POP_VCF}"

    bcftools index -f -t "${POP_VCF}"

    plink \
        --vcf "${POP_VCF}" \
        --indep-pairwise 50 5 0.2 \
        --threads "${THREADS}" \
        --out "${OUTDIR}/${pop}.pruning"

    plink \
        --vcf "${POP_VCF}" \
        --extract "${OUTDIR}/${pop}.pruning.prune.in" \
        --pca 10 \
        --threads "${THREADS}" \
        --out "${OUTDIR}/${pop}.Global_PCA"
done
