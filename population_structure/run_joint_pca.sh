#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 9 || $# -gt 10 ]]; then
    cat >&2 <<USAGE
Usage: $0 <assemblies.tsv> <hg38.regions.bed> <GRCh38.regions.fa> <graph_dir> \
          <1KGP.chr11.vcf.gz> <PGA.copies.tsv> <assemblies.new_superpop.list> \
          <1KGP.metadata.txt> <output_dir> [threads]
USAGE
    exit 1
fi

MANIFEST=$1
REGIONS_BED=$2
REF_REGIONS_FA=$3
GRAPH_DIR=$4
KG_CHR11_VCF=$5
ASSEMBLY_COPIES=$6
ASSEMBLY_SUPERPOP=$7
KG_METADATA=$8
OUTDIR=$9
THREADS=${10:-32}
N_REGIONS=20

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

mkdir -p \
    "${OUTDIR}/assembly" \
    "${OUTDIR}/kg" \
    "${OUTDIR}/ref" \
    "${OUTDIR}/shared" \
    "${OUTDIR}/concordance_filter" \
    "${OUTDIR}/plink" \
    "${OUTDIR}/logs"

# -----------------------------------------------------------------------------
# 1. Define the diploid assembly sample set from paired haplotypes.
# -----------------------------------------------------------------------------
MASTER_HAPS="${OUTDIR}/assembly/master_haps.txt"
MASTER_DIPLOID="${OUTDIR}/assembly/master_diploid.samples.txt"

awk -F'\t' 'NR > 1 {
    if ($2 == "hap1") print $1 ".1";
    else if ($2 == "hap2") print $1 ".2";
    else if ($2 == "hap0") print $1;
}' "${MANIFEST}" > "${MASTER_HAPS}"

python "${SCRIPT_DIR}/make_master_diploid_samples.py" \
    --master-haps "${MASTER_HAPS}" \
    --out "${MASTER_DIPLOID}"

# -----------------------------------------------------------------------------
# 2. Standardize the 20 assembly-derived regional VCFs.
# -----------------------------------------------------------------------------
assembly_vcfs=()
for i in $(seq 1 "${N_REGIONS}"); do
    region="region${i}"
    in_vcf="${GRAPH_DIR}/${region}/${region}.SNP.vcf.gz"
    failed="${GRAPH_DIR}/${region}/failed.list"
    tmp_vcf="${OUTDIR}/assembly/${region}.standardized.vcf"
    out_vcf="${tmp_vcf}.gz"

    [[ -s "${in_vcf}" ]] || { echo "Missing assembly VCF: ${in_vcf}" >&2; exit 1; }

    python "${SCRIPT_DIR}/standardize_assembly_region_diploid.py" \
        --vcf "${in_vcf}" \
        --samples "${MASTER_DIPLOID}" \
        --failed-list "${failed}" \
        --expected-chrom "${region}" \
        --out "${tmp_vcf}" \
        > "${OUTDIR}/logs/${region}.standardize.log" 2>&1

    bgzip -@ "${THREADS}" -f "${tmp_vcf}"
    tabix -f -p vcf "${out_vcf}"
    assembly_vcfs+=("${out_vcf}")
done

bcftools concat \
    --threads "${THREADS}" \
    -Oz \
    -o "${OUTDIR}/assembly/assembly.all20.standardized.vcf.gz" \
    "${assembly_vcfs[@]}"
tabix -f -p vcf "${OUTDIR}/assembly/assembly.all20.standardized.vcf.gz"

# -----------------------------------------------------------------------------
# 3. Prepare local GRCh38 reference FASTAs and the matching 1KGP SNPs.
# -----------------------------------------------------------------------------
kg_region_vcfs=()
for i in $(seq 1 "${N_REGIONS}"); do
    region="region${i}"

    seqkit grep -n -p "GRCh38.${region}" "${REF_REGIONS_FA}" \
        | seqkit replace -p '^.+$' -r "${region}" \
        > "${OUTDIR}/ref/${region}.fa"
    samtools faidx "${OUTDIR}/ref/${region}.fa"

    bed_line=$(awk -v n="${i}" 'NR == n {print; exit}' "${REGIONS_BED}")
    [[ -n "${bed_line}" ]] || { echo "Missing region ${i} in ${REGIONS_BED}" >&2; exit 1; }

    chrom=$(printf '%s\n' "${bed_line}" | cut -f1)
    start0=$(printf '%s\n' "${bed_line}" | cut -f2)
    end0=$(printf '%s\n' "${bed_line}" | cut -f3)

    [[ "${chrom}" == "chr11" ]] || { echo "region${i} is not on chr11: ${chrom}" >&2; exit 1; }
    start1=$((start0 + 1))

    raw_vcf="${OUTDIR}/kg/${region}.raw.vcf.gz"
    norm_vcf="${OUTDIR}/kg/${region}.norm.vcf.gz"

    {
        bcftools view \
            -r "${chrom}:${start1}-${end0}" \
            -v snps -m2 -M2 \
            --threads "${THREADS}" \
            -h "${KG_CHR11_VCF}" \
        | awk -v c="${region}" '
            /^##contig=/ {next}
            /^#CHROM/ {
                print "##contig=<ID=" c ",length=1000000>"
                print
                next
            }
            {print}
        '

        bcftools view \
            -r "${chrom}:${start1}-${end0}" \
            -v snps -m2 -M2 \
            --threads "${THREADS}" \
            -H "${KG_CHR11_VCF}" \
        | awk -v c="${region}" -v s="${start0}" 'BEGIN {OFS="\t"} {
            $1 = c
            $2 = $2 - s
            print
        }'
    } | bcftools sort -Oz -o "${raw_vcf}" --temp-dir "${OUTDIR}/kg/tmp.${region}"
    tabix -f -p vcf "${raw_vcf}"

    bcftools norm \
        -f "${OUTDIR}/ref/${region}.fa" \
        -c w \
        --threads "${THREADS}" \
        -Ou "${raw_vcf}" \
    | bcftools annotate \
        --set-id '%CHROM:%POS:%REF:%ALT' \
        --threads "${THREADS}" \
        -Oz -o "${norm_vcf}"
    tabix -f -p vcf "${norm_vcf}"

    kg_region_vcfs+=("${norm_vcf}")
done

bcftools concat \
    --threads "${THREADS}" \
    -Oz \
    -o "${OUTDIR}/kg/kg.all20.norm.vcf.gz" \
    "${kg_region_vcfs[@]}"
tabix -f -p vcf "${OUTDIR}/kg/kg.all20.norm.vcf.gz"

# -----------------------------------------------------------------------------
# 4. Retain exact shared SNPs (CHROM/POS/REF/ALT).
# -----------------------------------------------------------------------------
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    "${OUTDIR}/assembly/assembly.all20.standardized.vcf.gz" \
    | sort -k1,1V -k2,2n -k3,3 -k4,4 -u \
    > "${OUTDIR}/shared/assembly.keys.tsv"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    "${OUTDIR}/kg/kg.all20.norm.vcf.gz" \
    | sort -k1,1V -k2,2n -k3,3 -k4,4 -u \
    > "${OUTDIR}/shared/kg.keys.tsv"

awk 'NR == FNR {
        key=$1"\t"$2"\t"$3"\t"$4
        a[key]=1
        next
     }
     {
        key=$1"\t"$2"\t"$3"\t"$4
        if (key in a) print
     }' \
    "${OUTDIR}/shared/assembly.keys.tsv" \
    "${OUTDIR}/shared/kg.keys.tsv" \
    | sort -k1,1V -k2,2n -k3,3 -k4,4 -u \
    > "${OUTDIR}/shared/shared.keys.tsv"

n_shared=$(wc -l < "${OUTDIR}/shared/shared.keys.tsv")
echo "Shared exact SNPs: ${n_shared}" >&2
[[ "${n_shared}" -gt 0 ]] || { echo "No exact shared SNPs found" >&2; exit 1; }

python "${SCRIPT_DIR}/filter_vcf_by_sitekey.py" \
    --vcf "${OUTDIR}/assembly/assembly.all20.standardized.vcf.gz" \
    --sites "${OUTDIR}/shared/shared.keys.tsv" \
    --out "${OUTDIR}/shared/assembly.shared.vcf"
bgzip -@ "${THREADS}" -f "${OUTDIR}/shared/assembly.shared.vcf"
tabix -f -p vcf "${OUTDIR}/shared/assembly.shared.vcf.gz"

python "${SCRIPT_DIR}/filter_vcf_by_sitekey.py" \
    --vcf "${OUTDIR}/kg/kg.all20.norm.vcf.gz" \
    --sites "${OUTDIR}/shared/shared.keys.tsv" \
    --out "${OUTDIR}/shared/kg.shared.vcf"
bgzip -@ "${THREADS}" -f "${OUTDIR}/shared/kg.shared.vcf"
tabix -f -p vcf "${OUTDIR}/shared/kg.shared.vcf.gz"

# -----------------------------------------------------------------------------
# 5. Prefix sample IDs and merge assembly and 1KGP panels.
# -----------------------------------------------------------------------------
bcftools query -l "${OUTDIR}/shared/assembly.shared.vcf.gz" \
    | awk '{print $1"\tASM_"$1}' \
    > "${OUTDIR}/shared/assembly.rename.map"

bcftools reheader \
    -s "${OUTDIR}/shared/assembly.rename.map" \
    -o "${OUTDIR}/shared/assembly.shared.renamed.vcf.gz" \
    "${OUTDIR}/shared/assembly.shared.vcf.gz"
tabix -f -p vcf "${OUTDIR}/shared/assembly.shared.renamed.vcf.gz"

bcftools query -l "${OUTDIR}/shared/kg.shared.vcf.gz" \
    | awk '{print $1"\tKG_"$1}' \
    > "${OUTDIR}/shared/kg.rename.map"

bcftools reheader \
    -s "${OUTDIR}/shared/kg.rename.map" \
    -o "${OUTDIR}/shared/kg.shared.renamed.vcf.gz" \
    "${OUTDIR}/shared/kg.shared.vcf.gz"
tabix -f -p vcf "${OUTDIR}/shared/kg.shared.renamed.vcf.gz"

JOINT_VCF="${OUTDIR}/joint.sharedSNP.diploid.vcf.gz"
bcftools merge \
    -m none \
    --threads "${THREADS}" \
    -Oz \
    -o "${JOINT_VCF}" \
    "${OUTDIR}/shared/assembly.shared.renamed.vcf.gz" \
    "${OUTDIR}/shared/kg.shared.renamed.vcf.gz"
tabix -f -p vcf "${JOINT_VCF}"

# -----------------------------------------------------------------------------
# 6. Build metadata and filter panel-biased sites using duplicated individuals.
# -----------------------------------------------------------------------------
bcftools query -l "${JOINT_VCF}" > "${OUTDIR}/joint.samples.txt"

python "${SCRIPT_DIR}/make_joint_pca_metadata.py" \
    --joint-vcf-samples "${OUTDIR}/joint.samples.txt" \
    --assembly-copies "${ASSEMBLY_COPIES}" \
    --assembly-new-superpop "${ASSEMBLY_SUPERPOP}" \
    --kg-metadata "${KG_METADATA}" \
    --out "${OUTDIR}/joint_pca.metadata.tsv"

FILTER_PREFIX="${OUTDIR}/concordance_filter/joint.sharedSNP.diploid.concord98.missdiff005"
FILTER_VCF="${FILTER_PREFIX}.vcf"

python "${SCRIPT_DIR}/filter_panel_biased_sites_by_overlap.py" \
    --vcf "${JOINT_VCF}" \
    --out-prefix "${FILTER_PREFIX}" \
    --out-vcf "${FILTER_VCF}" \
    --min-overlap-called 10 \
    --min-concordance 0.98 \
    --max-missing-diff 0.05 \
    --max-panel-missing 0.20 \
    > "${OUTDIR}/concordance_filter/filter.log" 2>&1

cat "${OUTDIR}/concordance_filter/filter.log" >&2
bgzip -@ "${THREADS}" -f "${FILTER_VCF}"
tabix -f -p vcf "${FILTER_VCF}.gz"

# -----------------------------------------------------------------------------
# 7. Final joint-panel filtering, LD pruning, and PCA.
# -----------------------------------------------------------------------------
PLINK_PREFIX="${OUTDIR}/plink/joint"

plink \
    --vcf "${FILTER_VCF}.gz" \
    --make-bed \
    --double-id \
    --allow-extra-chr \
    --threads "${THREADS}" \
    --out "${PLINK_PREFIX}.raw"

plink \
    --bfile "${PLINK_PREFIX}.raw" \
    --maf 0.05 \
    --geno 0.05 \
    --make-bed \
    --allow-extra-chr \
    --threads "${THREADS}" \
    --out "${PLINK_PREFIX}.filtered"

plink \
    --bfile "${PLINK_PREFIX}.filtered" \
    --indep-pairwise 50 5 0.2 \
    --allow-extra-chr \
    --threads "${THREADS}" \
    --out "${PLINK_PREFIX}.pruning"

plink \
    --bfile "${PLINK_PREFIX}.filtered" \
    --extract "${PLINK_PREFIX}.pruning.prune.in" \
    --pca 10 \
    --allow-extra-chr \
    --threads "${THREADS}" \
    --out "${PLINK_PREFIX}.PCA"
