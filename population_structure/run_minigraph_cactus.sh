#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
    echo "Usage: $0 <assemblies.tsv> <GRCh38.regions.fa> <assembly_regions_dir> <output_dir> [threads]" >&2
    exit 1
fi

MANIFEST=$1
REF_REGIONS_FA=$2
ASSEMBLY_REGIONS=$3
OUTDIR=$4
THREADS=${5:-32}
N_REGIONS=20

# External repositories used in the original analysis.
PREPARE_VCF_MC_DIR="/path/to/genotyping-pipelines/prepare-vcf-MC"
COLLAPSE_BUBBLE_DIR="/path/to/collapse-bubble"
SAMPLE_INFO="/path/to/Samples.Sex.tsv"

# Minigraph-Cactus resource settings used in the original analysis.
MG_MEMORY="200G"
CONS_MEMORY="200G"
INDEX_MEMORY="200G"

mkdir -p "${OUTDIR}"

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

for i in $(seq 1 "${N_REGIONS}"); do
    region="region${i}"
    workdir="${OUTDIR}/${region}"
    mkdir -p "${workdir}"

    # Local GRCh38 reference sequence with a contig name matching the VCF.
    seqkit grep -n -p "GRCh38.${region}" "${REF_REGIONS_FA}" \
        | seqkit replace -p '^.+$' -r "${region}" \
        > "${workdir}/GRCh38.${region}.fa"

    seqfile="${workdir}/seq.file"
    : > "${seqfile}"

    while IFS=$'\t' read -r sample hap source fasta; do
        [[ "${sample}" == "sample" ]] && continue
        [[ -n "${sample}" && -n "${hap}" ]] || continue

        graph_id=$(sample_to_graph_id "${sample}" "${hap}")
        region_fa="${ASSEMBLY_REGIONS}/${region}/${graph_id}.fa"
        if [[ -s "${region_fa}" ]]; then
            printf '%s\t%s\n' "${graph_id}" "${region_fa}" >> "${seqfile}"
        fi
    done < "${MANIFEST}"

    printf 'GRCh38\t%s\n' "${workdir}/GRCh38.${region}.fa" >> "${seqfile}"

    cactus-pangenome \
        "${workdir}/js" \
        "${seqfile}" \
        --outDir "${workdir}" \
        --outName "${region}" \
        --reference GRCh38 CHM13v2 CN1v1 \
        --vcf \
        --vcfReference GRCh38 \
        --odgi clip \
        --mapCores "${THREADS}" \
        --mgMemory "${MG_MEMORY}" \
        --consMemory "${CONS_MEMORY}" \
        --indexMemory "${INDEX_MEMORY}"

    # prepare-vcf-MC configuration for this regional graph.
    cat > "${workdir}/prepare-vcf-MC.config.yaml" <<CFG
results: "${workdir}"
callsets:
 ${region}:
  vcf: "${workdir}/${region}.vcf.gz"
  gfa: "${workdir}/${region}.gfa.gz"
  sample_info: "${SAMPLE_INFO}"
  reference_prefix: "chr"
  reference_version: "GRCh38"
CFG

    snakemake \
        -s "${PREPARE_VCF_MC_DIR}/workflow/Snakefile" \
        --configfile "${workdir}/prepare-vcf-MC.config.yaml" \
        -j "${THREADS}" \
        --directory "${workdir}"

    biallelic_vcf="${workdir}/vcf/${region}/${region}_filtered_ids_biallelic.vcf.gz"
    uniqid_vcf="${workdir}/${region}.biallelic.uniqid.vcf.gz"
    wave_vcf="${workdir}/${region}.biallelic.uniqid.vcfwave.vcf.gz"
    norm_vcf="${workdir}/${region}.biallelic.uniqid.vcfwave.norm.vcf.gz"
    merged_vcf="${workdir}/${region}.biallelic.uniqid.vcfwave.norm.merge_dup.vcf.gz"
    final_vcf="${workdir}/${region}.SNP.vcf.gz"

    python "${COLLAPSE_BUBBLE_DIR}/annotate_var_id.py" \
        -i "${biallelic_vcf}" \
        -o "${uniqid_vcf}"

    bcftools annotate -x INFO/AT -Ov "${uniqid_vcf}" \
        | bcftools +fill-tags -Ov -- -t AC,AN,AF \
        | vcfwave -t "${THREADS}" -I 1000 \
        | bgzip -@ "${THREADS}" -c \
        > "${wave_vcf}"

    bcftools norm \
        --threads "${THREADS}" \
        -f "${workdir}/GRCh38.${region}.fa" \
        "${wave_vcf}" \
        | bcftools sort -Oz -o "${norm_vcf}"
    tabix -f -p vcf "${norm_vcf}"

    python "${COLLAPSE_BUBBLE_DIR}/merge_duplicates.py" \
        -i "${norm_vcf}" \
        -o "${merged_vcf}.tmp.vcf.gz" \
        -c repeat \
        -t ID

    bcftools +fill-tags \
        "${merged_vcf}.tmp.vcf.gz" \
        -Oz -o "${merged_vcf}" \
        -- -t AC,AN,AF
    tabix -f -p vcf "${merged_vcf}"

    bcftools view -v snps "${merged_vcf}" \
        | bcftools view -q 0.05:minor -Oz -o "${final_vcf}"
    tabix -f -p vcf "${final_vcf}"
done
