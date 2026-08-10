#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 || $# -gt 6 ]]; then
    echo "Usage: $0 <graph.vcf.gz> <GRCh38.fa> <region> <chimpanzee.fa> <output_prefix> [threads]" >&2
    exit 1
fi

VCF=$1
REF=$2
REGION=$3
CHIMP=$4
PREFIX=$5
THREADS=${6:-16}

for cmd in bcftools tabix bgzip vcfcreatemulti samtools mafft trimal iqtree; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "Missing required command: $cmd" >&2; exit 1; }
done

for file in "$VCF" "$REF" "$CHIMP"; do
    [[ -s "$file" ]] || { echo "Input not found or empty: $file" >&2; exit 1; }
done

workdir=$(mktemp -d "${TMPDIR:-/tmp}/pga_flank.XXXXXX")
trap 'rm -rf "$workdir"' EXIT

regional_vcf="$workdir/region.vcf.gz"
regional_multi="$workdir/region.multi.vcf.gz"
regional_ref="$workdir/reference.fa"
human_fa="${PREFIX}.human.fa"
all_fa="${PREFIX}.with_chimpanzee.fa"
aln_fa="${PREFIX}.mafft.fa"
trim_fa="${PREFIX}.mafft.trim.fa"

# Keep variants up to 1 kb in the selected flanking interval.
bcftools view \
    -r "$REGION" \
    -i 'INFO/LEN<=1000' \
    "$VCF" \
    -Oz -o "$regional_vcf"
tabix -f -p vcf "$regional_vcf"

# Merge overlapping indel records before consensus reconstruction.
vcfcreatemulti "$regional_vcf" | bgzip -@ "$THREADS" > "$regional_multi"
tabix -f -p vcf "$regional_multi"

# bcftools consensus expects a reference FASTA matching the selected interval.
samtools faidx "$REF" "$REGION" > "$regional_ref"

: > "$human_fa"
while read -r sample; do
    sample_vcf="$workdir/${sample}.vcf.gz"

    bcftools view -s "$sample" "$regional_multi" \
        | bcftools +fill-tags -- -t AC,AN,AF \
        | bcftools view -a -c 1 -Oz -o "$sample_vcf"
    tabix -f -p vcf "$sample_vcf"

    bcftools consensus -f "$regional_ref" "$sample_vcf" -H 1 \
        | sed "s/^>.*/>${sample}/" \
        >> "$human_fa"
done < <(bcftools query -l "$regional_multi")

cat "$human_fa" "$CHIMP" > "$all_fa"

mafft --auto --thread -1 "$all_fa" > "$aln_fa"
trimal -in "$aln_fa" -out "$trim_fa" -automated1 -keepheader

iqtree \
    -s "$trim_fa" \
    --prefix "$PREFIX" \
    -T auto \
    -B 1000 \
    -bnni \
    -safe \
    -o chimpanzee

echo "Tree: ${PREFIX}.treefile"
