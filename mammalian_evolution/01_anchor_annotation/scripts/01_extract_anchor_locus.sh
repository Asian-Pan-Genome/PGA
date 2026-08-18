#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 || $# -gt 7 ]]; then
  cat >&2 <<'USAGE'
Usage:
  01_extract_anchor_locus.sh <assembly> <toga_assembly_dir> <genome.2bit> <output_root> \
    [gepard.jar] [gepard_matrix] [legacy_cluster_script.py]

The optional Gepard files enable the locus self-alignment QC image. The optional
legacy clustering script enables the historical transcript-cluster/protein
export; downstream local-unit refinement is otherwise handled by step 02.
USAGE
  exit 2
fi

species=$1
toga_dir=$2
genome=$3
output_root=$4
gepard_jar=${5:-}
gepard_matrix=${6:-}
legacy_cluster_script=${7:-}

bed="${toga_dir}/query_annotation.bed"
protein="${toga_dir}/protein.fa.gz"
outdir="${output_root}/${species}"

[[ -s "${bed}" ]] || { echo "ERROR: missing TOGA BED: ${bed}" >&2; exit 1; }
[[ -s "${genome}" ]] || { echo "ERROR: missing 2bit genome: ${genome}" >&2; exit 1; }
mkdir -p "${outdir}"

awk -F'\t' '{split($4,a,"#"); if(a[2]=="VPS37C") print}' "${bed}" > "${outdir}/VPS37C.bed"
awk -F'\t' '{split($4,a,"#"); if(a[2]=="VWCE") print}' "${bed}" > "${outdir}/VWCE.bed"

mapfile -t vps_contigs < <(cut -f1 "${outdir}/VPS37C.bed" | sort -u)
mapfile -t vwce_contigs < <(cut -f1 "${outdir}/VWCE.bed" | sort -u)
if [[ ${#vps_contigs[@]} -ne 1 || ${#vwce_contigs[@]} -ne 1 ]]; then
  echo "ERROR: each anchor must map to exactly one contig" >&2
  exit 1
fi
contig1=${vps_contigs[0]}
contig2=${vwce_contigs[0]}
if [[ "${contig1}" != "${contig2}" ]]; then
  echo "ERROR: anchors are on different contigs" >&2
  exit 1
fi

start1=$(cut -f2 "${outdir}/VPS37C.bed" | sort -n | head -1)
end1=$(cut -f3 "${outdir}/VPS37C.bed" | sort -n | tail -1)
start2=$(cut -f2 "${outdir}/VWCE.bed" | sort -n | head -1)
end2=$(cut -f3 "${outdir}/VWCE.bed" | sort -n | tail -1)

intergenic_start_pos=$(printf "%d\n%d\n" "${end1}" "${end2}" | sort -n | head -1)
intergenic_end_pos=$(printf "%d\n%d\n" "${start1}" "${start2}" | sort -n | tail -1)
if [[ ${intergenic_start_pos} -ge ${intergenic_end_pos} ]]; then
  echo "ERROR: invalid internal anchor interval" >&2
  exit 1
fi

echo "Interval: ${contig1}:${intergenic_start_pos}-${intergenic_end_pos}"

awk -v chr="${contig1}" -v s="${intergenic_start_pos}" -v e="${intergenic_end_pos}" '
$1==chr && $2>=s && $3<=e
' "${bed}" > "${outdir}/all.VPS37C-VWCE.locus_annotations.bed"

awk '$10 == 9' "${outdir}/all.VPS37C-VWCE.locus_annotations.bed"   > "${outdir}/all.VPS37C-VWCE.locus_annotations.full_exon.bed"

printf "%s\t%s\t%s\t%s\n" "${species}" "${contig1}"   "${intergenic_start_pos}" "${intergenic_end_pos}" > "${outdir}/pga.anchor.locus.tsv"
awk '{print $2"\t"$3"\t"$4"\t"$1}' "${outdir}/pga.anchor.locus.tsv"   > "${outdir}/pga.anchor.locus.bed"

twoBitToFa -bed="${outdir}/pga.anchor.locus.bed" -noMask   "${genome}" "${outdir}/pga.anchor.locus.fa"

if seqkit seq -w 0 "${outdir}/pga.anchor.locus.fa" |
  grep -v "^>" | grep -qiE 'N{50,}'; then
  echo "GAP warning: contains a continuous gap of N >= 50 bp" >&2
  exit 1
fi

n_ratio=$(grep -v "^>" "${outdir}/pga.anchor.locus.fa" | awk '
{total += length($0); n += gsub(/[Nn]/,"")}
END{if (total == 0) print 100; else print (n/total)*100}
')
if (( $(echo "${n_ratio} > 5.0" | bc -l) )); then
  echo "GAP warning: N ratio (${n_ratio}%) exceeds 5%" >&2
  exit 1
fi

if [[ -n "${gepard_jar}" || -n "${gepard_matrix}" ]]; then
  [[ -s "${gepard_jar}" && -s "${gepard_matrix}" ]] || {
    echo "ERROR: provide both Gepard JAR and matrix, or neither" >&2
    exit 1
  }
  java -Djava.awt.headless=true -cp "${gepard_jar}" \
    org.gepard.client.cmdline.CommandLine \
    -seq "${outdir}/pga.anchor.locus.fa" "${outdir}/pga.anchor.locus.fa" \
    -matrix "${gepard_matrix}" \
    -outfile "${outdir}/${species}.synteny.bmp" -format bmp -word 25
else
  echo "INFO: Gepard paths not supplied; skipping self-alignment QC image"
fi

if [[ -n "${legacy_cluster_script}" ]]; then
  [[ -s "${legacy_cluster_script}" ]] || {
    echo "ERROR: missing legacy clustering script: ${legacy_cluster_script}" >&2
    exit 1
  }
  [[ -s "${protein}" ]] || { echo "ERROR: missing protein FASTA: ${protein}" >&2; exit 1; }
  python3 "${legacy_cluster_script}" \
    "${outdir}/all.VPS37C-VWCE.locus_annotations.full_exon.bed" \
    "${outdir}/all.VPS37C-VWCE.transcript_clusters.tsv"
  sed '1d' "${outdir}/all.VPS37C-VWCE.transcript_clusters.tsv" | cut -f5 \
    > "${outdir}/all.VPS37C-VWCE.transcript_clusters.id.txt"
  seqkit grep -w 0 -f "${outdir}/all.VPS37C-VWCE.transcript_clusters.id.txt" "${protein}" |
    seqkit grep -s -r -p "[^*]+\*$" |
    seqkit replace -p "(.*)" -r "\$1#${species}" \
      -o "${outdir}/all.VPS37C-VWCE.transcript_clusters.final_fixed.fasta"
else
  echo "INFO: legacy clustering script not supplied; proceed with 02_run_refine_toga_local_units.sh"
fi

echo "Anchor locus extraction completed: ${outdir}"
