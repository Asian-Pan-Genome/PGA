#!/usr/bin/env bash
set -euo pipefail

# Build the 27-locus primate PGGB graph and project graph paths onto the
# single-copy Ateles hybridus and Rhinopithecus bieti references.
#
# Usage:
#   run_primate_pggb_untangle.sh [apes_owms.fa] [output_dir]
#
# Environment:
#   THREADS   number of threads (default: 16)
#
# Required executables in PATH:
#   samtools, pggb, odgi, bedtools

fasta=${1:-apes_owms.fa}
outdir=${2:-apes_owms}
threads=${THREADS:-16}

if [[ ! -s "$fasta" ]]; then
    echo "ERROR: input FASTA not found or empty: $fasta" >&2
    exit 1
fi

n_paths=$(grep -c '^>' "$fasta")
if [[ "$n_paths" -ne 27 ]]; then
    echo "ERROR: expected 27 FASTA records, found $n_paths in $fasta" >&2
    exit 1
fi

mkdir -p "$outdir"
samtools faidx "$fasta"

pggb \
    -i "$fasta" \
    -o "$outdir" \
    -n "$n_paths" \
    -p 80 \
    -c 2 \
    -t "$threads" \
    --skip-viz

mapfile -t graphs < <(find "$outdir" -maxdepth 1 -type f -name '*.smooth.final.og' -print)
if [[ ${#graphs[@]} -ne 1 ]]; then
    echo "ERROR: expected exactly one *.smooth.final.og in $outdir, found ${#graphs[@]}" >&2
    printf '  %s\n' "${graphs[@]:-}" >&2
    exit 1
fi

graph=${graphs[0]}
odgi stepindex -i "$graph" -t "$threads"
stepindex="${graph}.stpidx"
if [[ ! -s "$stepindex" ]]; then
    echo "ERROR: odgi step index was not created: $stepindex" >&2
    exit 1
fi

run_untangle() {
    local reference=$1
    local merge_distance=$2
    local output="$outdir/${reference}.${merge_distance}.untangle.txt"

    {
        printf 'query.name\tquery.start\tquery.end\tref.name\tref.start\tref.end\tscore\tinv\tself.cov\tn.th\n'
        odgi untangle \
            -i "$graph" \
            -r "$reference" \
            --threads "$threads" \
            -P \
            -a "$stepindex" \
            -m "$merge_distance" \
        | bedtools sort -i -
    } \
    | awk 'BEGIN { OFS="\t" }
           NR == 1 { print; next }
           $8 == "-" { tmp=$6; $6=$5; $5=tmp }
           { print }' \
    > "$output"
}

for reference in Ateles_hybridus Rhinopithecus_bieti; do
    for merge_distance in 128 256 500 1000 2000; do
        run_untangle "$reference" "$merge_distance"
    done
done

echo "PGGB graph: $graph"
echo "Untangle outputs: $outdir/{Ateles_hybridus,Rhinopithecus_bieti}.{128,256,500,1000,2000}.untangle.txt"
