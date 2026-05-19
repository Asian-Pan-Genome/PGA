# Selecting easy regions
```shell
bedtools merge -d 5000 -i GRCh38_Easyregions.bed | awk '$3-$2>=1000000' | awk 'OFS="\t" {print $1,$2,$2+1000000}' | rg '^chr11' | sort -R | head -20 > hg38.regions.bed
bedtools getfasta -fi GRCh38.p14.genome.fa -bed hg38.regions.bed | seqkit replace -p '(.+)' -r 'GRCh38.region{nr}' -w 0 > hg38.regions.fa
```

# Extracting aligned sequence
```shell
minimap2 -x asm5 $sample.$hap.chr11.fa $home/hg38.regions.fa > $sample.$hap.paf -c --secondary=no -t8
paftools.js liftover -q 0 $sample.$hap.paf <(echo -e "GRCh38.region$i\t0\t1000000") | bedtools sort -i - > tmp.bed

chrom=`head -1 tmp.bed | cut -f 1`
start=`head -1 tmp.bed | cut -f 2`
end=`tail -1 tmp.bed | cut -f 3`
fxTools getseq -r <(echo -e "$chrom\t$start\t$end") $sample.$hap.chr11.fa | seqkit replace -p '(.+)' -r "${b[0]}.${b[1]}.region$i" > $sample.$hap.region$i.fa
```

# Constructing Minigraph-Cactus graph and decomposing variants

# LD pruning and PCA extraction
```shell
plink --vcf all_regions.SNP.hap.vcf.gz \
      --make-bed --out all_regions.SNP.hap \
      --const-fid --allow-extra-chr \
      --vcf-half-call haploid \
      --threads 32
plink --bfile all_regions.SNP.hap \
      --indep-pairwise 50 5 0.2 \
      --out all_regions.SNP.hap.pruning \
      --allow-extra-chr \
      --threads 32
plink --bfile all_regions.SNP.hap \
      --extract all_regions.SNP.hap.pruning.prune.in \
      --pca 10 \
      --out all_regions.SNP.hap.PCA \
      --allow-extra-chr \
      --threads 32
```
