# Population structure

Joint PCA of diploid assembly genotypes and the 1KGP high-coverage panel using biallelic SNPs from 20 1-Mbp GIAB easy regions on chromosome 11.

## Assembly-derived variants

Project the 20 GRCh38 regions to each assembly and extract the corresponding sequences:

```bash
bash extract_assembly_regions.sh \
    assemblies.tsv \
    GRCh38.regions.fa \
    assembly_regions \
    8
```

For each region, `run_minigraph_cactus.sh` constructs a Minigraph–Cactus graph and processes graph-derived variants using the prepare-vcf-MC/collapse-bubble workflow and `vcfwave`:

```bash
bash run_minigraph_cactus.sh \
    assemblies.tsv \
    GRCh38.regions.fa \
    assembly_regions \
    graph \
    32
```

## Joint assembly–1KGP PCA

```bash
bash run_joint_pca.sh \
    assemblies.tsv \
    hg38.regions.bed \
    GRCh38.regions.fa \
    graph \
    1KGP.chr11.vcf.gz \
    PGA.copies.tsv \
    assemblies.new_superpop.list \
    1KGP.metadata.txt \
    joint_pca \
    32
```

Assembly-derived and 1KGP SNPs are matched by region, position, REF and ALT. Sites with strong cross-panel discordance are removed using individuals represented in both datasets. The joint panel is filtered at MAF ≥ 0.05 and genotype missingness ≤ 0.05, LD-pruned with `--indep-pairwise 50 5 0.2`, and analysed with PLINK v1.9.

Plot the first two PCs with:

```bash
python plot_joint_pca.py \
    --eigenvec joint_pca/plink/joint.PCA.eigenvec \
    --metadata joint_pca/joint_pca.metadata.tsv \
    --out joint_pca/joint.PCA
```
