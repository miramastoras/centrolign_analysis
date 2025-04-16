### Creating centrolign flank distances with MC multi sample vcf

CHM13 chr 12 HOR array coordinates (merged, produced during [simulations analysis](https://github.com/miramastoras/centrolign_analysis/blob/main/analysis_notes/simulations/chm13_simulations.md#1-identify-cyclic-hor-arrays-without-svs))
```
/private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/chm13v2.0.labels.as_hor.chr12.active.mrg.bed
```
Create bedfiles of the p and q arm of this array, with varying sizes
```sh
cd /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/chm13/

# 50 kb downstream (q arm)
awk 'BEGIN{OFS="\t"} {print $1, $3, $3+50000}' \
  /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/chm13v2.0.labels.as_hor.chr12.active.mrg.bed \
  > chm13v2.0.as_hor.chr12.active.50kb.q.bed

# 50 kb upsteam (p arm)
awk 'BEGIN{OFS="\t"} {start=$2-50000; if(start<0) start=0; print $1, start, $2}' \
    /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/chm13v2.0.labels.as_hor.chr12.active.mrg.bed \
    > chm13v2.0.as_hor.chr12.active.50kb.p.bed

# 20 kb downstream (q arm)
awk 'BEGIN{OFS="\t"} {print $1, $3, $3+20000}' \
  /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/chm13v2.0.labels.as_hor.chr12.active.mrg.bed \
  > chm13v2.0.as_hor.chr12.active.20kb.q.bed

# 20 kb upsteam (p arm)
awk 'BEGIN{OFS="\t"} {start=$2-20000; if(start<0) start=0; print $1, start, $2}' \
    /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/chm13v2.0.labels.as_hor.chr12.active.mrg.bed \
    > chm13v2.0.as_hor.chr12.active.20kb.p.bed

# 10 kb downstream (q arm)
awk 'BEGIN{OFS="\t"} {print $1, $3, $3+10000}' \
  /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/chm13v2.0.labels.as_hor.chr12.active.mrg.bed \
  > chm13v2.0.as_hor.chr12.active.10kb.q.bed

# 10 kb upsteam (p arm)
awk 'BEGIN{OFS="\t"} {start=$2-10000; if(start<0) start=0; print $1, start, $2}' \
    /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/chm13v2.0.labels.as_hor.chr12.active.mrg.bed \
    > chm13v2.0.as_hor.chr12.active.10kb.p.bed

# 5 kb downstream (q arm)
awk 'BEGIN{OFS="\t"} {print $1, $3, $3+5000}' \
  /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/chm13v2.0.labels.as_hor.chr12.active.mrg.bed \
  > chm13v2.0.as_hor.chr12.active.5kb.q.bed

# 5 kb upsteam (p arm)
awk 'BEGIN{OFS="\t"} {start=$2-5000; if(start<0) start=0; print $1, start, $2}' \
    /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/chm13v2.0.labels.as_hor.chr12.active.mrg.bed \
    > chm13v2.0.as_hor.chr12.active.5kb.p.bed
```
Subset MC vcf to these regions
```sh
cd /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/chm13

for dist in 5kb 10kb 20kb 50kb ;
    do bedtools intersect -header \
    -a /private/groups/cgl/hprc-graphs/hprc-v2.0-feb28/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.wave.vcf.gz \
    -b chm13v2.0.as_hor.chr12.active.${dist}.p.bed \
    > /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_vcfs/chr12.HOR.${dist}.p.hprc-v2.0-mc-chm13.wave.vcf
done

for dist in 5kb 10kb 20kb 50kb ;
    do bedtools intersect -header \
    -a /private/groups/cgl/hprc-graphs/hprc-v2.0-feb28/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.wave.vcf.gz \
    -b chm13v2.0.as_hor.chr12.active.${dist}.q.bed \
    > /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_vcfs/chr12.HOR.${dist}.q.hprc-v2.0-mc-chm13.wave.vcf
done
```
Run bcftools stats on each vcf
```sh
cd /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_vcfs

mkdir -p bcftools_stats

ls chr12* | while read line ;
    do bcftools stats $line > bcftools_stats/$line.txt
  done
```

Use plink to generate distance matrix for each vcf
```sh
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
    dnastack/plink:1.9 plink \
    --vcf /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_vcfs/chr12.HOR.10kb.p.hprc-v2.0-mc-chm13.wave.vcf \
    --distance square \
    --vcf-half-call m \
    --out /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_vcfs/plink_test
```
Plink doesn't like genotypes `0|.` where one hap is missing
📌 Options for --vcf-half-call:
m: Treat half-calls as missing

h: Treat half-calls as homozygous for the specified allele

r: Treat half-calls as reference homozygous

a: Treat half-calls as alternate homozygous


sci kit allele
```
import allel
import numpy as np
from scipy.spatial.distance import pdist, squareform

# Load VCF and extract genotype calls (0, 1, 2 for allele counts)
callset = allel.read_vcf('input.vcf')
gt = allel.GenotypeArray(callset['calldata/GT'])

# Convert to number of alternate alleles per sample per variant
geno_alt = gt.to_n_alt()

# Optional: Filter missing data (mask or impute as needed)
# Here, we'll drop variants with any missing data
mask = np.all(geno_alt >= 0, axis=0)
geno_alt_clean = geno_alt[:, mask]

# Transpose: samples as rows, variants as columns
geno_matrix = geno_alt_clean.T

# Compute pairwise Hamming distances
dist_matrix = squareform(pdist(geno_matrix, metric='hamming'))

# Convert to number of differing positions (instead of fraction)
num_snps = geno_matrix.shape[1]
dist_matrix_counts = dist_matrix * num_snps

# Print or save
print(dist_matrix_counts)
```
Check missingness of vcfs
```sh
ls *.vcf | while read line ; do  
echo $line
grep -v '^#' $line | \
awk -F'\t' '
{
  for (i=10; i<=NF; i++) {
    split($i, gt, ":");
    if (gt[1] == "./.") {
      full_missing++;
    } else if (gt[1] ~ /^[\.0-9]\/\.$/ || gt[1] ~ /^\.[\/|][\.0-9]$/) {
      half_missing++;
    }
    total++;
  }
}
END {
  print "Total genotype calls:", total;
  print "Fully missing (./.):", full_missing;
  print "Half-missing (e.g. 0/., ./1):", half_missing;
}'
done
```
