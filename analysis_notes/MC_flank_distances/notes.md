## Creating centrolign flank distances with MC multi sample vcf

### Testing with chr12 to begin

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

# first remove grch38 from vcf
bcftools view -s ^GRCh38 -Oz /private/groups/cgl/hprc-graphs/hprc-v2.0-feb28/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.wave.vcf.gz > /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/hprc-v2.0-mc-chm13.wave.no38.vcf.gz

for dist in 5kb 10kb 20kb 50kb ;
    do docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups \
    pegi3s/bedtools bedtools intersect -header \
    -a /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/hprc-v2.0-mc-chm13.wave.no38.vcf.gz \
    -b /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/chm13/chm13v2.0.as_hor.chr12.active.${dist}.p.bed \
    > /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_vcfs/chr12.HOR.${dist}.p.hprc-v2.0-mc-chm13.wave.vcf
done

for dist in 5kb 10kb 20kb 50kb ;
    do docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups \
    pegi3s/bedtools bedtools intersect -header \
    -a /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/hprc-v2.0-mc-chm13.wave.no38.vcf.gz \
    -b /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/chm13/chm13v2.0.as_hor.chr12.active.${dist}.q.bed \
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

Merge p and q arm vcfs
```sh
cd /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_vcfs

for dist in 5kb 10kb 20kb 50kb ;
  do bcftools concat chr12.HOR.${dist}.q.hprc-v2.0-mc-chm13.wave.vcf chr12.HOR.${dist}.p.hprc-v2.0-mc-chm13.wave.vcf > chr12.HOR.${dist}.pq.hprc-v2.0-mc-chm13.wave.vcf
done
```
Separate the haplotypes and convert them to homozygous
```sh
cd /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances

ls MC_vcfs/*.vcf | while read line ; do
  name=$(basename "$line")
  name="${name%.*}"

  bcftools view $line | \
  awk '
  BEGIN {
      OFS = "\t";
    }
    {
      if ($0 ~ /^##/) {
          print;  # meta lines unchanged
        } else if ($0 ~ /^#CHROM/) {
          # Header line: print fixed columns
          for (i = 1; i <= 9; i++) {
              printf("%s%s", $i, (i < 9 ? OFS : OFS));
            }
            # Print sample names with .1 and .2
            for (i = 10; i <= NF; i++) {
              printf("%s.1%s%s.2", $i, OFS, $i);
              if (i < NF) {
                  printf("%s", OFS);
                } else {
                  printf("\n");
                }
              }
            } else {
              # Variant data lines
              for (i = 1; i <= 9; i++) {
                printf("%s%s", $i, (i < 9 ? OFS : OFS));
              }
              # Process genotypes
              for (i = 10; i <= NF; i++) {
                split($i, fields, ":");
                split(fields[1], alleles, /[\/|]/);
                # First haplotype (.1)
                fields[1] = alleles[1] "|" alleles[1];
                printf("%s", fields[1]);
                for (j = 2; j <= length(fields); j++) {
                  printf(":%s", fields[j]);
                }
                printf("%s", OFS);
                # Second haplotype (.2)
                fields[1] = alleles[2] "|" alleles[2];
                printf("%s", fields[1]);
                for (j = 2; j <= length(fields); j++) {
                  printf(":%s", fields[j]);
                }
                if (i < NF) {
                  printf("%s", OFS);
              } else {
                printf("\n");
              }
            }
          }
        }' > MC_hap_separated_vcfs/${name}.hap_separated.vcf
done
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
done > missingness.txt
```

Use plink to generate distance matrix for each vcf
```sh
cd /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_hap_separated_vcfs

mkdir -p /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_hap_separated_vcfs/plink_dists

ls *.vcf | while read line ; do
  VCF=`realpath $line`
  NAME=`basename $line .vcf`
  docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
      dnastack/plink:1.9 plink \
      --vcf ${VCF} \
      --distance square \
      --vcf-half-call m \
      --out /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_hap_separated_vcfs/plink_dists/${NAME}
done
```
Reformat plink distance matrix so that first row and first col are sample IDs
```sh
cd /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_hap_separated_vcfs/plink_dists/

mkdir -p work

ls *.dist | while read line ; do

  # convert dist matrix to csv
  sed 's/\t/,/g' ${line} > work/${line}.csv

  # convert row names to column names, set as first line in reformatted csv
  cut -f1 ${line}.id | paste -sd',' |  sed 's/^/0,/g' > ${line}.formatted.csv

  # paste sample names as col 1, and then dist values
  cut -f 1 ${line}.id  | paste -d "," - work/${line}.csv >> ${line}.formatted.csv

done
```
Combine with centrolign all pairs distances and create a new tree

> Ran into an error because all of the q matrix files have some NAN files due to the small number of SNPs.
```sh
cd /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_hap_separated_vcfs/plink_dists
ls *.dist | while read line ; do echo $line ; grep "nan" $line | wc -l ; done
ls *.dist | while read line ; do echo $line ; awk '{ for (i=1; i<=NF; i++) if ($i == "nan") col[i]++ } END { for (i in col) print i, col[i] }' $line | wc -l ; done
```

Prepare input files
```sh
cd /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_hap_separated_vcfs

# Get intersection of MC vcf and chr12 all pairs
bcftools query -l chr12.HOR.10kb.p.hprc-v2.0-mc-chm13.wave.hap_separated.vcf > MC_samples.txt
comm -23 <(sort /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_chr12.txt) <(sort MC_samples.txt )

# HG00272.1 is the only sample in the centrolign list that is missing from MC vcf
echo "sample1,sample2,dist" > chr12_direct_pairwise_distance_excl_HG00272.1.csv
grep -v HG00272.1 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_distance.csv >> chr12_direct_pairwise_distance_excl_HG00272.1.csv
```
Run combined distance script
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_hap_separated_vcfs/combine_HOR_flank_trees

cd /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_hap_separated_vcfs/plink_dists

ls *.p*.formatted.csv | while read line ; do
    MATRIX=`realpath $line`
    SAMPLE=`basename $line .dist.formatted.csv`
    docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
        miramastoras/centromere_scripts:v0.1.4 \
        python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/combine_HOR_flank_dist.py \
        -c /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_hap_separated_vcfs/chr12_direct_pairwise_distance_excl_HG00272.1.csv \
        -f ${MATRIX} \
        -o /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_hap_separated_vcfs/combine_HOR_flank_trees/${SAMPLE}
      done
```

Plot a random sampling of 130 samples against the chr12 all pairs distances in the pairwise heatmap
```sh
shuf -n 130 /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_chr12.txt > /private/groups/patenlab/mira/centrolign/guide_tree_testing/MC_flank_distances/MC_hap_separated_vcfs/combine_HOR_flank_trees/chr12_0402225_random130_samples.txt
```
Plot heatmaps - local computer
```sh
# change format of input tree to 1
conda activate tree_python
mkdir -p /Users/miramastoras/Desktop/MC_vcf_heatmaps_0402225/heatmaps

ls /Users/miramastoras/Desktop/MC_vcf_heatmaps_0402225/combine_HOR_flank_trees/*nwk | while read line ; do
  SAMPLE=`basename $line hprc-v2.0-mc-chm13.wave.hap_separated_HOR_flank_dist_weighted.nwk`
  TITLE=`echo Refined Tree $SAMPLE`
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
          -t ${line} \
          -s /Users/miramastoras/Desktop/MC_vcf_heatmaps_0402225/chr12_0402225_random130_samples.txt \
          -p /Users/miramastoras/Desktop/MC_vcf_heatmaps_0402225/chr12_direct_pairwise_distance_excl_HG00272.1.csv \
          -m "Centrolign all pairs" \
          -n ${TITLE} \
          -d "Combined distances" \
          -o /Users/miramastoras/Desktop/MC_vcf_heatmaps_0402225/heatmaps/${SAMPLE}_
  done
```
For p and pq 5kb, plot tanglegram with the chr12 refined tree
```sh
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/MC_vcf_heatmaps_0402225/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/MC_vcf_heatmaps_0402225/combine_HOR_flank_trees/chr12.HOR.5kb.p.hprc-v2.0-mc-chm13.wave.hap_separated_HOR_flank_dist_weighted.nwk")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/MC_vcf_heatmaps_0402225/tanglegrams/chr12.HOR.5kb.p_vs_chr12_cenhap_refined_tree.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
For p and pq 5kb, plot tanglegram with the chr12 refined tree
```sh
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/MC_vcf_heatmaps_0402225/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/MC_vcf_heatmaps_0402225/combine_HOR_flank_trees/chr12.HOR.5kb.pq.hprc-v2.0-mc-chm13.wave.hap_separated_HOR_flank_dist_weighted.nwk")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/MC_vcf_heatmaps_0402225/tanglegrams/chr12.HOR.5kb.pq_vs_chr12_cenhap_refined_tree.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
Tree from distances from just MC vs sasha's cenhap tree

sasha's cenhap tree vs MC distances 
