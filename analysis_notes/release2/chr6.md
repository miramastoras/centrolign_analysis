## Centrolign analysis of chr 6 for release 2

### 0. Prepare data

#### Cenhap distances

Cenhap guide trees and distance matrix provided by Chuck and Sasha Langley on 4/11/2025. These 7 samples are removed due to large deletions:
```
NA20827.1
HG00642.2
HG00344.2
HG00639.1
HG02622.1
HG02145.1
NA20799.2
```
distance matrix from p and q
```sh
/private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr6_58200000_58286706_61058390_61121735_het68_m_mira_dgp_rnj.m
```
Reformat matrix file to enable reading into python script
```sh
# convert spaces to comma, remove trailing commas
sed 's/ /,/g' HPRC_chr6_58200000_58286706_61058390_61121735_het68_m_mira_dgp_rnj.m | sed 's/,$//' | sed 's/\t/,/' | tail -n +2 > matrix_no_header.csv

# convert row names to column names
N=184
cut -f1 -d"," matrix_no_header.csv | sed 's/$/,/g' | tr -d '\n' | cut -f1-$N -d"," | sed 's/^/0,/g' > rowname_to_colname.csv

cat rowname_to_colname.csv matrix_no_header.csv > HPRC_chr6_58200000_58286706_61058390_61121735_het68_m_mira_dgp_rnj.reformatted.m

rm matrix_no_header.csv rowname_to_colname.csv
```
#### Centrolign all pairs distances

Convert cigar to distance matrix
```sh
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/cigar_to_distance.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_cigar/

# header required for next script
echo "sample1,sample2,distance" > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_distance.csv

# remove sample not in tree
grep -v "NA20827.1" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_cigar/pairwise_distance.csv | grep -v "HG00642.2" | grep -v "HG00344.2" | grep -v "HG00639.1" | grep -v "HG02622.1" | grep -v "HG02145.1" | grep -v "NA20799.2" >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_distance.csv
```
### 1. Create refined tree combining centrolign all pairs distances and cenhap distances

Run script to combine distances with weighted sum and compute NJ tree
```sh
conda activate skbio

python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/combine_HOR_flank_dist.py \
    -c /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_distance.csv \
    -f /private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr6_58200000_58286706_61058390_61121735_het68_m_mira_dgp_rnj.reformatted.m \
    -o /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr6_cenhap_20250414_centrolign_all_pairs
```
Convert tree to ete3 format 5 which is accepted by centrolign

```py
from ete3 import Tree

tree = Tree("/private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr6_cenhap_20250414_centrolign_all_pairs_HOR_flank_dist_weighted.nwk", format=1)
tree.write(outfile="/private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr6_cenhap_20250414_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk", format=5)
```
Remove 7 samples from fasta file
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2

cat /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HPRC_release2_contiguous_HORs_chr6.fasta_list.txt | grep -v "HG00642.2" | grep -v "HG00344.2" | grep -v "HG00639.1" | grep -v "HG02622.1" | grep -v "HG02145.1" | grep -v "NA20799.2" | grep -v "NA20827.1" | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr6.excl_7_del_samples.fasta
samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr6.excl_7_del_samples.fasta
```
Run centrolign on refined tree
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_release2_chr6_refined_tree
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_refined_tree/chr6/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr6_cenhap_20250414_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr6.excl_7_del_samples.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_refined_tree/chr6/HPRC_r2.chr6.refined_tree.centrolign.gfa
```

Infer a tree from centrolign all pairs alignments of HOR
```sh
conda activate skbio

python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_cigar/ > /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr6_r2_centrolign_all_pairs_nj_tree.nwk
```

Reformat for centrolign
```py
from ete3 import Tree

tree = Tree("/private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr6_r2_centrolign_all_pairs_nj_tree.nwk", format=1)
tree.write(outfile="/private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr6_r2_centrolign_all_pairs_nj_tree.format5.nwk", format=5)
```
Run centrolign
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_release2_chr6_allpairs_NJ_HOR_tree
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_HOR_all_pairs_NJ_tree/chr6/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr6_r2_centrolign_all_pairs_nj_tree.format5.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr6.excl_7_del_samples.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_HOR_all_pairs_NJ_tree/chr6/HPRC_r2.chr6.allpairs_HOR_NJ.centrolign.gfa
```

Plot tanglegram between refined tree and HOR tree:
```R
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/chr6_r2_centrolign_all_pairs_nj_tree.format5.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/HPRC_r2_chr6_cenhap_20250414_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/chr6_r2_all_pairs_vs_refined_tree.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```

##### Plot pairwise distances

Select 130 random samples from the list for the heatmap
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs

cat HPRC_release2_contiguous_HORs_chr6.txt | grep -v "HG00642.2" | grep -v "HG00344.2" | grep -v "HG00639.1" | grep -v "HG02622.1" | grep -v "HG02145.1" | grep -v "NA20799.2" | grep -v "NA20827.1" | shuf -n 130 > HPRC_release2_contiguous_HORs_chr6.shuf_130.txt
```

Plot pairwise heatmap for the HOR tree
```sh
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_r2_centrolign_all_pairs_nj_tree.format5.nwk \
        -s /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/HPRC_release2_contiguous_HORs_chr6.shuf_130.txt \
        -p /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_pairwise_distance.csv \
        -m "Centrolign all pairs distances" \
        -n "Centrolign all pairs HOR NJ Tree" \
        -d "All pairs Distances" \
        -o /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_centrolign_all_pairs.shuf130.HOR_NJ_tree_

```
Plot pairwise heatmap for the refined tree
```sh
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/HPRC_r2_chr6_cenhap_20250414_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk \
        -s /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/HPRC_release2_contiguous_HORs_chr6.shuf_130.txt \
        -p /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_pairwise_distance.csv \
        -m "Centrolign all pairs distances" \
        -n "Centrolign all pairs refined tree" \
        -d "weighted distances" \
        -o /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_centrolign_all_pairs.shuf130.refined_tree_
```
### 2.2 Run linear regression on all the trees

```sh
cat /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/HPRC_release2_contiguous_HORs_chr6.txt | grep -v "HG00642.2" | grep -v "HG00344.2" | grep -v "HG00639.1" | grep -v "HG02622.1" | grep -v "HG02145.1" | grep -v "NA20799.2" | grep -v "NA20827.1" > tmp ; mv tmp /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/HPRC_release2_contiguous_HORs_chr6.txt
```

HOR tree
```sh
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/regression_tree_eval.py \
    -t /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_r2_centrolign_all_pairs_nj_tree.format5.nwk \
    -s /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/HPRC_release2_contiguous_HORs_chr6.txt \
    -d /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_pairwise_distance.csv \
    -f col \
    -o /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_all_pairs_NJ_tree_vs_centrolign_all_pairs_

# results
Number of internal nodes: 183
Number of leaves: 184
Root Mean Squared Error: 0.03914495287746881

```
Plot linear regression
```sh
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_r2_centrolign_all_pairs_nj_tree.format5.nwk \
        -s /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/HPRC_release2_contiguous_HORs_chr6.shuf_130.txt \
        -p /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_all_pairs_NJ_tree_vs_centrolign_all_pairs_residuals_pairwise.csv \
        -m "Centrolign all pairs residuals" \
        -n "Centrolign all pairs HOR NJ Tree" \
        -d "NJ Distances" \
        -o /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_centrolign_all_pairs_residuals.HOR_nj_tree_
```

HOR tree
```sh
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/regression_tree_eval.py \
    -t /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/HPRC_r2_chr6_cenhap_20250414_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk \
    -s /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/HPRC_release2_contiguous_HORs_chr6.txt \
    -d /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_pairwise_distance.csv \
    -f col \
    -o /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_all_refined_tree_vs_centrolign_all_pairs_

# results
Number of internal nodes: 183
Number of leaves: 184
Root Mean Squared Error: 0.03493871852999413
```
Plot linear regression
```sh
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/HPRC_r2_chr6_cenhap_20250414_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk \
        -s /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/HPRC_release2_contiguous_HORs_chr6.shuf_130.txt \
        -p /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_all_refined_tree_vs_centrolign_all_pairs_residuals_pairwise.csv \
        -m "Centrolign all pairs residuals" \
        -n "Centrolign refined Tree" \
        -d "weighted distances" \
        -o /Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_centrolign_all_pairs_residuals.refined_tree_
```
Tanglegram for the 130 samples
```R
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_centrolign_all_pairs.shuf130.HOR_NJ_tree_pruned_tree.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_centrolign_all_pairs.shuf130.refined_tree_pruned_tree.nwk")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_r2_all_pairs_vs_refined_tree.shuf130.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```

#### Pairwise consistency
Run pairwise consistency script between induced pairwise cigars and direct pairwise cigars
```sh
python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_HOR_all_pairs_NJ_tree/chr6/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/benchmarking/release2/pairwise_consistency/chr6_HOR_all_pairs_tree_pairwise_consistency.txt
```
Tree distance:
```
/private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr6_r2_centrolign_all_pairs_nj_tree.nwk.pair_dists.tsv
```
```R
library(ggplot2)

# read in data for refined tree
dat = read.table("/Users/miramastoras/Desktop/chr6_tree_consistency_69/chr6_69_refined_tree_pairwise_consistency.txt", header = T)
dists = read.tsv("/private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr6_r2_centrolign_all_pairs_nj_tree.nwk.pair_dists.tsv", header = T)

key1 = paste(dists$sample1, dists$sample2, sep = "_")
key2 = paste(dists$sample2, dists$sample1, sep = "_")

dists = rbind(dists, dists)
row.names(dists) = c(key1, key2)

sample_dists = dists[paste(dat$sample1, dat$sample2, sep = "_"), "distance"]

dat[["dist"]] = sample_dists

#plot(density(dat$jaccard))

plot(hist(dat$jaccard, breaks = 100), xlim = c(0, 1.1), ylim = c(0, 250))

plot(dat$dist, dat$jaccard, pch = 19, col = alpha("black", 0.1), xlim = c(0, 1.1), ylim = c(0, 1.1),
     xlab = "Patristic distance", ylab = "Jaccard similarity", main = "All pairs (chr12 refined tree)")

plot(hist(dat$aligned_jaccard, breaks = 100),xlim = c(0, 1.1), ylim = c(0, 250))

plot(dat$dist, dat$aligned_jaccard, pch = 19, col = alpha("black", 0.1),  xlim = c(0, 1.1), ylim = c(0, 1.1),
     xlab = "Patristic distance", ylab = "Jaccard similarity", main = "Only aligned pairs (chr12 refined tree)")

```
