## Tanglegram between HOR NJ tree and refined tree

```R
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/chr12_r2_centrolign_all_pairs_nj_tree.format5.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk")


obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/chr12_r2_all_pairs_vs_refined_tree.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
https://docs.google.com/presentation/d/1T6ln7REQrexSORkSfWkdfYGD8N7yVyZaUz5BNHESh2I/edit?slide=id.p#slide=id.p

Based off the tanglegram, splitting the tree into two pieces and taking the top half where the majority of switching happens
```sh
python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/split_fasta_by_tree.py \
    -t /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr12_r2_centrolign_all_pairs_nj_tree.format5.nwk \
    -f /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr12.excl_HG00741.1.fasta \
    -n 2 \
    -o /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/split_fasta_by_tree/
```
Plotting tanglegram with these samples to make sure its the correct subset
```R
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/chr12_r2_centrolign_all_pairs_nj_tree.format5.subgroup0.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.subgroup0.nwk")

samples <- readLines("/Users/miramastoras/Desktop/subgroup0.samples.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/chr12_r2_all_pairs_vs_refined_tree.subgroup0.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
Split this subgroup into 3 now
```sh
python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/split_fasta_by_tree.py \
    -t /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/split_fasta_by_tree/subgroup_0_tree.nwk \
    -f /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/split_fasta_by_tree/subgroup_0_seqs.fasta \
    -n 11 \
    -o /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/split_fasta_by_tree/split_subgroup0/

# identify subgroup
wc -l *.fai
 122 subgroup_0_seqs.fasta.fai
  17 subgroup_1_seqs.fasta.fai
  53 subgroup_2_seqs.fasta.fai
```
Plotting subgroup 2
```sh
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/chr12_r2_centrolign_all_pairs_nj_tree.format5.shuf50.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.shuf50.nwk")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/chr12_r2_all_pairs_vs_refined_tree.shuf50.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
```
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/chr12_r2_centrolign_all_pairs_nj_tree.format5.90.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.90.nwk")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/chr12_r2_all_pairs_vs_refined_tree.shuf50.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```

Take subgroup 6 and cut it into 3
```
python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/split_fasta_by_tree.py \
    -t /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/split_fasta_by_tree/split_subgroup0/subgroup_6_tree.nwk \
    -f /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/split_fasta_by_tree/split_subgroup0/subgroup_6_seqs.fasta \
    -n 4 \
    -o /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/split_fasta_by_tree/split_subgroup0/split_subgroup6/
```

Run centrolign using HOR tree, and using refined tree
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_r2_chr12_allpairs_NJ_HOR_tree_79
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=72:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/centrolign_HOR_all_pairs_NJ/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr12_r2_centrolign_all_pairs_nj_tree.format5.nwk \
    -A /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/centrolign_HOR_all_pairs_NJ/pairwise_cigars/pairwise_cigar \
    -R \
    --threads 32 \
    /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/subset_tree_79_samples.fasta \
    > /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/centrolign_HOR_all_pairs_NJ/HPRC_r2.chr12.allpairs_HOR_NJ.79.centrolign.gfa
```
Refined tree
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_r2_chr12_refined_tree_79
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=72:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/centrolign_refined_tree/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk \
    -A /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/centrolign_refined_tree/pairwise_cigars/pairwise_cigar \
    -R \
    --threads 32 \
    /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/subset_tree_79_samples.fasta \
    > /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/centrolign_refined_tree/HPRC_r2.chr12.allpairs_HOR_NJ.79.centrolign.gfa
```

Run pairwise consistency script
```sh
# refined tree
python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/centrolign_refined_tree/pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/centrolign_refined_tree/chr12_79_refined_tree_pairwise_consistency.txt

# all pairs tree
python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/centrolign_HOR_all_pairs_NJ/pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/centrolign_HOR_all_pairs_NJ/chr12_79_HOR_all_pairs_NJ_pairwise_consistency.txt
```

Plot using self_consistency.R, and run a paired t test to compare the jaccard distributions  

HOR NJ tree
```R
library(ggplot2)

# read in data for refined tree
dat = read.table("/Users/miramastoras/Desktop/chr12_tree_consistency_79/chr12_79_refined_tree_pairwise_consistency.txt", header = T)
dists = read.csv("/Users/miramastoras/Desktop/chr12_tree_consistency_79/pairwise_distance_excl_HG00741.1.csv", header = T)

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

# read in data for all pairs tree
dat2 = read.table("/Users/miramastoras/Desktop/chr12_tree_consistency_79/chr12_79_HOR_all_pairs_NJ_pairwise_consistency.txt", header = T)

sample_dists2 = dists[paste(dat2$sample1, dat2$sample2, sep = "_"), "distance"]

dat2[["dist"]] = sample_dists2

plot(hist(dat2$jaccard, breaks = 100), xlim = c(0, 1.1), ylim = c(0, 250))

plot(hist(dat2$aligned_jaccard, breaks = 100),xlim = c(0, 1.1), ylim = c(0, 250))

dat_sorted <- dat[order(dat[[1]], dat[[2]]), ]
dat2_sorted <- dat2[order(dat2[[1]], dat2[[2]]), ]
# Paired t-test
t.test(dat$jaccard, dat2$jaccard, paired = TRUE)

t.test(dat$aligned_jaccard, dat2$aligned_jaccard, paired = TRUE)

# find biggest difference to look at their cigar strings
# Compute the absolute difference
diffs <- abs(dat$aligned_jaccard - dat2$aligned_jaccard)

# Find the index of the max difference
max_diff_index <- which.max(diffs)

# View the index and the actual values (optional)
max_diff_index
dat[max_diff_index, ]
dat2[max_diff_index, ]
```

## Subset samples from chr 6

Split the 130 random samples into 2 subtrees

```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2

cat /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_chr6.shuf_130.txt | while read line ; do
    grep $line HPRC_release2_contiguous_HORs_chr6.fasta_list.txt ; done | while read line ; do cat $line ; done >> /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/HPRC_release2_HORs_chr6_shuf_130.fasta

samtools faidx /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/HPRC_release2_HORs_chr6_shuf_130.fasta
```
```sh
python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/split_fasta_by_tree.py \
    -t /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr6_r2_centrolign_all_pairs_nj_tree.format5.nwk \
    -f /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/HPRC_release2_HORs_chr6_shuf_130.fasta \
    -n 2 \
    -o /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/split_fasta_by_tree/
```

Quickly ran a tanglegram to confirm that subgroup 0 is the subtree with more switching between the two trees. https://docs.google.com/presentation/d/1T6ln7REQrexSORkSfWkdfYGD8N7yVyZaUz5BNHESh2I/edit?slide=id.g34c9dee2485_0_104#slide=id.g34c9dee2485_0_104

Run centrolign on subgroup 0 with both trees
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_release2_chr6_69_allpairs_NJ_HOR_tree
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/subset_69_MSA_HOR_all_pairs_NJ/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr6_r2_centrolign_all_pairs_nj_tree.format5.nwk \
    -A /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/subset_69_MSA_HOR_all_pairs_NJ/pairwise_cigars/pairwise_cigar \
    -R \
    --threads 32 \
    /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/split_fasta_by_tree/subgroup_0_seqs.fasta \
    > /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/subset_69_MSA_HOR_all_pairs_NJ/HPRC_r2.chr6_sub69.allpairs_HOR_NJ.centrolign.gfa
```
Refined tree
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_r2_chr6_69_refined_tree
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/subset_69_MSA_refined_tree/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr6_cenhap_20250414_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk \
    -A /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/subset_69_MSA_refined_tree/pairwise_cigars/pairwise_cigar \
    -R \
    --threads 32 \
    /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/split_fasta_by_tree/subgroup_0_seqs.fasta \
    > /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/subset_69_MSA_refined_tree/HPRC_r2.chr6_sub69.refined_tree.centrolign.gfa
```

Run pairwise consistency script
```sh
# all pairs tree
time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/subset_69_MSA_HOR_all_pairs_NJ/pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/subset_69_MSA_HOR_all_pairs_NJ/chr6_69_all_pairs_tree_pairwise_consistency.txt

#real    597m38.658s
#user    546m35.449s
#sys     35m41.605s

# refined tree
time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/subset_69_MSA_refined_tree/pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr6_shuf_130/subset_69_MSA_refined_tree/chr6_69_refined_tree_pairwise_consistency.txt

#real    597m17.650s
#user    545m14.004s
#sys     35m45.879s
```

Plot using self_consistency.R, and run a paired t test to compare the jaccard distributions  

HOR NJ tree
```R
library(ggplot2)

# read in data for refined tree
dat = read.table("/Users/miramastoras/Desktop/chr6_tree_consistency_69/chr6_69_refined_tree_pairwise_consistency.txt", header = T)
dists = read.csv("/Users/miramastoras/Desktop/chr6_tree_consistency_69/chr6_all_pairs_pairwise_distance.csv", header = T)

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

# read in data for all pairs tree
dat2 = read.table("/Users/miramastoras/Desktop/chr6_tree_consistency_69/chr6_69_all_pairs_tree_pairwise_consistency.txt", header = T)

sample_dists2 = dists[paste(dat2$sample1, dat2$sample2, sep = "_"), "distance"]

dat2[["dist"]] = sample_dists2

plot(hist(dat2$jaccard, breaks = 100), xlim = c(0, 1.1), ylim = c(0, 250))

plot(hist(dat2$aligned_jaccard, breaks = 100),xlim = c(0, 1.1), ylim = c(0, 250))

dat_sorted <- dat[order(dat[[1]], dat[[2]]), ]
dat2_sorted <- dat2[order(dat2[[1]], dat2[[2]]), ]
# Paired t-test
t.test(dat$jaccard, dat2$jaccard, paired = TRUE)

t.test(dat$aligned_jaccard, dat2$aligned_jaccard, paired = TRUE)

# find biggest difference to look at their cigar strings
# Compute the absolute difference
diffs <- abs(dat$aligned_jaccard - dat2$aligned_jaccard)

# Find the index of the max difference
max_diff_index <- which.max(diffs)

# View the index and the actual values (optional)
max_diff_index
dat[max_diff_index, ]
dat2[max_diff_index, ]
```
