## Tanglegram between HOR NJ tree and refined tree

```R
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/chr12_r2_centrolign_all_pairs_nj_tree.format5.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk")

samples <- readLines("/Users/miramastoras/Desktop/HPRC_release2_contiguous_HORs_chr12.excl_HG00741.1.txt")

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
# all pairs tree
python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/centrolign_refined_tree/pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/centrolign_refined_tree/chr12_79_refined_tree_pairwise_consistency.txt

# refined tree
python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/centrolign_HOR_all_pairs_NJ/pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/top_subtree/centrolign_HOR_all_pairs_NJ/chr12_79_HOR_all_pairs_NJ_pairwise_consistency.txt
```

## How different are the centrolign alignments when using just the HOR tree vs the refined tree?

Prepare a test set of 50 randomly selected samples from chr12 r2

```sh
shuf -n 50 /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HPRC_release2_contiguous_HORs_chr12.fasta_list.txt > /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/HPRC_release2_contiguous_HORs_chr12.fasta_list.shuf50.txt

cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2

cat /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/HPRC_release2_contiguous_HORs_chr12.fasta_list.shuf50.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/HPRC_release2_HORs_chr12.shuf50.fasta
samtools faidx /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/HPRC_release2_HORs_chr12.shuf50.fasta
```
Run centrolign using HOR tree, and using refined tree
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_r2_chr12_allpairs_NJ_HOR_tree_shuf50
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr12_shuf50_HOR_all_pairs_NJ/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr12_r2_centrolign_all_pairs_nj_tree.format5.nwk \
    -A /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr12_shuf50_HOR_all_pairs_NJ/pairwise_cigars/pairwise_cigar \
    -R \
    --threads 32 \
    /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/HPRC_release2_HORs_chr12.shuf50.fasta \
    > /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr12_shuf50_HOR_all_pairs_NJ/HPRC_r2.chr12.allpairs_HOR_NJ.shuf50.centrolign.gfa
```
Refined tree
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_r2_chr12_refined_tree_shuf50
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr12_shuf50_refined_tree/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk \
    -A /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr12_shuf50_refined_tree/pairwise_cigars/pairwise_cigar \
    -R \
    --threads 32 \
    /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/HPRC_release2_HORs_chr12.shuf50.fasta \
    > /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr12_shuf50_refined_tree/HPRC_r2.chr12.allpairs_HOR_NJ.shuf50.centrolign.gfa
```
Run pairwise consistency script
```sh
python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/benchmarking/HOR_vs_refined_tree/chr12_shuf50_HOR_all_pairs_NJ/pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_cigar/pairwise_cigar_ \
    > chr12_shuf50_HOR_all_pairs_NJ_pairwise_consistency.txt
```

```
python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/test/pairwise_consistency/induced/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/test/pairwise_consistency/direct/pairwise_cigar_ \
    > pairwise_consistency.txt
```
