## Testing refining the guide tree by combining SNP and HOR distances with a weighted sum

### Jordan's formula 1:

HOR distance = h
flank distance = f
```
d = (1 - (1 - h)^2 + f^2) / 2
```

Steps:
- compute HOR and flank distance weighted sum
- pass these distances to infer tree
- run pairwise heatmap

### Compute HOR and flank distance weighted sum, reinfer tree

File locations
```
# flank distance
/Users/miramastoras/Desktop/distance_compare/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_hap.hD.txt

# initial HOR distance from centrolign run with initial guide tree
# local file
/Users/miramastoras/Desktop/tree_heatmap_chr12/pairwise_distance.csv

# phoenix file location
/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/pairwise_distance.csv
```
New script to compute weighted sum and infer tree from it using neighbor joining
```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/combine_HOR_flank_dist.py \
    -c /Users/miramastoras/Desktop/tree_heatmap_chr12/pairwise_distance.csv \
    -f /Users/miramastoras/Desktop/distance_compare/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_hap.hD.txt \
    -o /Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12
```
plot pairwise distances with original centrolign distances
```
conda activate tree_python

python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_HOR_flank_dist_weighted.nwk \
        -s /Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt \
        -p /Users/miramastoras/Desktop/tree_heatmap_chr12/pairwise_distance.csv \
        -o /Users/miramastoras/Desktop/combine_HOR_flank_dist/combine_HOR_flank_dist_chr12_

# plot just trios
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_HOR_flank_dist_weighted.nwk \
        -s /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt  \
        -p /Users/miramastoras/Desktop/tree_heatmap_chr12/pairwise_distance.csv \
        -o /Users/miramastoras/Desktop/plot_all_tree_methods_just_trios/combine_HOR_flank_dist_chr12_tree_original_distances
```
Plot pairwise heatmap with new calculated distances
```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_HOR_flank_dist_weighted.nwk \
        -s /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt \
        -p /Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_HOR_flank_dist_weighted.txt \
        -o /Users/miramastoras/Desktop/plot_all_tree_methods_just_trios/combine_HOR_flank_dist_chr12_tree_jordans_distances
```

Plot tree tanglegram against original UPGMA tree
```
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/replace_subtrees/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.chr12_hprc_samples.nwk.txt")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_HOR_flank_dist_weighted.nwk")

samples <- readLines("/Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/combine_HOR_flank_dist/cophylo_tree_method2_against_upgma_tree.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
Plot tree tanglegram against refined tree with method 1
```
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/replace_subtrees/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.all_samples.tree_method1.all_subgroups.nwk.txt")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_HOR_flank_dist_weighted.nwk")

samples <- readLines("/Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/combine_HOR_flank_dist/cophylo_tree_method2_against_tree_method1.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
### rerun centrolign with new tree

Convert tree to format 5, accepted by centrolign
```py
from ete3 import Tree

tree = Tree("/Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_HOR_flank_dist_weighted.nwk", format=1)
tree.write(outfile="/Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_HOR_flank_dist_weighted.formatted.nwk.txt", format=5)
```

```
#!/bin/bash
#SBATCH --job-name=centrolign_chr12_tree_method2
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
    -S /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/rerun_centrolign_new_tree/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/chr12_HOR_flank_dist_weighted.formatted.nwk.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/initial_test_no_gaps_chr12.fasta \
    > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/rerun_centrolign_new_tree/initial_test_no_gaps_chr12.tree_method2.centrolign.gfa
```
restart and output pairwise alignments
```
#!/bin/bash
#SBATCH --job-name=centrolign_chr12_tree_method2
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
    -S /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/rerun_centrolign_new_tree/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/chr12_HOR_flank_dist_weighted.formatted.nwk.txt \
    -A /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/rerun_centrolign_new_tree/pairwise_cigars/pairwise_cigar \
    -R \
    --threads 32 \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/initial_test_no_gaps_chr12.fasta \
    > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/rerun_centrolign_new_tree/initial_test_no_gaps_chr12.tree_method2.centrolign.gfa
```

Plot pairwise heatmap
```
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/cigar_to_distance.py /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/rerun_centrolign_new_tree/pairwise_cigars/

# on personal computer
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_HOR_flank_dist_weighted.nwk \
        -s /Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt \
        -p /Users/miramastoras/Desktop/combine_HOR_flank_dist/pairwise_distance_centrolign_rerun_tree_method2.csv \
        -o /Users/miramastoras/Desktop/combine_HOR_flank_dist/combine_HOR_flank_dist_chr12_rerun_centrolign

# only trios
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_HOR_flank_dist_weighted.nwk \
        -s /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt \
        -p /Users/miramastoras/Desktop/combine_HOR_flank_dist/pairwise_distance_centrolign_rerun_tree_method2.csv \
        -o /Users/miramastoras/Desktop/plot_all_tree_methods_just_trios/combine_HOR_flank_dist_chr12_rerun_centrolign
```
### Test out another round of inferring the tree and running centrolign

Combine snp and centrolign distances from rerun
```
conda activate skbio

python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/combine_HOR_flank_dist.py \
    -c /Users/miramastoras/Desktop/combine_HOR_flank_dist/pairwise_distance_centrolign_rerun_tree_method2.csv \
    -f /Users/miramastoras/Desktop/distance_compare/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_hap.hD.txt \
    -o /Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_round2
```

PLot pairwise heatmap with weighted distances
```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_round2_HOR_flank_dist_weighted.nwk \
        -s /Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt \
        -p /Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_round2_HOR_flank_dist_weighted.txt \
        -o /Users/miramastoras/Desktop/combine_HOR_flank_dist/combine_HOR_flank_dist_chr12_round2_weighted_dists
```
PLot pairwise heatmap with centrolign distances from previous run
```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_round2_HOR_flank_dist_weighted.nwk \
        -s /Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt \
        -p /Users/miramastoras/Desktop/combine_HOR_flank_dist/pairwise_distance_centrolign_rerun_tree_method2.csv \
        -o /Users/miramastoras/Desktop/combine_HOR_flank_dist/combine_HOR_flank_dist_chr12_round2_centrolign_dists
```
compare round 1 vs round 2 trees
```R
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_HOR_flank_dist_weighted.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_round2_HOR_flank_dist_weighted.nwk")

samples <- readLines("/Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/combine_HOR_flank_dist/cophylo_tree_method2_round1_vs_round2.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```

Rerun centrolign
```
#!/bin/bash
#SBATCH --job-name=centrolign_chr12_tree_method2
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
    -S /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/rerun_centrolign_new_tree_round2/jobstore/ \
    -R \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/chr12_round2_HOR_flank_dist_weighted.formatted.nwk.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/initial_test_no_gaps_chr12.fasta \
    > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/rerun_centrolign_new_tree_round2/initial_test_no_gaps_chr12.tree_method2_round2.centrolign.gfa
```

Re-start and output pairwise alignments
```
time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/rerun_centrolign_new_tree_round2/jobstore/ \
    -A /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/rerun_centrolign_new_tree_round2/pairwise_cigars/pairwise_cigar \
    -R \
    --threads 32 \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/chr12_round2_HOR_flank_dist_weighted.formatted.nwk.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/initial_test_no_gaps_chr12.fasta \
    > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/rerun_centrolign_new_tree_round2/initial_test_no_gaps_chr12.tree_method2_round2.centrolign.gfa
```
Plot pairwise heatmap
```
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/cigar_to_distance.py /private/groups/patenlab/mira/centrolign/guide_tree_testing/method2_weighted_sum/rerun_centrolign_new_tree_round2/pairwise_cigars/

# on personal computer
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_round2_HOR_flank_dist_weighted.formatted.nwk.txt \
        -s /Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt \
        -p /Users/miramastoras/Desktop/combine_HOR_flank_dist/pairwise_distance_centrolign_rerun_tree_method2_round2.csv \
        -o /Users/miramastoras/Desktop/combine_HOR_flank_dist/combine_HOR_flank_dist_chr12_rerun_centrolign_round2_
```
```
# trios only
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/combine_HOR_flank_dist/chr12_round2_HOR_flank_dist_weighted.formatted.nwk.txt \
        -s /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt \
        -p /Users/miramastoras/Desktop/combine_HOR_flank_dist/pairwise_distance_centrolign_rerun_tree_method2_round2.csv \
        -o /Users/miramastoras/Desktop/plot_all_tree_methods_just_trios/combine_HOR_flank_dist_chr12_rerun_centrolign_round2_
```
