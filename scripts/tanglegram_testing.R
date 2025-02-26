install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)
# create trees
# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/chr12_HOR_flank_dist_weighted.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_50kb_flanks_all_pairs_combined_dist_HOR_flank_dist_weighted.nwk")

samples <- readLines("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/centrolign_all_pairs_flanks/HOR_vs_50kb_tree.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()