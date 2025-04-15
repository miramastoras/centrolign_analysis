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