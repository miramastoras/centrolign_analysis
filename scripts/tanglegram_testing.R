install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/chr6_pairwise_heatmaps/chr6_centrolign_all_pairs.shuf130.HOR_NJ_tree_pruned_tree.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/subgroup_0_tree.nwk")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/chr6_pairwise_heatmaps/subgroup0.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()