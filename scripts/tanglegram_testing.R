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