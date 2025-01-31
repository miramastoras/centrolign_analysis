install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
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