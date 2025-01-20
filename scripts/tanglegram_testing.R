install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/debug_branch_lengths/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.subgroups_0_and_4.nwk.txt")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/debug_branch_lengths/replace_subtree_subgroup_0_4.nwk.txt")

samples <- readLines("/Users/miramastoras/Desktop/debug_branch_lengths/subgroups_0_and_4.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/debug_branch_lengths/cophylo_0_4.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()