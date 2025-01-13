#install.packages("ape")
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/replace_subtrees/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.chr12_hprc_samples.nwk.txt")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/replace_subtrees/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.all_samples.tree_method1.all_subgroups.nwk.txt")

samples <- readLines("/Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt")

# cophyloplot
# # create associations matrix 
# a <- as.character(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40))
# b <- as.character(c(14,1,4,1,9,12,2,10,6,3,13,5,14,15,18,19,19,7,14,9,10,11,25,22,21,16,23,24,26,17,1,12,12,21,15,16,21,8,20,21)) 
# association <- cbind(samples, samples)
# 
# png(filename="/Users/miramastoras/Desktop/cophylo.png")
# 
# # plot
# cophyloplot(orginal_tree, new_tree, assoc = association,length.line = 1000, space = 300, gap = 1000,show.tip.label=FALSE)
# dev.off()

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/cophylo_2.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
