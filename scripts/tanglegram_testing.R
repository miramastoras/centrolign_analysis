install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/replace_subtrees/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.chr12_hprc_samples.HPRC_naming.nwk.txt")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.upgma.nwk")

samples <- readLines("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.HPRC_naming.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/KGP4_TRIOS_MAC5_chr12_vs_HPRC_chr12_P_Q_mira.2.25.25.upgma.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()