install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/MC_vcf_heatmaps_0402225/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/MC_vcf_heatmaps_0402225/combine_HOR_flank_trees/chr12.HOR.5kb.pq.hprc-v2.0-mc-chm13.wave.hap_separated_HOR_flank_dist_weighted.nwk")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/MC_vcf_heatmaps_0402225/tanglegrams/chr12.HOR.5kb.pq_vs_chr12_cenhap_refined_tree.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()