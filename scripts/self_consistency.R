library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
pairwise <- args[1]
patristic <- args[2]
chr <- args[3]
outPNG <- args[4]


# pairwise consistency table
dat = read.table(pairwise, header = T)

# pairwise patristic distances from tree (using centrolign tree pair dist)
dists = read.csv(patristic, header = T, sep = "\t")

key1 = paste(dists$sample1, dists$sample2, sep = "_")
key2 = paste(dists$sample2, dists$sample1, sep = "_")

dists = rbind(dists, dists)
row.names(dists) = c(key1, key2)

sample_dists = dists[paste(dat$sample1, dat$sample2, sep = "_"), "distance"]

dat[["dist"]] = sample_dists

#plot(density(dat$jaccard))

png(paste(outPNG, chr, "jaccard_hist.png", sep = "_"), width = 800, height = 600)
plot(hist(dat$jaccard, breaks = 100), xlim = c(0, 1.1),xlab = "Jaccard")
dev.off()

png(paste(outPNG, chr, "jaccard_vs_dist.png", sep = "_"), width = 800, height = 600)
plot(dat$dist, dat$jaccard, pch = 19, col = alpha("black", 0.1), xlim = c(0, 1.1), ylim = c(0, 1.1),
     xlab = "Patristic distance", ylab = "Jaccard similarity")
dev.off()

png(paste(outPNG, chr, "aligned_jaccard_hist.png", sep = "_"), width = 800, height = 600)
plot(hist(dat$aligned_jaccard, breaks = 100),xlim = c(0, 1.1),xlab = "Aligned Jaccard",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.8)
dev.off()

png(paste(outPNG, chr, "aligned_jaccard_vs_dist.png", sep = "_"), width = 800, height = 600)
plot(dat$dist, dat$aligned_jaccard, pch = 19, col = alpha("black", 0.1),  xlim = c(0, 1.1), ylim = c(0, 1.1),
     xlab = "Patristic distance", ylab = "Jaccard similarity", main = "Only aligned pairs",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.8)
dev.off()
