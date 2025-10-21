library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

# Read in file paths from command line
pairwise <- args[1]
patristic <- args[2]
chr <- args[3]
alignment_distance_file <- args[4]
outPNG <- args[5]


# Read pairwise consistency table
dat = read.table(pairwise, header = TRUE)

# Read pairwise patristic distances
dists = read.csv(patristic, header = TRUE, sep = "\t")
key1 = paste(dists$sample1, dists$sample2, sep = "_")
key2 = paste(dists$sample2, dists$sample1, sep = "_")
dists = rbind(dists, dists)
row.names(dists) = c(key1, key2)
print(length(row.names))
sample_dists = dists[paste(dat$sample1, dat$sample2, sep = "_"), "distance"]
dat[["dist"]] = sample_dists

# Read alignment distance file and merge
aln_dists = read.csv(alignment_distance_file, header = FALSE, sep = ",")
print(head(aln_dists))
key1 = paste(aln_dists$sample1, aln_dists$sample2, sep = "_")
key2 = paste(aln_dists$sample2, aln_dists$sample1, sep = "_")
aln_dists = rbind(aln_dists, aln_dists)
row.names(aln_dists) = c(key1, key2)
print(length(row.names))
sample_aln_dists = aln_dists[paste(dat$sample1, dat$sample2, sep = "_"), "distance"]
dat[["aln_dist"]] = sample_aln_dists  # add alignment distance column

# Plotting
png(paste(outPNG, chr, "jaccard_hist.png", sep = "_"), width = 800, height = 600)
plot(hist(dat$jaccard, breaks = 100), xlim = c(0, 1.1), xlab = "Jaccard")
dev.off()

png(paste(outPNG, chr, "jaccard_vs_dist.png", sep = "_"), width = 800, height = 600)
plot(dat$dist, dat$jaccard, pch = 19, col = alpha("black", 0.1), xlim = c(0, 1.1), ylim = c(0, 1.1),
     xlab = "Patristic distance", ylab = "Jaccard similarity")
dev.off()

png(paste(outPNG, chr, "aligned_jaccard_hist.png", sep = "_"), width = 800, height = 600)
plot(hist(dat$aligned_jaccard, breaks = 100), xlim = c(0, 1.1), xlab = "Aligned Jaccard",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.8)
dev.off()

png(paste(outPNG, chr, "aligned_jaccard_vs_dist.png", sep = "_"), width = 800, height = 600)
plot(dat$dist, dat$aligned_jaccard, pch = 19, col = alpha("black", 0.1), xlim = c(0, 1.1), ylim = c(0, 1.1),
     xlab = "Patristic distance", ylab = "Jaccard similarity", main = "Only aligned pairs",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.8)
dev.off()

## plot centrolign pairwise direct distances vs jaccard
png(paste(outPNG, chr, "aligned_jaccard_vs_centrolign_dist.png", sep = "_"), width = 800, height = 600)
plot(dat$aln_dist, dat$aligned_jaccard, pch = 19, col = alpha("black", 0.1),  xlim = c(0, 1.1), ylim = c(0, 1.1),
     xlab = "Centrolign direct pairwise distance", ylab = "Jaccard similarity", main = "Only aligned pairs",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.8)
dev.off()

## plot centrolign pairwise direct distances vs jaccard
png(paste(outPNG, chr, "jaccard_vs_centrolign_dist.png", sep = "_"), width = 800, height = 600)
plot(dat$aln_dist, dat$jaccard, pch = 19, col = alpha("black", 0.1),  xlim = c(0, 1.1), ylim = c(0, 1.1),
     xlab = "Centrolign direct pairwise distance", ylab = "Jaccard similarity",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.8)
dev.off()
