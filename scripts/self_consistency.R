library(ggplot2)

# Read in files
args <- commandArgs(trailingOnly = TRUE)
pairwise <- args[1]
patristic <- args[2]
chr <- args[3]
alignment_distance_file <- args[4]  # centrolign alignment distances file
outPNG <- args[5]

# Read pairwise consistency table
dat = read.table(pairwise, header = TRUE)

# === Patristic distances ===
dists = read.csv(patristic, header = TRUE, sep = "\t")

# Build symmetric keys for lookup
key1 = paste(dists$sample1, dists$sample2, sep = "_")
key2 = paste(dists$sample2, dists$sample1, sep = "_")

# Duplicate data to ensure symmetric lookup
dists_reversed = dists
dists_reversed$sample1 = dists$sample2
dists_reversed$sample2 = dists$sample1

# Combine both directions
dists_full = rbind(dists, dists_reversed)
rownames(dists_full) = paste(dists_full$sample1, dists_full$sample2, sep = "_")

# Match and extract patristic distances for each row in dat
sample_dists = dists_full[paste(dat$sample1, dat$sample2, sep = "_"), "distance"]
dat[["dist"]] = sample_dists

# === Alignment distances ===
aln_dists = read.csv(alignment_distance_file, header = TRUE, sep = "\t")

# Create symmetric keys
aln_key1 = paste(aln_dists$sample1, aln_dists$sample2, sep = "_")
aln_key2 = paste(aln_dists$sample2, aln_dists$sample1, sep = "_")

# Duplicate data for reverse pairs
aln_dists_reversed = aln_dists
aln_dists_reversed$sample1 = aln_dists$sample2
aln_dists_reversed$sample2 = aln_dists$sample1

# Combine both directions
aln_dists_full = rbind(aln_dists, aln_dists_reversed)
rownames(aln_dists_full) = paste(aln_dists_full$sample1, aln_dists_full$sample2, sep = "_")

# Match and extract alignment distances
sample_aln_dists = aln_dists_full[paste(dat$sample1, dat$sample2, sep = "_"), "distance"]
dat[["aln_dist"]] = sample_aln_dists

# === Plotting ===

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
