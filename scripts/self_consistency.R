library(ggplot2)
library(scales)  # for alpha()

args <- commandArgs(trailingOnly = TRUE)

# Command-line arguments:
# 1 = pairwise consistency table
# 2 = patristic distance file
# 3 = chromosome/label
# 4 = alignment distance file
# 5 = sample combinations file
# 6 = output prefix

pairwise <- args[1]
patristic <- args[2]
chr <- args[3]
alignment_distance_file <- args[4]
sample_combos_file <- args[5]
outPNG <- args[6]


# ---- Read input data ----
# Pairwise consistency table
dat <- read.table(pairwise, header = TRUE)

# Pairwise patristic distances
dists <- read.csv(patristic, header = TRUE, sep = "\t")
key1 <- paste(dists$sample1, dists$sample2, sep = "_")
key2 <- paste(dists$sample2, dists$sample1, sep = "_")
dists <- rbind(dists, dists)
row.names(dists) <- c(key1, key2)

sample_dists <- dists[paste(dat$sample1, dat$sample2, sep = "_"), "distance"]
dat[["dist"]] <- sample_dists

# Alignment distance file
aln_dists <- read.csv(alignment_distance_file, header = FALSE, sep = ",")
colnames(aln_dists) <- c("sample1", "sample2", "distance")

key1 <- paste(aln_dists$sample1, aln_dists$sample2, sep = "_")
key2 <- paste(aln_dists$sample2, aln_dists$sample1, sep = "_")
aln_dists <- rbind(aln_dists, aln_dists)
row.names(aln_dists) <- c(key1, key2)

sample_aln_dists <- aln_dists[paste(dat$sample1, dat$sample2, sep = "_"), "distance"]
dat[["aln_dist"]] <- sample_aln_dists

# ---- Filter by specified sample combinations ----
combos <- read.csv(sample_combos_file, header = FALSE)
colnames(combos) <- c("sample1", "sample2","dist")

# Create keys for filtering (include both directions)
combo_keys <- c(paste(combos$sample1, combos$sample2, sep = "_"),
                paste(combos$sample2, combos$sample1, sep = "_"))

dat_keys <- paste(dat$sample1, dat$sample2, sep = "_")

# Filter dat to only include specified combinations
dat <- dat[dat_keys %in% combo_keys, ]

cat("Number of pairs after filtering:", nrow(dat), "\n")

# ---- Plotting ----
png(paste(outPNG, chr, "jaccard_hist.png", sep = "_"), width = 800, height = 600)
plot(hist(dat$jaccard, breaks = 100), xlim = c(0, 1.1), xlab = "Jaccard Similarity",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.8, main = paste(chr, "All Cigar Positions (Filtered)", sep=" "))
dev.off()

png(paste(outPNG, chr, "jaccard_vs_dist.png", sep = "_"), width = 800, height = 600)
plot(dat$dist, dat$jaccard, pch = 19, col = alpha("black", 0.1),
     xlab = "Patristic distance (HOR NJ Tree)", ylab = "Jaccard similarity", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.8,
     main = paste(chr, "All Cigar Positions (Filtered)", sep=" "))
dev.off()

png(paste(outPNG, chr, "aligned_jaccard_hist.png", sep = "_"), width = 800, height = 600)
plot(hist(dat$aligned_jaccard, breaks = 100), xlim = c(0, 1.1), xlab = "Jaccard similarity",
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.8, main = paste(chr, "Only aligned pairs (Filtered)", sep=" "))
dev.off()

png(paste(outPNG, chr, "aligned_jaccard_vs_dist.png", sep = "_"), width = 800, height = 600)
plot(dat$dist, dat$aligned_jaccard, pch = 19, col = alpha("black", 0.1),
     xlab = "Patristic distance (HOR NJ Tree)", ylab = "Jaccard similarity",
     main = paste(chr, "Only aligned pairs (Filtered)", sep=" "), cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.8)
dev.off()

png(paste(outPNG, chr, "aligned_jaccard_vs_centrolign_dist.png", sep = "_"), width = 800, height = 600)
plot(dat$aln_dist, dat$aligned_jaccard, pch = 19, col = alpha("black", 0.1),
     xlab = "Centrolign direct pairwise distance", ylab = "Jaccard similarity",
     main = paste(chr, "Only aligned pairs (Filtered)", sep=" "), cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.8)
dev.off()

png(paste(outPNG, chr, "jaccard_vs_centrolign_dist.png", sep = "_"), width = 800, height = 600)
plot(dat$aln_dist, dat$jaccard, pch = 19, col = alpha("black", 0.1), 
     xlab = "Centrolign direct pairwise distance", ylab = "Jaccard similarity",
     main = paste(chr, "All Cigar Positions (Filtered)", sep=" "), cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.8)
dev.off()
