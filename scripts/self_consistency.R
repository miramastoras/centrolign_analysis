library(ggplot2)

# read in data for refined tree
dat = read.table("/Users/miramastoras/Desktop/chr6_tree_consistency_69/chr6_69_refined_tree_pairwise_consistency.txt", header = T)
dists = read.csv("/Users/miramastoras/Desktop/chr6_tree_consistency_69/chr6_all_pairs_pairwise_distance.csv", header = T)

key1 = paste(dists$sample1, dists$sample2, sep = "_")
key2 = paste(dists$sample2, dists$sample1, sep = "_")

dists = rbind(dists, dists)
row.names(dists) = c(key1, key2)

sample_dists = dists[paste(dat$sample1, dat$sample2, sep = "_"), "distance"]

dat[["dist"]] = sample_dists

#plot(density(dat$jaccard))

plot(hist(dat$jaccard, breaks = 100), xlim = c(0, 1.1), ylim = c(0, 250))

plot(dat$dist, dat$jaccard, pch = 19, col = alpha("black", 0.1), xlim = c(0, 1.1), ylim = c(0, 1.1),
     xlab = "Patristic distance", ylab = "Jaccard similarity", main = "All pairs (chr12 refined tree)")

plot(hist(dat$aligned_jaccard, breaks = 100),xlim = c(0, 1.1), ylim = c(0, 250))

plot(dat$dist, dat$aligned_jaccard, pch = 19, col = alpha("black", 0.1),  xlim = c(0, 1.1), ylim = c(0, 1.1),
     xlab = "Patristic distance", ylab = "Jaccard similarity", main = "Only aligned pairs (chr12 refined tree)")

# read in data for all pairs tree
dat2 = read.table("/Users/miramastoras/Desktop/chr6_tree_consistency_69/chr6_69_all_pairs_tree_pairwise_consistency.txt", header = T)

sample_dists2 = dists[paste(dat2$sample1, dat2$sample2, sep = "_"), "distance"]

dat2[["dist"]] = sample_dists2

plot(hist(dat2$jaccard, breaks = 100), xlim = c(0, 1.1), ylim = c(0, 250))

plot(hist(dat2$aligned_jaccard, breaks = 100),xlim = c(0, 1.1), ylim = c(0, 250))

dat_sorted <- dat[order(dat[[1]], dat[[2]]), ]
dat2_sorted <- dat2[order(dat2[[1]], dat2[[2]]), ]
# Paired t-test
t.test(dat$jaccard, dat2$jaccard, paired = TRUE)

t.test(dat$aligned_jaccard, dat2$aligned_jaccard, paired = TRUE)

# find biggest difference to look at their cigar strings
# Compute the absolute difference
diffs <- abs(dat$aligned_jaccard - dat2$aligned_jaccard)

# Find the index of the max difference
max_diff_index <- which.max(diffs)

# View the index and the actual values (optional)
max_diff_index
dat[max_diff_index, ]
dat2[max_diff_index, ]

## Plot two histograms on same axis 
h1 = hist(dat$aligned_jaccard, breaks = 100, plot = FALSE)
h2 = hist(dat2$aligned_jaccard, breaks = 100, plot = FALSE)

# Calculate means
mean1 <- mean(dat$aligned_jaccard)
mean2 <- mean(dat2$aligned_jaccard)

# Paired t-test
dat_sorted <- dat[order(dat[[1]], dat[[2]]), ]
dat2_sorted <- dat2[order(dat2[[1]], dat2[[2]]), ]

t_test = t.test(dat$aligned_jaccard, dat2$aligned_jaccard, paired = TRUE)
p_val <- t_test$p.value


plot(h1, col = rgb(0, 0, 1, 0.25), xlim = c(0, 1.1), ylim = c(0, 120), main = "Chr 6 pairwise consistency, 69 samples", xlab = "Aligned Jaccard")

plot(h2, col = rgb(1, 1, 0, 0.25), add = TRUE)

# Add vertical lines for means
abline(v = mean1, col = "orange", lwd = 2, lty = 2)
abline(v = mean2, col = "blue", lwd = 2, lty = 2)

legend("topright",
       legend = c(
         paste0("HOR all pairs tree (mean = ", round(mean1, 3), ")"),
         paste0("Refined Tree (mean = ", round(mean2, 3), ")"),
         paste0("p-value = ", format.pval(p_val, digits = 3, eps = .001))
       ),
       fill = c(rgb(1, 1, 0, 0.5), rgb(0, 0, 1, 0.5), NA),
       border = NA,
       bty = "n",
       cex = 1.2)



