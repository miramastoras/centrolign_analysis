library(ggplot2)

# read in data for refined tree 
dat = read.table("/Users/miramastoras/Desktop/chr12_tree_consistency_79/chr12_79_refined_tree_pairwise_consistency.txt", header = T)
dists = read.csv("/Users/miramastoras/Desktop/chr12_tree_consistency_79/pairwise_distance_excl_HG00741.1.csv", header = T)

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
dat2 = read.table("/Users/miramastoras/Desktop/chr12_tree_consistency_79/chr12_79_HOR_all_pairs_NJ_pairwise_consistency.txt", header = T)

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
