dat = read.table("/Users/miramastoras/Desktop/pairwise_consistency.txt", header = T)
dists = read.csv("/Users/miramastoras/Desktop/pairwise_distance_excl_HG00741.1.csv", header = T)

key1 = paste(dists$sample1, dists$sample2, sep = "_")
key2 = paste(dists$sample2, dists$sample1, sep = "_")

dists = rbind(dists, dists)
row.names(dists) = c(key1, key2)

sample_dists = dists[paste(dat$sample1, dat$sample2, sep = "_"), "distance"]

dat[["dist"]] = sample_dists

#plot(density(dat$jaccard))

plot(hist(dat$jaccard, breaks = 100))

plot(dat$dist, dat$jaccard, pch = 19, col = alpha("black", 0.1), 
     xlab = "Patristic distance", ylab = "Jaccard similarity", main = "All pairs")

plot(hist(dat$aligned_jaccard, breaks = 100))

plot(dat$dist, dat$aligned_jaccard, pch = 19, col = alpha("black", 0.1), 
     xlab = "Patristic distance", ylab = "Jaccard similarity", main = "Only aligned pairs")


# subtree_samples = c("HG01358.1", "HG01071.1", "HG01361.1", "HG03540.2", "HG02055.1", "HG03540.1",
#                     "HG00438.1", "HG01978.1", "HG00438.2", "HG01928.2")
# 
# aln_in_subtree = (dat$sample1 %in% subtree_samples) & (dat$sample2 %in% subtree_samples)
# dat_subtree = dat[aln_in_subtree, ]
# 
# plot(dat_subtree$dist, dat_subtree$aligned_jaccard, pch = 19, col = alpha("black", 0.1), 
#      xlab = "Patristic distance", ylab = "Jaccard similarity", main = "Only aligned pairs")
