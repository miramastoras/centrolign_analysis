library(ggplot2)

args = commandArgs(trailingOnly=TRUE) # Rscript pairwise_simulations.R input_file chr output_prefix

dat = read.table(args[1])
chr=args[2]

outRecallPNG = paste(args[3], "recall.png", sep="_")
outPrecisionPNG = paste(args[3], "precision.png", sep="_")
outF1PNG = paste(args[3], "f1.png", sep="_")

colnames(dat) = c("case", "aligner", "distance", "truth_matches", "truth_match_rate", "matches", "match_rate", 
                  "mismatches", "mismatch_rate", "recall", "precision")

# add f1 score columns
dat$f1 <- with(dat, ifelse((precision + recall) == 0, NA, 
                           2 * precision * recall / (precision + recall)))

# Calculate mean recall for selected aligners
mean_recalls <- aggregate(recall ~ aligner, 
                          data = dat[dat$aligner %in% c("centrolign", "unialigner", "rama"), ],
                          FUN = mean, na.rm = TRUE)
colnames(mean_recalls)[2] <- "mean_recall"

# Step 2: Create a named vector to relabel aligners with mean recall in legend
label_map <- setNames(
  paste0(mean_recalls$aligner, " (", round(mean_recalls$mean_recall, 2), ")"),
  mean_recalls$aligner
)

# Step 3: Update aligner labels in the dataset
dat$aligner_labeled <- ifelse(dat$aligner %in% names(label_map),
                              label_map[dat$aligner],
                              dat$aligner)

# Step 4: Update mean_recalls to match new labels
mean_recalls$aligner_labeled <- label_map[mean_recalls$aligner]

# Step 5: Plot
recall <- ggplot(data = dat) +
  aes(x = truth_match_rate, y = recall, color = aligner_labeled) +
  geom_point() +
  labs(x = "True match rate", y = "Recall", color = "Aligner (mean recall)") +
  ggtitle(paste("Direct pairwise alignment", chr, sep = " ")) +
  ylim(c(0, 1)) +
  geom_hline(data = mean_recalls, 
             aes(yintercept = mean_recall, color = aligner_labeled),
             linetype = "dashed", show.legend = FALSE) +
  theme(text = element_text(size = 16))

png(outRecallPNG)
print(recall)
dev.off()

mean_precisions <- aggregate(precision ~ aligner, 
                          data = dat[dat$aligner %in% c("centrolign", "unialigner", "rama"), ],
                          FUN = mean, na.rm = TRUE)
colnames(mean_precisions)[2] <- "mean_precision"

# Step 2: Create a named vector to relabel aligners with mean recall in legend
label_map <- setNames(
  paste0(mean_precisions$aligner, " (", round(mean_precisions$mean_precision, 2), ")"),
  mean_precisions$aligner
)

# Step 3: Update aligner labels in the dataset
dat$aligner_labeled <- ifelse(dat$aligner %in% names(label_map),
                              label_map[dat$aligner],
                              dat$aligner)

# Step 4: Update mean_recalls to match new labels
mean_precisions$aligner_labeled <- label_map[mean_precisions$aligner]

# Step 5: Plot
precision <- ggplot(data = dat) +
  aes(x = truth_match_rate, y = precision, color = aligner_labeled) +
  geom_point() +
  labs(x = "True match rate", y = "Precision", color = "Aligner (mean precision)") +
  ggtitle(paste("Direct pairwise alignment", chr, sep = " ")) +
  ylim(c(0, 1)) +
  geom_hline(data = mean_precisions, 
             aes(yintercept = mean_precision, color = aligner_labeled),
             linetype = "dashed", show.legend = FALSE) +
  theme(text = element_text(size = 16))

png(outPrecisionPNG)
print(precision)
dev.off()

### F1 plot
mean_f1s <- aggregate(f1 ~ aligner, 
                             data = dat[dat$aligner %in% c("centrolign", "unialigner", "rama"), ],
                             FUN = mean, na.rm = TRUE)
colnames(mean_f1s)[2] <- "mean_f1"

# Step 2: Create a named vector to relabel aligners with mean recall in legend
label_map <- setNames(
  paste0(mean_f1s$aligner, " (", round(mean_f1s$mean_f1, 2), ")"),
  mean_f1s$aligner
)

# Step 3: Update aligner labels in the dataset
dat$aligner_labeled <- ifelse(dat$aligner %in% names(label_map),
                              label_map[dat$aligner],
                              dat$aligner)

# Step 4: Update mean_recalls to match new labels
mean_f1s$aligner_labeled <- label_map[mean_f1s$aligner]

# Step 5: Plot
f1 <- ggplot(data = dat) +
  aes(x = truth_match_rate, y = f1, color = aligner_labeled) +
  geom_point() +
  labs(x = "True match rate", y = "F1", color = "Aligner (mean F1)") +
  ggtitle(paste("Direct pairwise alignment", chr, sep = " ")) +
  ylim(c(0, 1)) +
  geom_hline(data = mean_f1s, 
             aes(yintercept = mean_f1, color = aligner_labeled),
             linetype = "dashed", show.legend = FALSE) +
  theme(text = element_text(size = 16))

png(outF1PNG)
print(f1)
dev.off()