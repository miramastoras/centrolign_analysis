library(ggplot2)

args = commandArgs(trailingOnly=TRUE) # Rscript pairwise_simulations.R input_file chr output_prefix

dat = read.table(args[1])
chr=args[2]

outRecallPNG = paste(args[3], "recall.png", sep="_")
outPrecisionPNG = paste(args[3], "precision.png", sep="_")

colnames(dat) = c("case", "aligner", "distance", "truth_matches", "truth_match_rate", "matches", "match_rate", 
                  "mismatches", "mismatch_rate", "recall", "precision")

# Calculate mean recall for selected aligners
mean_recalls <- aggregate(recall ~ aligner, 
                          data = dat[dat$aligner %in% c("centrolign", "unialigner", "rama"), ],
                          FUN = mean, na.rm = TRUE)
colnames(mean_recalls)[2] <- "mean_recall"

recall = (ggplot(data = dat) + aes(x = truth_match_rate, y = recall, color = aligner) 
     + geom_point()
     + labs(x = "True match rate", y = "Recall", color = "Aligner")
     + ggtitle(paste("Direct pairwise alignment", chr, sep=" "))
     + ylim(c(0,1)) 
     + geom_hline(data = mean_recalls, aes(yintercept = mean_recall, color = aligner),
                  linetype = "dashed", show.legend = FALSE) 
     + geom_text(data = mean_recalls, aes(x = Inf, y = mean_recall, label = round(mean_recall, 2)),
                 hjust = 1.1, vjust = -0.5, color = "black", inherit.aes = FALSE)
)
png(outRecallPNG)
print(recall)
dev.off()

mean_precisions <- aggregate(precision ~ aligner, 
                          data = dat[dat$aligner %in% c("centrolign", "unialigner", "rama"), ],
                          FUN = mean, na.rm = TRUE)
colnames(mean_precisions)[2] <- "mean_precision"

precision = (ggplot(data = dat) + aes(x = truth_match_rate, y = precision, color = aligner) 
          + geom_point()
          + labs(x = "True match rate", y = "Precision", color = "Aligner")
          + ggtitle(paste("Direct pairwise alignment", chr, sep=" "))
          + ylim(c(0,1))
          + geom_hline(data = mean_precisions, aes(yintercept = mean_precision, color = aligner),
                       linetype = "dashed", show.legend = FALSE) 
          + geom_text(data = mean_precisions, aes(x = Inf, y = mean_precision, label = round(mean_precision, 2)),
                      hjust = 1.1, vjust = -0.5, color = "black", inherit.aes = FALSE)
)

png(outPrecisionPNG)
print(precision)
dev.off()

