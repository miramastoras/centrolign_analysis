library(ggplot2)

args = commandArgs(trailingOnly=TRUE) # Rscript pairwise_simulations.R input_file chr output_prefix

dat = read.table(args[1])
chr=args[2]

outRecallPNG = paste(args[3], "recall.png", sep="_")
outPrecisionPNG = paste(args[3], "precision.png", sep="_")

colnames(dat) = c("case", "aligner", "distance", "truth_matches", "truth_match_rate", "matches", "match_rate", 
                  "mismatches", "mismatch_rate", "recall", "precision")

dat = dat[dat$aligner != "wfa",]

recall = (ggplot(data = dat) + aes(x = truth_match_rate, y = recall, color = aligner) 
     + geom_point()
     + labs(x = "True match rate", y = "Recall", color = "Aligner")
     + ggtitle(paste("Direct pairwise alignment", chr, sep=" "))
     + ylim(c(0,1))
)
png(outRecallPNG)
print(recall)
dev.off()

precision = (ggplot(data = dat) + aes(x = truth_match_rate, y = precision, color = aligner) 
          + geom_point()
          + labs(x = "True match rate", y = "Precision", color = "Aligner")
          + ggtitle(paste("Direct pairwise alignment", chr, sep=" "))
          + ylim(c(0,1))
)

png(outPrecisionPNG)
print(precision)
dev.off()

