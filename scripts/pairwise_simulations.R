library(ggplot2)

setwd("/Users/Jordan/Documents/Research/Pangenomics/Centromeres/working/")

dat = read.table("pairwise_chrX_summary_tables.txt")
colnames(dat) = c("case", "aligner", "distance", "truth_matches", "truth_match_rate", "matches", "match_rate", 
                  "mismatches", "mismatch_rate", "recall", "precision")

dat = dat[dat$aligner != "wfa",]

recall = (ggplot(data = dat) + aes(x = truth_match_rate, y = recall, color = aligner) 
     + geom_point()
     + labs(x = "True match rate", y = "Recall", color = "Aligner")
     + ggtitle("Direct pairwise alignment")
     + ylim(c(0,1))
)
recall

precision = (ggplot(data = dat) + aes(x = truth_match_rate, y = precision, color = aligner) 
          + geom_point()
          + labs(x = "True match rate", y = "Precision", color = "Aligner")
          + ggtitle("Direct pairwise alignment")
          + ylim(c(0,1))
)
precision
