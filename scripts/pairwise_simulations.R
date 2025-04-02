library(ggplot2)


dat = read.table("/Users/miramastoras/Desktop/pair_chr12_sim_cases_20250331_aln_summary_table.txt")
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
