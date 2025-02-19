library(ggplot2)

setwd("/Users/Jordan/Documents/Research/Pangenomics/Centromeres/working/")

dat = read.table("msa_chrX_summary_tables_max_scale_20240308.txt")
colnames(dat) = c("case", "distance", "truth_matches", "truth_match_rate", "matches", "match_rate", 
                  "mismatches", "mismatch_rate", "recall", "precision")

recall = (ggplot(data = dat) + aes(x = truth_match_rate, y = recall) 
          + geom_point(alpha = 0.25) 
          + labs(x = "True match rate (pairwise upper bound)", y = "Precision")
          + ggtitle("Induced pairwise alignments from MSA")
          + ylim(c(0,1))
)
recall

precision = (ggplot(data = dat) + aes(x = truth_match_rate, y = precision) 
             + geom_point(alpha = 0.25)
             + labs(x = "True match rate (pairwise upper bound)", y = "Recall")
             + ggtitle("Induced pairwise alignments from MSA")
             + ylim(c(0,1))
)
precision


pair_dat = read.table("pairwise_chrX_summary_tables.txt")
colnames(pair_dat) = c("case", "aligner", "distance", "truth_matches", "truth_match_rate", "matches", "match_rate", 
                      "mismatches", "mismatch_rate", "recall", "precision")

joined_dat = rbind(dat[,c("truth_match_rate", "distance", "recall", "precision")],
                   pair_dat[pair_dat$aligner == "centrolign" ,c("truth_match_rate", "distance", "recall", "precision")])
joined_dat[["pairwise"]] = as.factor(c(rep.int("induced", nrow(dat)), 
                                       rep.int("direct", sum(pair_dat$aligner == "centrolign"))))

recall = (ggplot(data = joined_dat) + aes(x = truth_match_rate, y = recall, color = pairwise) 
          + geom_point(alpha = .8) 
          + labs(x = "True match rate", y = "Precision", color = "Pair method")
          + ggtitle("Direct and induced pairwise alignments")
          + ylim(c(0,1))
)
recall

precision = (ggplot(data = joined_dat) + aes(x = truth_match_rate, y = precision, color = pairwise) 
             + geom_point(alpha = .8)
             + labs(x = "True match rate", y = "Recall", color = "Pair method")
             + ggtitle("Direct and induced pairwise alignments")
             + ylim(c(0,1))
)
precision



