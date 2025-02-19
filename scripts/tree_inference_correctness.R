library(ggplot2)
library(ggbeeswarm)

setwd("/Users/Jordan/Documents/Research/Pangenomics/Centromeres/working/")
dat = read.table("msa_chrX_tree_comparisons.tsv")
colnames(dat) = c("case", "subtree_height", "subtree_size", "correct")
dat$correct = as.factor(ifelse(dat$correct, "Correct", "Incorrect"))

mean(dat$correct == "Correct")
table(dat[,c("case", "correct")])

p = (ggplot(data = dat)
    + aes(x = correct, y = subtree_height, color = correct)
    + labs(x = "Inferred partitions", y = "Internal node height (generations)")
    + ggtitle("Inferred partitions of leaves in neighbor-joining trees from MSAs")
    + theme(legend.position = "none")
    + geom_quasirandom()
)
p
