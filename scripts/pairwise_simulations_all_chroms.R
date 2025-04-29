library(ggplot2)

args = commandArgs(trailingOnly=TRUE) # Rscript pairwise_simulations_all_chroms.R input_file output_prefix

dat = read.table(args[1])

outF1PNG = paste(args[2], "f1.png", sep="_")

colnames(dat) = c("case", "aligner", "distance", "truth_matches", "truth_match_rate", "matches", "match_rate", 
                  "mismatches", "mismatch_rate", "recall", "precision","chr")

# add f1 score columns
dat$f1 <- with(dat, ifelse((precision + recall) == 0, NA, 
                           2 * precision * recall / (precision + recall)))

# Step 1: Sort chromosomes by numeric order
dat$chr <- factor(dat$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))

# Calculate average F1 score per aligner and chromosome
avg_f1_per_align_chr <- aggregate(f1 ~ chr + aligner, data = dat, FUN = mean, na.rm = TRUE)

blue_shades <- c("#AE8799", "#0F8B8D", "#EC9A29")  # Modify as needed

# Step 2: Plot the average F1 score per aligner per chromosome
f1_plot <- ggplot(avg_f1_per_align_chr, aes(x = chr, y = f1, fill = aligner)) +
  geom_bar(stat = "identity", position = "dodge") +  # Separate bars by aligner
  labs(x = "Chromosome", y = "Average F1 Score", title = "Average F1 Score per Aligner per Chromosome") +
  scale_fill_manual(values = blue_shades) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.background = element_blank(),  # Remove the background panel color
        plot.background = element_blank())   # Remove background around the plot

# Save the plot with a wider aspect ratio
png(outF1PNG, width = 12, height = 6, units = "in", res = 300)  # Save as PNG
print(f1_plot)  # Render the plot
dev.off()  # Close the device to save the file
