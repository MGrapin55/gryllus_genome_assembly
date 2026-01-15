# R script from https://github.com/kataokaklab/Gryllus_assimilis_genome/blob/main/figures/fig5/plot_coverage.R
# Adpated By Michael Grapin

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)

# Global Parameters
wkdir = "~/Downloads/MOORE_LAB_UNL/GRYLLUS_GENOME_ASSEMBLY/GIT_REPO/gryllus_genome_assembly/Paper/Validate_X"
result_dir = "./results"
coverage_file = "GFIRM_LZ.regions.bed"
key_file = "../Repeat_Features/Gfirm.asm.key.tsv"

setwd(wkdir)

# Load the data
cov <- read.csv(coverage_file, header = FALSE, sep = "\t")
key <- read_tsv(key_file)

# Create a vector of scaffold names
scaffolds  <- key[[1]][1:15]
old_labels <- key[[2]][1:15] 

new_labels <- case_when(
  old_labels == "chr_1" ~ "X_chr",
  old_labels %in% paste0("chr_", 2:15) ~ 
    paste0("chr_", as.numeric(str_remove(old_labels, "chr_")) - 1),
  TRUE ~ "Unplaced"
)

new_order <-  c("X_chr", str_c("chr_", 1:14))

# Filter rows corresponding to specific scaffolds
cov_filtered <- filter(cov, V1 %in% scaffolds)

# Update the V1 column with new labels and set the order
cov_filtered$V1 <- factor(cov_filtered$V1, 
                               levels = scaffolds,  # Original scaffold order
                               labels = new_labels)  # Replace with new labels
cov_filtered$V1 <- factor(cov_filtered$V1, levels = new_order)  # Set the desired order

# Remove outliers using the interquartile range (IQR)
cov_filtered <- cov_filtered %>%
  group_by(V1) %>%
  mutate(
    Q1 = quantile(V4, 0.25),  # First quartile
    Q3 = quantile(V4, 0.75),  # Third quartile
    IQR = Q3 - Q1  # Interquartile range
  ) %>%
  filter(V4 >= (Q1 - 1.5 * IQR) & V4 <= (Q3 + 1.5 * IQR)) %>%  # Keep values within the IQR range
  ungroup() %>%
  select(-Q1, -Q3, -IQR)  # Remove temporary columns for quartiles and IQR

# Get the maximum value after filtering outliers
max_value <- max(cov_filtered$V4)

# Create a boxplot
plot <- ggplot(cov_filtered, aes(x = V1, y = V4)) +
  geom_boxplot(width = 0.7, lwd = 0.5, color = "black", fill = NA, outlier.shape = NA) +  # Customize boxplot appearance
  ylim(0, max_value) +  # Set Y-axis range from 0 to the maximum value
  labs(x = "Chromosomes", y = "Coverage") +  # Add axis labels
  theme(
    text = element_text(family = "Helvetica", size = 20),  # Set font and size
    axis.text = element_text(size = 20, face = "bold", colour = "black"),  # Customize axis text
    axis.title = element_text(size = 20, face = "bold", colour = "black"),  # Customize axis titles
    axis.line = element_line(linewidth = 1.2, colour = "black"),  # Customize axis lines
    panel.background = element_blank(),  # Remove background
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate X-axis labels
    legend.position = "none"  # Hide legend
  )

# Make results directory
dir.create(result_dir)

# Save the plot in SVG format
ggsave(result_dir, plot = plot, device = "svg")

# Save the plot in PNG format
ggsave(result_dir, plot = plot, device = "png")

# Save the plot in PDF format
ggsave(result_dir, plot = plot, device = "pdf")
