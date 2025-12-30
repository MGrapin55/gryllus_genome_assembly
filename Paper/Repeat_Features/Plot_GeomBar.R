
# Plot Number of repetitive features on a bar graph
# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(ggplot2)
library(scales)
#########################################
##          Preprocessing              ##
#########################################
# Data cleaned (Make * a new columb)
Gpenn_rm <- read.table("Gpenn.clean.out",
                       skip = 3,
                       header = FALSE,
                       sep = "",
                       stringsAsFactors = FALSE,
                       fill = TRUE,
                       quote = "")

Gfirm_rm <- read.table("Gfirm.clean.out",
                       skip = 3,
                       header = FALSE,
                       sep = "",
                       stringsAsFactors = FALSE,
                       fill = TRUE,
                       quote = "")

# Column Names
columns <- c(
  "SW_score", "perc_div", "perc_del", "perc_ins",
  "query_sequence", "query_start", "query_end", "query_left",
  "strand",
  "repeat_name", "repeat_classfam",
  "repeat_start", "repeat_end", "repeat_left",
  "ID", "star"
)

colnames(Gpenn_rm) <- columns
colnames(Gfirm_rm) <- columns


# Add a species column 
Gpenn_rm$species <- "Gpenn"
Gfirm_rm$species <- "Gfirm"

# Combine Rows
combined <- bind_rows(Gpenn_rm, Gfirm_rm)

# Want to calculated the center of each repeat feature
# feature_center = abs(query_start + query_end / 2) (its the arthimetic mean)
combined <- combined %>% mutate(
  repeat_center = abs(query_start + query_end / 2)
)
# Now we want to split up repeat_classfam
combined <- combined %>%
  mutate(
    repeat_class = sub("/.*", "", repeat_classfam)
  )

# Subset chromosomes
chromosomes <- str_c("chr_", 1:15)

# Define desired order: chr_X first, then chr_1:14, then Unplaced
group_levels <- c("chr_X", str_c("chr_", 1:14), "Unplaced")

combined <- combined %>%
  mutate(
    group = case_when(
      query_sequence == "chr_1" ~ "chr_X",
      query_sequence %in% paste0("chr_", 2:15) ~ paste0("chr_", as.numeric(str_remove(query_sequence, "chr_")) - 1),
      TRUE ~ "Unplaced"
    ),
    group = factor(group, levels = group_levels)
  )

# Summarize by repeat_classfam
summary_feature <- combined %>%
  count(species, group, repeat_classfam, name = "n")

summary_feature <- summary_feature %>%
  mutate(species = factor(species, levels = c("Gpenn", "Gfirm")))

# Open PDF
pdf("class_bar.pdf", width = 8, height = 6)

classes <- unique(summary_feature$repeat_classfam)

for (c in classes) {
  
  df <- summary_feature %>%
    filter(repeat_classfam == c)
  
  if (nrow(df) == 0) next
  
  p <- ggplot(df, aes(x = group, y = n, fill = species)) +
    geom_col(
      position = position_dodge(width = 0.9),
      color = "black",
      linewidth = 0.4
    ) +
    scale_fill_manual(values = c("#618B4A", "#AFBC88")) +
    theme_classic() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      y = "Count of Repeats",
      x = "Chromosome",
      title = paste("Repeat Class:", c)
    )
  
  print(p)
}
# Close PDF
dev.off()


# Have these could summarized by chr 
# Summarize by feature per chromosome
summary_chr <- combined %>%
  group_by(species, group) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )
summary_chr <- summary_chr %>%
  mutate(species = factor(species, levels = c("Gpenn", "Gfirm")))

pdf("TotalPerRepeat_bar.pdf", width = 8, height = 6)

ggplot(summary_chr, aes(x = group, y = n, fill = species)) +
  geom_col(
    position = position_dodge(width = 0.9),
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_manual(values = c("#618B4A", "#AFBC88")) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    y = "Count of Repeats",
    x = "Chromosome",
    title = paste("Total Counts of Repeat Content")
  ) +
  scale_y_continuous(labels = comma)
dev.off()



# Summary by just the class
# Summarize by repeat_classfam
summary_class <- combined %>%
  count(species, group, repeat_class, name = "n")

summary_class <- summary_class %>%
  mutate(species = factor(species, levels = c("Gpenn", "Gfirm")))

# Open PDF
pdf("class_bar.pdf", width = 8, height = 6)

classes <- unique(summary_class$repeat_class)

for (c in classes) {
  
  df <- summary_class %>%
    filter(repeat_class == c)
  
  if (nrow(df) == 0) next
  
  p <- ggplot(df, aes(x = group, y = n, fill = species)) +
    geom_col(
      position = position_dodge(width = 0.9),
      color = "black",
      linewidth = 0.4
    ) +
    scale_fill_manual(values = c("#618B4A", "#AFBC88")) +
    theme_classic() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      y = "Count of Repeats",
      x = "Chromosome",
      title = paste("Repeat Class:", c)
    )
  
  print(p)
}
# Close PDF
dev.off()
