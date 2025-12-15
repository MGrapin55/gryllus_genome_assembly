# Repeat Distribution Gpenn/Gfirm 
# Author: Michael Grapin -- Moore Lab Research Technician 

library(tidyverse)
######################################################
##                    Analysis                      ##
######################################################

# Prepocessing the data 
# Bash command 
#awk '{
#   star="FALSE";             # default False
#    if($NF ~ /\*$/){
#        star="TRUE";          # set True if * present
#         sub(/\*$/,"",$NF)
#     }
#     print $0, star
# }' *.out > *.clean.out

# Read in data
Gpenn_rm <- read.table("Gpenn.clean.out",
                 skip = 3,
                 header = FALSE,
                 sep = "",
                 stringsAsFactors = FALSE,
                 fill = TRUE,
                 quote = "")
colnames(Gpenn_rm) <- c(
  "SW_score", "perc_div", "perc_del", "perc_ins",
  "query_sequence", "query_start", "query_end", "query_left",
  "strand",
  "repeat_name", "repeat_classfam",
  "repeat_start", "repeat_end", "repeat_left",
  "ID", "star"
)


# Join with Renamed Scaffold name (Super-Scaffolded###)
Gpenn_key <- read_tsv("Gpenn.asm.key.tsv")
colnames(Gpenn_key)[2] <- "query_sequence"
# vector with key,value pairs
lookup_vec <- setNames(Gpenn_key$Renamed, Gpenn_key$query_sequence)
Gpenn_rm$Renamed <- lookup_vec[Gpenn_rm$query_sequence]
rm(lookup_vec, Gpenn_key)


# Summarize repeats by chromosome 
# First subset the data by chromosomes (super-scaffold1:15) and all unplaced scaffolds (super-scaffold15:)
# Group by Renamed and repeat_classfam 
# combine the sumarized rows together

# Subset chromosomes (super-scaffold1 to super-scaffold15) 
chromosomes <- str_c("Super-Scaffold_", 1:15)

# Create a new column for "Unplaced" scaffolds
Gpenn_rm <- Gpenn_rm %>%
  mutate(
    scaffold_group = ifelse(Renamed %in% chromosomes, Renamed, "Unplaced")
  )

# Group by scaffold and repeat family
summary_df <- Gpenn_rm %>%
  group_by(scaffold_group, repeat_classfam) %>%
  summarise(
    n_repeats = n()
  ) %>%
  ungroup()

# Define the desired order
scaffold_levels <- c(str_c("Super-Scaffold_", 1:15), "Unplaced")

# Convert scaffold_group to factor with these levels
summary_df <- summary_df %>%
  mutate(
    scaffold_group = factor(scaffold_group, levels = scaffold_levels)
  )
# Compute total per scaffold
summary_df <- summary_df %>%
  group_by(scaffold_group) %>%
  mutate(
    total_scaffold = sum(n_repeats),       # total repeats per scaffold
    percent = n_repeats / sum(n_repeats) * 100   # fraction per repeat class
  ) %>%
  ungroup()             

# Make Plot of repeat features per chromosome
# Stacked barplot
ggplot(summary_df, aes(x = scaffold_group, y = percent, fill = repeat_classfam)) +
  geom_bar(stat = "identity", color = "black") +   # add black border
  theme_minimal() +
  labs(
    x = "Chromosomes",
    y = "Percent of Repeats",
    fill = "Repeat Classification",
    title = "Repeat Classes per Scaffold"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(summary_df, aes(x = scaffold_group, y = n_repeats, fill = repeat_classfam)) +
  geom_bar(stat = "identity", color = "black") +   # add black border
  theme_minimal() +
  labs(
    x = "Chromosomes",
    y = "Percent of Repeats",
    fill = "Repeat Classification",
    title = "Repeat Classes per Scaffold"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(summary_df, aes(x = scaffold_group, y = total_scaffold * percent, fill = repeat_classfam)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() +
  labs(
    x = "Scaffold",
    y = "Total Repeats",
    fill = "Repeat Class",
    title = "Repeat Composition per Scaffold"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#####################################################
##                Gfirm                           ##
#####################################################
# Read in data
Gfirm_rm <- read.table("Gfirm.clean.out",
                       skip = 3,
                       header = FALSE,
                       sep = "",
                       stringsAsFactors = FALSE,
                       fill = TRUE,
                       quote = "")
colnames(Gfirm_rm) <- c(
  "SW_score", "perc_div", "perc_del", "perc_ins",
  "query_sequence", "query_start", "query_end", "query_left",
  "strand",
  "repeat_name", "repeat_classfam",
  "repeat_start", "repeat_end", "repeat_left",
  "ID", "star"
)


# Join with Renamed Scaffold name (Super-Scaffolded###)
Gfirm_key <- read_tsv("Gfirm.asm.key.tsv")
colnames(Gfirm_key)[2] <- "query_sequence"
# vector with key,value pairs
lookup_vec <- setNames(Gfirm_key$Renamed, Gfirm_key$query_sequence)
Gfirm_rm$Renamed <- lookup_vec[Gfirm_rm$query_sequence]
rm(lookup_vec, Gfirm_key)


# Summarize repeats by chromosome 
# First subset the data by chromosomes (super-scaffold1:15) and all unplaced scaffolds (super-scaffold15:)
# Group by Renamed and repeat_classfam 
# combine the sumarized rows together

# Subset chromosomes (super-scaffold1 to super-scaffold15) 
chromosomes <- str_c("Super-Scaffold_", 1:15)

# Create a new column for "Unplaced" scaffolds
Gfirm_rm <- Gfirm_rm %>%
  mutate(
    scaffold_group = ifelse(Renamed %in% chromosomes, Renamed, "Unplaced")
  )

# Group by scaffold and repeat family
summary_df2 <- Gfirm_rm %>%
  group_by(scaffold_group, repeat_classfam) %>%
  summarise(
    n_repeats = n()
  ) %>%
  ungroup()

# Define the desired order
scaffold_levels <- c(str_c("Super-Scaffold_", 1:15), "Unplaced")

# Convert scaffold_group to factor with these levels
summary_df2 <- summary_df2 %>%
  mutate(
    scaffold_group = factor(scaffold_group, levels = scaffold_levels)
  )
# Compute total per scaffold
summary_df2 <- summary_df2 %>%
  group_by(scaffold_group) %>%
  mutate(
    total_scaffold = sum(n_repeats),       # total repeats per scaffold
    percent = n_repeats / sum(n_repeats) * 100   # fraction per repeat class
  ) %>%
  ungroup()             

# Make Plot of repeat features per chromosome
# Stacked barplot
ggplot(summary_df2, aes(x = scaffold_group, y = percent, fill = repeat_classfam)) +
  geom_bar(stat = "identity", color = "black") +   # add black border
  theme_minimal() +
  labs(
    x = "Chromosomes",
    y = "Percent of Repeats",
    fill = "Repeat Classification",
    title = "Repeat Classes per Scaffold"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(summary_df2, aes(x = scaffold_group, y = n_repeats, fill = repeat_classfam)) +
  geom_bar(stat = "identity", color = "black") +   # add black border
  theme_minimal() +
  labs(
    x = "Chromosomes",
    y = "Percent of Repeats",
    fill = "Repeat Classification",
    title = "Repeat Classes per Scaffold"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(summary_df, aes(x = scaffold_group, y = total_scaffold * percent, fill = repeat_classfam)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() +
  labs(
    x = "Scaffold",
    y = "Total Repeats",
    fill = "Repeat Class",
    title = "Repeat Composition per Scaffold"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################################################
##            Combined Analysis               ##
################################################

# join by repeat_classfam 
summary_df$species <- "Gpenn"
summary_df2$species <- "Gfirm"

species_joined <- bind_rows(summary_df, summary_df2)

diff_df <- species_joined %>%
  select(species, scaffold_group, repeat_classfam, n_repeats, percent) %>%
  pivot_wider(
    names_from = species,
    values_from = c(n_repeats, percent)
  ) %>%
  mutate(
    diff_n       = n_repeats_Gpenn - n_repeats_Gfirm,
    diff_percent = percent_Gpenn   - percent_Gfirm
  )

library(dplyr)
library(tidyr)
library(pheatmap)

# -------------------------------
# 1. Prepare heatmap matrix
# -------------------------------
heatmap_mat <- diff_df %>%
  select(scaffold_group, repeat_classfam, diff_percent) %>%
  drop_na(diff_percent) %>%        # remove rows where diff_percent is NA
  pivot_wider(
    names_from = scaffold_group,
    values_from = diff_percent,
    values_fill = 0                # optional: fill remaining missing values with 0
  ) %>%
  column_to_rownames("repeat_classfam") %>%
  as.matrix()

scaffold_levels <- c(str_c("Super-Scaffold_", 1:15), "Unplaced")
heatmap_mat <- heatmap_mat[,scaffold_levels]

# -------------------------------
# 2. Define page parameters
# -------------------------------
rows_per_page <- 20
n_rows <- nrow(heatmap_mat)
pages <- ceiling(n_rows / rows_per_page)

# -------------------------------
# 3. Define fixed color scale
# -------------------------------
global_min <- min(heatmap_mat, na.rm = TRUE)
global_max <- max(heatmap_mat, na.rm = TRUE)

n_colors <- 50
my_colors <- colorRampPalette(c("blue", "white", "red"))(n_colors)
my_breaks <- seq(global_min, global_max, length.out = n_colors + 1)

# -------------------------------
# 4. Plot heatmaps page by page
# -------------------------------
pdf("heatmap_pages_fixed_scale.pdf", width = 12, height = 8)

for (i in 1:pages) {
  # Rows for this page
  row_idx <- ((i - 1) * rows_per_page + 1):min(i * rows_per_page, n_rows)
  heatmap_chunk <- heatmap_mat[row_idx, ]
  
  # Plot
  pheatmap(
    heatmap_chunk,
    color = my_colors,
    breaks = my_breaks,         # fixed color scale
    border_color = "black",     # outlines tiles
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 12,
    fontsize_col = 10,
    angle_col = 45,
    cellheight = 20,
    cellwidth = 18,
    main = paste("Change In Percent (Gpenn - Gfirm) per Repeat Family - Page", i)
  )
}

dev.off()



## Make a histogram 
Gpenn_rm <- Gpenn_rm %>% mutate(
  length = abs(query_end - query_start)
)
summary(Gpenn_rm$length)

Gfirm_rm <- Gfirm_rm %>% mutate(
  length = abs(query_end - query_start)
)
summary(Gfirm_rm$length)

Gpenn_rm$species <- "Gpenn"
Gfirm_rm$species <- "Gfirm"

combined <- bind_rows(Gpenn_rm, Gfirm_rm)


library(ggplot2)

ggplot(Gpenn_rm, aes(x = log10(length))) +
  geom_histogram(
    fill = "steelblue",     # fill color of bars
    color = "black"         # border color of bars
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Lenth",
    y = "Count",
    title = "Length of Repeats In Query"
  )


ggplot(Gfirm_rm %>% filter(length > 0), aes(x = length)) +
  geom_histogram(
    fill = "steelblue",
    color = "black",
    bins = 50
  ) +
  scale_x_log10() +
  theme_minimal(base_size = 14) +
  labs(
    x = "Length (log10 scale)",
    y = "Count",
    title = "Length of Repeats In Query"
  )


# Suppose you have a 'species' column in Gpenn_rm
ggplot(combined %>% filter(length > 0), aes(x = length, color = species, fill = species)) +
  geom_density(alpha = 0.4) +     # alpha controls transparency
  scale_x_log10() +               # log10 scale for wide ranges
  theme_minimal(base_size = 14) +
  labs(
    x = "Repeat Length",
    y = "Density",
    title = "Kernel Density of Repeat Lengths by Species",
    color = "Species",
    fill = "Species"
  )


