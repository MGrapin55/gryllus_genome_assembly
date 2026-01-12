setwd("~/Downloads/MOORE_LAB_UNL/GRYLLUS_GENOME_ASSEMBLY/GIT_REPO/gryllus_genome_assembly/Paper/Repeat_Features")

# Load some handy functions in
source(file = "functions.R")

# Reading in the data with a function read_rm_out() from functions.R 
p <- read_rm_out("Gpenn.clean.out")
f <- read_rm_out("Gfirm.clean.out")

# Add a species column 
p$species <- "Gpenn"
f$species <- "Gfirm"

# Combine dataframes
rm_df <- bind_rows(p, f)

# Set factor levels and get new columns
rm_df <- rm_df %>%
  mutate(species = factor(species, levels = c("Gpenn", "Gfirm")),
         midpoint = abs((begin + end) / 2),
         midpoint_mb = midpoint / 1e6)
#####################################################################
# Shared parameters
window_size <- 1e6  # 1 Mb windows
chr <- c("X_chr", str_c("chr_", 1:14))

# Plot colors 
colors <- c("Gpenn" = "#618B4A", "Gfirm" = "#AFBC88")

####################################################################
# Remove class_family values that aren't in both 
rm_df <- rm_df %>%
  group_by(class_family) %>%
  filter(n_distinct(species) == 2) %>%
  ungroup()

# formatting the data 
df_elem_density <- rm_df %>% 
  filter(group %in% chr) %>%
  mutate(window = floor(midpoint / window_size)) %>%
  group_by(species, window, group, class_family) %>%
  summarize(
    pos_mb = mean(midpoint_mb),
    density_per_mb = n() / (window_size / 1e6),
    .groups = "drop"
  )

shared_units <- df_elem_density %>%
  distinct(species, group, window, class_family) %>%
  count(group, window, class_family) %>%
  filter(n == 2) %>%          # present in both species
  select(group, window, class_family)


df_diff <- df_elem_density %>%
  semi_join(shared_units, by = c("group", "window", "class_family")) %>%
  pivot_wider(
    id_cols    = c(group, window, class_family),
    names_from = species,
    values_from = density_per_mb
  ) %>%
  mutate(density_diff = Gfirm - Gpenn)

df_diff <- df_diff %>% mutate(
  class = sub("/.*", "", class_family)
)


# the plot
# x = genomic positions (1 mb windows)
# y = net change in species counts by window

ggplot(df_diff, aes(x = window, y = density_diff, color = class)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linewidth = 0.6)



df_diff <- df_diff %>% arrange(group, window) %>% 
  mutate(geneome_length_of_windown_MB = row_number()- 1)

groups <- unique(df_diff$group)

for (g in groups) {
  p <- ggplot(
    df_diff %>% filter(group == g),
    aes(x = window, y = density_diff, color = class)
  ) +
    geom_point() +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.6) +
    labs(title = paste("Net Variation in Repeat Features on:",g), 
         subtitle = "Net = Gfirm.feature - Gpenn.feature") 
  
  print(p)
}



# Define chromosome order
chrom_order <- c("X_chr", paste0("chr_", 1:14))  

# Compute max window per chromosome
chr_lengths <- df_diff %>%
  group_by(group) %>%
  summarize(max_window = max(window), .groups = "drop") %>%
  mutate(group = factor(group, levels = chrom_order)) %>%
  arrange(group)

chr_lengths <- chr_lengths %>%
  mutate(
    chr_start = lag(cumsum(max_window + 1), default = 0)  # +1 because windows start at 0
  )

df_diff <- df_diff %>%
  left_join(chr_lengths %>% select(group, chr_start), by = "group") %>%
  mutate(
    genome_pos = chr_start + window  # continuous genome-wide coordinate
  )



ggplot(df_diff, aes(x = genome_pos, y = density_diff, color = class)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linewidth = 0.6) +
  scale_x_continuous(
    breaks = chr_lengths$chr_start + chr_lengths$max_window / 2,
    labels = chr_lengths$group
  ) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Genome (windows 1 Mb)", y = "Density difference", 
       title = "Net Density Variation", 
       subtitle = ("Net = Gfirm.elementDensity - Gpenn.elementDensity"))

+ facet_wrap(~class)
## this is useful and we can break this down even more among class_family to see how the density changes
## This will look nice to then break out to smaller subgroups that we can test for with enrichment tests or 
## chi-squared. 


# Plots every combination at the Class level
groups   <- unique(df_diff$group)
features <- unique(df_diff$class)

pdf("density_diff_by_group_and_class.pdf", width = 10, height = 4)

for (g in groups) {
  for (f in features) {
    
    df_sub <- df_diff %>%
      filter(group == g, class == f)
    
    if (nrow(df_sub) == 0) next   # skip empty combinations
    
    p <- ggplot(df_sub, aes(x = genome_pos, y = density_diff)) +
      geom_point() +
      geom_hline(yintercept = 0, color = "red", linewidth = 0.6) +
      scale_x_continuous(
        breaks = chr_lengths$chr_start + chr_lengths$max_window / 2,
        labels = chr_lengths$group
      ) +
      theme_bw() +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      labs(
        x = "Genomic Window (1 Mb)",
        y = "Density difference",
        title = paste("Repeat Elements:", f, "along", g)
      )
    
    print(p)   # one page per plot
  }
}

dev.off()
# This is useful to look at and we can make a figure that has the whole genome and then branches out to
# specific parts of the genome. 



# Vizually this will be easter to see in a heatmap.
library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)

# ----------------------------
# Parameters
# ----------------------------
rows_per_page <- 20

# ----------------------------
# Open PDF
# ----------------------------
pdf("density_diff_heatmaps_by_chromosome.pdf", width = 11, height = 8)

# ----------------------------
# Loop over chromosomes
# ----------------------------
for (g in groups) {
  
  df_chr <- df_diff %>%
    filter(group == g) %>%
    arrange(class_family, window)
  
  if (nrow(df_chr) == 0) next
  
  # ----------------------------
  # Build class_family × window matrix
  # ----------------------------
  heat_mat <- df_chr %>%
    select(class_family, window, density_diff) %>%
    pivot_wider(
      names_from  = window,
      values_from = density_diff
    ) %>%
    column_to_rownames("class_family") %>%
    as.matrix()
  
  # Ensure windows are ordered
  heat_mat <- heat_mat[, order(as.numeric(colnames(heat_mat))), drop = FALSE]
  
  # Color scale limits
  lim <- max(abs(heat_mat), na.rm = TRUE)
  
  # ----------------------------
  # Split rows INTO PAGES (FIX)
  # ----------------------------
  row_ids <- rownames(heat_mat)
  
  class_chunks <- split(
    row_ids,
    ceiling(seq_along(row_ids) / rows_per_page)
  )
  
  # ----------------------------
  # Plot one heatmap per page
  # ----------------------------
  for (i in seq_along(class_chunks)) {
    
    sub_mat <- heat_mat[class_chunks[[i]], , drop = FALSE]
    
    pheatmap(
      sub_mat,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      scale = "none",
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-lim, lim, length.out = 101),
      show_colnames = FALSE,
      fontsize_row = 8,
      main = paste0(
        "Density difference (Gfirm − Gpenn)\n",
        "Chromosome: ", g,
        " | Rows ", ((i - 1) * rows_per_page + 1),
        "–", min(i * rows_per_page, nrow(heat_mat))
      )
    )
  }
}

# ----------------------------
# Close PDF
# ----------------------------
dev.off()
## This is useful to have but it is hard as a vizualization to actually have in a figure. We can always
## show specific rows that we think are biologically relevant. 




