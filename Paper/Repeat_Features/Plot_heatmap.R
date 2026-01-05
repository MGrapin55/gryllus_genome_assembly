library(dplyr)
library(tidyr)
library(pheatmap)
library(stringr)

# Reading in the data with a function read_rm_out() from functions.R 
p <- read_rm_out("Gpenn.clean.out")
f <- read_rm_out("Gfirm.clean.out")

# Add a species column 
p$species <- "Gpenn"
f$species <- "Gfirm"

# Combine dataframes
combined <- bind_rows(p, f)

# Set factor levels 
combined <- combined %>%
  mutate(species = factor(species, levels = c("Gpenn", "Gfirm")))

# -------------------------------
# 1. Summarize differences (counts)
# -------------------------------
summary_feature <- combined %>%
  count(species, group, class_family, name = "n")

# Pivot to wide format
diff_df <- summary_feature %>%
  pivot_wider(
    names_from = species,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    diff_n = Gpenn - Gfirm
  )

# -------------------------------
# 2. Prepare heatmap matrix
# -------------------------------
heatmap_mat <- diff_df %>%
  select(group, class_family, diff_n) %>%
  pivot_wider(
    names_from = group,
    values_from = diff_n,
    values_fill = 0
  ) %>%
  column_to_rownames("class_family") %>%
  as.matrix()

# Reorder columns
scaffold_levels <- c("X_chr", paste0("chr_", 1:14), "Unplaced")
existing_cols <- intersect(scaffold_levels, colnames(heatmap_mat))
heatmap_mat <- heatmap_mat[, existing_cols, drop = FALSE]

# -------------------------------
# 3. Page parameters
# -------------------------------
rows_per_page <- 20
n_rows <- nrow(heatmap_mat)
pages <- ceiling(n_rows / rows_per_page)

# -------------------------------
# 4. Color scale
# -------------------------------
global_min <- min(heatmap_mat, na.rm = TRUE)
global_max <- max(heatmap_mat, na.rm = TRUE)

n_colors <- 50
my_colors <- colorRampPalette(c("blue", "white", "red"))(n_colors)
my_breaks <- seq(global_min, global_max, length.out = n_colors + 1)

# -------------------------------
# 5. Plot heatmaps page by page
# -------------------------------
pdf("heatmap_Repeat_counts.pdf", width = 12, height = 8)

for (i in 1:pages) {
  row_idx <- ((i - 1) * rows_per_page + 1):min(i * rows_per_page, n_rows)
  heatmap_chunk <- heatmap_mat[row_idx, , drop = FALSE]
  
  pheatmap(
    heatmap_chunk,
    color = my_colors,
    breaks = my_breaks,
    border_color = "black",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 12,
    fontsize_col = 10,
    angle_col = 45,
    cellheight = 20,
    cellwidth = 18,
    main = paste("Change in Count (Gpenn - Gfirm) per Repeat Family - Page", i)
  )
}

dev.off()



########################################################################################################
# Repeating with Percentages

# Set factor levels (Note: Moved this down to where 'summary_feature' is actually created)

# -------------------------------
# 1. Summarize differences (Percentages)
# -------------------------------

# First, get raw counts
summary_counts <- combined %>%
  count(species, group, class_family, name = "n") 

# Second, calculate percentages
# logic: (Count of specific repeat / Total repeats in that species) * 100
summary_feature <- summary_counts %>%
  group_by(species, group) %>%
  mutate(total_species_n = sum(n)) %>%
  mutate(pct = (n / total_species_n) * 100) %>%
  ungroup() %>%
  mutate(species = factor(species, levels = c("Gpenn", "Gfirm")))

# Pivot to wide format using Percentage
diff_df <- summary_feature %>%
  select(species, group, class_family, pct) %>% # Select pct, ignore raw 'n'
  pivot_wider(
    names_from = species,
    values_from = pct,
    values_fill = 0
  ) %>%
  mutate(
    diff_pct = Gpenn - Gfirm
  )

# -------------------------------
# 2. Prepare heatmap matrix
# -------------------------------
# We now select 'diff_pct' instead of 'diff_n'
heatmap_mat <- diff_df %>%
  select(group, class_family, diff_pct) %>%
  pivot_wider(
    names_from = group,
    values_from = diff_pct,
    values_fill = 0
  ) %>%
  column_to_rownames("class_family") %>%
  as.matrix()

# Reorder columns
scaffold_levels <- c("X_chr", paste0("chr_", 1:14), "Unplaced")
existing_cols <- intersect(scaffold_levels, colnames(heatmap_mat))
heatmap_mat <- heatmap_mat[, existing_cols, drop = FALSE]

# -------------------------------
# 3. Page parameters
# -------------------------------
rows_per_page <- 20
n_rows <- nrow(heatmap_mat)
pages <- ceiling(n_rows / rows_per_page)

# -------------------------------
# 4. Color scale
# -------------------------------
global_min <- min(heatmap_mat, na.rm = TRUE)
global_max <- max(heatmap_mat, na.rm = TRUE)

# Ensure 0 is centered in the color palette if you want a diverging scale (Blue=Negative, Red=Positive)
limit <- max(abs(global_min), abs(global_max))
my_breaks <- seq(-limit, limit, length.out = 51)
n_colors <- length(my_breaks) - 1
my_colors <- colorRampPalette(c("blue", "white", "red"))(n_colors)

# -------------------------------
# 5. Plot heatmaps page by page
# -------------------------------
pdf("heatmap_Repeat_Percentages.pdf", width = 12, height = 8)

for (i in 1:pages) {
  row_idx <- ((i - 1) * rows_per_page + 1):min(i * rows_per_page, n_rows)
  heatmap_chunk <- heatmap_mat[row_idx, , drop = FALSE]
  
  pheatmap(
    heatmap_chunk,
    color = my_colors,
    breaks = my_breaks,
    border_color = "black",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 12,
    fontsize_col = 10,
    angle_col = 45,
    cellheight = 20,
    cellwidth = 18,
    # Updated title to reflect unit change
    main = paste("Difference in Abundance (%) (Gpenn - Gfirm) - Page", i) 
  )
}

dev.off()
