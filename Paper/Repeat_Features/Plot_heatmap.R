library(dplyr)
library(tidyr)
library(pheatmap)
library(stringr)

# -------------------------------
# 1. Summarize differences
# -------------------------------
summary_feature <- combined %>%
  count(species, group, repeat_classfam, name = "n")

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
  select(group, repeat_classfam, diff_n) %>%
  pivot_wider(
    names_from = group,
    values_from = diff_n,
    values_fill = 0
  ) %>%
  column_to_rownames("repeat_classfam") %>%
  as.matrix()

# Reorder columns
scaffold_levels <- c("chr_X", paste0("chr_", 1:14), "Unplaced")
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

