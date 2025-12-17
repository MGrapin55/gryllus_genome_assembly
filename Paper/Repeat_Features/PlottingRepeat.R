# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(colorspace)

# ------------------------------------------------------------
# 1) Extract repeat superclass (before "/")
# ------------------------------------------------------------
Gpenn_summary <- Gpenn_summary %>%
  mutate(
    repeat_class = sub("/.*", "", repeat_classfam)
  )

# ------------------------------------------------------------
# 2) Base color per superclass
# ------------------------------------------------------------
classes <- sort(unique(Gpenn_summary$repeat_class))

base_cols <- setNames(
  qualitative_hcl(length(classes), palette = "Dark 3"),
  classes
)
# ------------------------------------------------------------
# 3) Build palette safely (group-wise)
# ------------------------------------------------------------
pal_df <- Gpenn_summary %>%
  distinct(repeat_class, repeat_classfam) %>%
  group_by(repeat_class) %>%
  group_modify(~ {
    n <- nrow(.x)
    base <- base_cols[.y$repeat_class]
    
    .x %>%
      mutate(
        fill_col = if (n == 1) {
          base
        } else {
          lighten(base, amount = seq(0.4, -0.3, length.out = n))
        }
      )
  }) %>%
  ungroup()

pal <- setNames(pal_df$fill_col, pal_df$repeat_classfam)
# ------------------------------------------------------------
# 4) Order factor levels to match palette (clean legend)
# ------------------------------------------------------------
Gpenn_summary$repeat_classfam <- factor(
  Gpenn_summary$repeat_classfam,
  levels = names(pal)
)

# ------------------------------------------------------------
# 5) Plot
# ------------------------------------------------------------
P <-  ggplot(
      Gpenn_summary,
      aes(x = scaffold_group, y = percent, fill = repeat_classfam)
    ) +
      geom_bar(stat = "identity", color = "black") +
      scale_fill_manual(values = pal) +
      theme_minimal() +
      labs(
        x = "Chromosomes",
        y = "Percent of Repeats",
        fill = "Repeat Classification",
        title = "G.pennsylvanicus"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )





# ------------------------------------------------------------
# Reuse an existing palette on a new species
# Missing classifications are automatically skipped
# ------------------------------------------------------------

# pal : named vector (repeat_classfam -> color) from species A
# Gfirm_summary  : summary dataframe for Gfirm

# Keep only colors that exist in the new species
pal_new <- pal[names(pal) %in% Gfirm_summary$repeat_classfam]

# Lock factor levels for consistent legend ordering
Gfirm_summary$repeat_classfam <- factor(
  Gfirm_summary$repeat_classfam,
  levels = intersect(names(pal), Gfirm_summary$repeat_classfam)
)

# Plot
F <- ggplot(
      Gfirm_summary,
      aes(x = scaffold_group, y = percent, fill = repeat_classfam)
    ) +
      geom_bar(stat = "identity", color = "black") +
      scale_fill_manual(values = pal_new) +
      theme_minimal() +
      labs(
        x = "Chromosomes",
        y = "Percent of Repeats",
        fill = "Repeat Classification",
        title = "G.firmus"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

# See what different names there are 
setdiff(names(pal), Gfirm_summary$repeat_classfam)





## Plot on the Same plot 


library(ggplot2)

repeats <- combine()

ggplot(test, aes(x=repeat_center, y=species, color=species, point_color=species, fill=species)) +
 geom_density_ridges(
   jittered_points=TRUE, scale = .95, rel_min_height = .01,
   point_shape = "|", point_size = 1, size = 0.25,
   position = position_points_jitter(height = 0), alpha = 0.5
 ) +
 scale_y_discrete(expand = c(.01, 0)) +
 scale_x_continuous(expand = c(0, 0), name = "bp") + theme(legend.position = "none")



# plotting for a density plot
# First have to combine Gpenn and Gfirm repeat masker into one dataframe

# Add a species column 
Gpenn_rm$species <- "Gpenn"
Gfirm_rm$species <- "Gfirm"

# Combine Rows
combined <- bind_rows(Gpenn_rm, Gfirm_rm)

# Want to caculated the center of each repeat feature
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

# Define desired order
group_levels <- c(chromosomes, "Unplaced")

combined <- combined %>%
  mutate(
    group = if_else(
      query_sequence %in% chromosomes,
      query_sequence,
      "Unplaced"
    ),
    group = factor(group, levels = group_levels)
  )

test <- combined %>%
  filter(group == "chr_1", repeat_classfam == "DNA/hAT")

# Determine the range of your data
x_min <- min(test$repeat_center, na.rm = TRUE)
x_max <- max(test$repeat_center, na.rm = TRUE)

ggplot(test, aes(x = repeat_center/1e6, y = species, fill = species)) +
  geom_density_ridges(
    aes(height = after_stat(density)),
    scale = 5,
    alpha = 0.7,
    linewidth = 0.8,                     # ridge outline thickness

  ) +
  scale_fill_manual(
    values = c("#AFBC88", "#618B4A"),
    labels = c(expression(italic("G. firmus")), 
               expression(italic("G. pennsylvanicus")))
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    name = "Genomic position (Mb)",
    breaks = seq(0, max(test$repeat_center)/1e6, by = 50),  # 0.03 Mb = 30 kb
    labels = scales::number_format(accuracy = 0.1) 
  ) +
  theme_ridges(
    grid = FALSE,
    center_axis_labels = TRUE
  ) +
  theme(
    legend.position = "none"
  ) +
  labs(
  y = "",
  title = "Chromosome 1 DNA/hAT "
  )


pdf("group_class_plots.pdf", width = 8, height = 6)

groups  <- unique(combined$group)
classes <- unique(combined$repeat_classfam)

# If using repeat_class and want to remove some categories 
#remove_vals <- c("SINE?", "")
#x <- x[!x %in% remove_vals]

for (g in groups) {
  for (c in classes) {
    group <- g
    class <- c
    
    df <- combined %>%
      filter(group == g, repeat_classfam == c)
    
    if (nrow(df) == 0) next
    
    p <- ggplot(df, aes(x = repeat_center/1e6, y = species, fill = species)) +
      geom_density_ridges(
        scale = 5,
        alpha = 0.7,
        linewidth = 0.8,                     # ridge outline thickness
        
      ) +
      scale_fill_manual(
        values = c("#AFBC88", "#618B4A"),
        labels = c(expression(italic("G. firmus")), 
                   expression(italic("G. pennsylvanicus")))
      ) +
      scale_x_continuous(
        expand = c(0, 0),
        name = "Genomic position (Mb)",
        breaks = seq(0, max(test$repeat_center)/1e6, by = 50),  # 0.03 Mb = 30 kb
        labels = scales::number_format(accuracy = 0.1) 
      ) +
      theme_ridges(
        grid = FALSE,
        center_axis_labels = TRUE
      ) +
      theme(
        legend.position = "none"
      ) +
      labs(
        y = "",
        title = paste(group, class)
      )
    
    print(p)   # ← each print() = new PDF page
  }
}

dev.off()



chr_lengths <- combined %>%
  group_by(species, group) %>%
  summarise(
    chr_length = max(query_end, na.rm = TRUE),
    .groups = "drop"
  )


combined <- combined %>%
  left_join(chr_lengths, by = c("species", "group")) %>%
  mutate(
    repeat_center = (query_start + query_end) / 2,
    rel_pos = repeat_center / chr_length
  )
#rel_pos = 0   #chromosome start
#rel_pos = 0.5 #midpoint
#rel_pos = 1   #chromosome end

ggplot(test, aes(x = rel_pos, y = species)) + 
  geom_density_ridges(
    scale = 5,
    alpha = 0.7,
    size = 0.8,
    aes(fill = species)
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    labels = scales::number
  ) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)




