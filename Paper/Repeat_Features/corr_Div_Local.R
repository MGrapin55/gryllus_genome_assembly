################################################################################################################
##
##                        Analysis for correlation between divergence and localization
##
################################################################################################################
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

rm_df <- rm_df %>%
  mutate(
    species = factor(species, levels = c("Gpenn", "Gfirm")),
    midpoint = abs((begin + end) / 2),
    midpoint_mb = midpoint / 1e6
  )

############################################################
## Setup and shared parameters
############################################################
window_size <- 1e6                 # 1 Mb windows
window_mb   <- window_size / 1e6
chr <- c("X_chr", str_c("chr_", 1:14))

colors <- c("Gpenn" = "#618B4A", "Gfirm" = "#AFBC88")

############################################################
## Compute young-repeat and total repeat densities per window
############################################################
# As repeat density increases continuously across windows, does young repeat density also tend to increase?
young_density <- rm_df %>%
  filter(group %in% chr) %>%
  mutate(
    window   = floor(midpoint / window_size),
    is_young = perc_div <= 5
  ) %>%
  group_by(species, group, window) %>%
  summarize(
    young_repeats_per_mb = sum(is_young) / window_mb,
    total_repeats_per_mb = n() / window_mb,
    young_old_ratio_per_mb = young_repeats_per_mb / total_repeats_per_mb,
    .groups = "drop"
  )

############################################################
## Visualization: young vs total repeat density
############################################################

ggplot(data = young_density, aes(x = young_repeats_per_mb, y = total_repeats_per_mb, color = species)) +
  geom_point()



ggplot(
  young_density,
  aes(
    x = young_repeats_per_mb,
    y = total_repeats_per_mb,
    color = species
  )
) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Young repeats (≤5% divergence) per Mb",
    y = "Total repeat density (repeats per Mb)",
    title = "Relationship between young repeat density and total repeat density"
  ) + facet_wrap(~ group)

############################################################
## Statistics: Spearman correlation by species
############################################################
# overall correleation 
res <-young_density %>%
  group_by(species) %>%
  summarize(
    spearman_r = cor(
      young_repeats_per_mb,
      total_repeats_per_mb,
      method = "spearman",
      use = "complete.obs"
    ),
    p_value = cor.test(
      young_repeats_per_mb,
      total_repeats_per_mb,
      method = "spearman"
    )$p.value,
    .groups = "drop"
  )


# Correlation by chromosome
res_chr <-young_density %>%
  group_by(species, group) %>%
  summarize(
    spearman_r = cor(
      young_repeats_per_mb,
      total_repeats_per_mb,
      method = "spearman",
      use = "complete.obs"
    ),
    p_value = cor.test(
      young_repeats_per_mb,
      total_repeats_per_mb,
      method = "spearman"
    )$p.value,
    .groups = "drop"
  )
# This just tells us that there is a signfigant relationship between repeat density and younger ages? 
# (i.e the more dense the repeats are the more of them are younger.)
# This is not really biologially informative because you would just expect this that as there are more repeats 
# in a window that more of then will be young. (I.e. drift is somewhat random process)


############################################
ggplot(young_density,
       aes(species, young_old_ratio_per_mb, fill = species)) +
  geom_boxplot(
    width = 0.6,
    outlier.size = 0.8,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.1,
    size = 2,
    alpha = 0.6
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "",
    y = "Ratio of Young to Total Repeat Elements Per Mb",
    title = ""
  ) + 
  facet_wrap(~ group) + scale_fill_manual(values = colors)
# This does not tell us anything. They are all basically the same. 



###########################################

values <- unique(rm_df$class_family)

# 2. Open the PDF device
# You can set the width and height in inches
pdf("Ratio_YoungToOld_Family_Distribution.pdf", width = 10, height = 8)

# 3. Loop through each family
for (fam in values) {
  
  # Filter data for the specific family
  plot_data <- subset(rm_df, class_family == fam)
  
  plot_data <- plot_data %>%
    filter(group %in% chr) %>%
    mutate(
      window   = floor(midpoint / window_size),
      is_young = perc_div <= 5
    ) %>%
    group_by(species, group, window) %>%
    summarize(
      young_repeats_per_mb = sum(is_young) / window_mb,
      total_repeats_per_mb = n() / window_mb,
      young_old_ratio_per_mb = young_repeats_per_mb / total_repeats_per_mb,
      .groups = "drop"
    )
  
  # Create the plot
  p <- ggplot(plot_data,
              aes(species, young_old_ratio_per_mb, fill = species)) +
    geom_boxplot(
      width = 0.6,
      outlier.size = 0.8,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.1,
      size = 2,
      alpha = 0.6
    ) +
    theme_minimal(base_size = 14) +
    labs(
      x = "",
      y = "Ratio of Young to Total Repeat Elements Per Mb",
      title = paste("Ratio Young to Old:", fam, "Distribution")
    ) + 
    facet_wrap(~ group) + scale_fill_manual(values = colors)
  
  # 4. Explicitly print the plot (required inside loops)
  print(p)
}

# 5. Close the PDF device
dev.off()


