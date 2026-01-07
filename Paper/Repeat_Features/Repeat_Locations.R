

repeat_density <- rm_df %>%
  filter(group %in% chr) %>%
  mutate(
    window = floor(midpoint / window_size)
  ) %>%
  group_by(species, group, window) %>%
  summarize(
    pos_mb = mean(midpoint) / 1e6,
    repeat_count = n(),
    repeat_bp = sum(end - begin + 1),
    .groups = "drop"
  )

# line plot 
ggplot(repeat_density,
       aes(pos_mb, repeat_bp / 1e6, color = species)) +
  geom_line() +
  facet_wrap(~ group, scales = "free_x") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Genomic position (Mb)",
    y = "Repeat bp per Mb",
    title = "Spatial distribution of repeats along the genome"
  ) + 
  scale_color_manual(values = colors)


# This is for all repeat features
#can't integrate to 1 becuase over such large distances it makes it infinitely small
ggplot(
  rm_df,
  aes(
    x = midpoint_mb,        
    y = group,      
    fill = species
  )
) +
  geom_density_ridges(
    scale = 1,                 # Controls how much ridges overlap
    alpha = 0.5,
    rel_min_height = 0.01      # Ignored if < 1% of max height
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Genomic Position (Mb)",
    y = "Chromosome",
    title = "Chromosome-wide Repeat Midpoint distributions"
  ) + 
  scale_fill_manual(values = colors)








# This is for specific repeat class_family
# This is useful to look at. Keep Point **We can not compare heights because it is not intergrated to 1**. 
# We would need to compute the densitys manually and then plot them. This is a question I would like to
# discuss with Dr. Moore. Because they plots can lead to drastially different conclusions. 

# This is all qualitative. We would need to run our stats to see which distributions are statistically 
# different. 

# 1. Define the unique families to iterate over
values <- unique(rm_df$class_family)

# 2. Open the PDF device
# You can set the width and height in inches
pdf("Repeat_Distributions_by_Family.pdf", width = 10, height = 8)

# 3. Loop through each family
for (fam in values) {
  
  # Filter data for the specific family
  plot_data <- subset(rm_df, class_family == fam)
  
  # Create the plot
  p <- ggplot(
    plot_data, 
    aes(
      x = midpoint_mb,        
      y = group,      
      fill = species
    )
  ) +
    geom_density_ridges(
      scale = 1,                 
      alpha = 0.5,
      rel_min_height = 0.01      
    ) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Genomic Position (Mb)",
      y = "Chromosome",
      # Dynamically update the title for each page
      title = paste("Distribution for family:", fam)
    ) + 
    scale_fill_manual(values = colors)
  
  # 4. Explicitly print the plot (required inside loops)
  print(p)
}

# 5. Close the PDF device
dev.off()


# Repeated with classes 
rm_df <- rm_df %>% mutate(class = sub("/.*", "", class_family))

# 1. Define the unique classes to iterate over
values <- unique(rm_df$class)

# 2. Open the PDF device
# You can set the width and height in inches
pdf("Repeat_Distributions_by_Class.pdf", width = 10, height = 8)

# 3. Loop through each family
for (fam in values) {
  
  # Filter data for the specific family
  plot_data <- subset(rm_df, class_family == fam)
  
  # Create the plot
  p <- ggplot(
    plot_data, 
    aes(
      x = midpoint_mb,        
      y = group,      
      fill = species
    )
  ) +
    geom_density_ridges(
      scale = 1,                 
      alpha = 0.5,
      rel_min_height = 0.01      
    ) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Genomic Position (Mb)",
      y = "Chromosome",
      # Dynamically update the title for each page
      title = paste("Distribution for family:", fam)
    ) + 
    scale_fill_manual(values = colors)
  
  # 4. Explicitly print the plot (required inside loops)
  print(p)
}

# 5. Close the PDF device
dev.off()



