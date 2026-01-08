################################################################################################################
##
##                                Analysis for location of repeat elements
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

# Set factor levels and get new columns
rm_df <- rm_df %>%
  mutate(species = factor(species, levels = c("Gpenn", "Gfirm")),
         midpoint = abs((begin + end) / 2),
         midpoint_mb = midpoint / 1e6)
############################################################
## Shared parameters
############################################################

window_size <- 1e6  # 1 Mb windows
chr <- c("X_chr", str_c("chr_", 1:14))


# Plot colors 
colors <- c("Gpenn" = "#618B4A", "Gfirm" = "#AFBC88")
#############################################################

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


###################################################
# Taking this plot and subseting it by class_family
# 1. Define the unique families
values <- unique(rm_df$class_family)

# 2. Open the PDF device
pdf("Repeat_SpatialLocations_by_Family.pdf", width = 12, height = 8)

# 3. Loop through each family
for (fam in values) {
  
  # Step A: Filter and calculate density FOR THIS SPECIFIC FAMILY
  plot_data <- rm_df %>%
    filter(class_family == fam, group %in% chr) %>%
    mutate(window = floor(midpoint / window_size)) %>%
    group_by(species, group, window) %>%
    summarize(
      pos_mb = mean(midpoint) / 1e6,
      repeat_bp_mb = sum(end - begin + 1) / 1e6, # Pre-calculate Y value
      .groups = "drop"
    )
  
  # Step B: Check if there is data to plot (prevents empty page errors)
  if (nrow(plot_data) == 0) next
  
  # Step C: Create the plot using the newly calculated plot_data
  p <- ggplot(plot_data, aes(x = pos_mb, y = repeat_bp_mb, color = species)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ group, scales = "free_x") +
    theme_minimal(base_size = 14) +
    labs(
      x = "Genomic position (Mb)",
      y = "Repeat bp per Mb",
      title = paste("Spatial Distribution:", fam),
      subtitle = "Calculated via sliding window"
    ) + 
    scale_color_manual(values = colors)
  
  # 4. Explicitly print
  print(p)
}

# 5. Close device
dev.off()

# This appears to be fairly informative. It tells us where the peaks of repeats are along each chromosome.
# We see that most of the repeats are similar between the two species but some are have different localization. 


# Future Goal: Look at which peaks are younger or older by there classification. I am curious if certain differences
# are between new repeat activity or whether they are ancestral. 




##########################################################################

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

########################################

# Per all repeat elements
ggplot(repeat_density,
       aes(species, repeat_bp / 1e6, fill = species)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~ group, scales = "free_y") +
  labs(y = "Repeat bp per Mb") +
  scale_fill_manual(values = colors)

# There is no differnece vizually in the overall repeat elements distribution. 

# Does species A have higher repeat density than species B on chr3?”

# Taking this plot and subseting it by class_family
# 1. Define the unique families
values <- unique(rm_df$class_family)

# 2. Open the PDF device
pdf("Repeat_DistributionsBoxlplotsBinned_by_Family.pdf", width = 12, height = 8)

# 3. Loop through each family
for (fam in values) {
  
  # Step A: Filter and calculate density FOR THIS SPECIFIC FAMILY
  plot_data <- rm_df %>%
    filter(class_family == fam, group %in% chr) %>%
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
    
  # Step B: Check if there is data to plot (prevents empty page errors)
  if (nrow(plot_data) == 0) next
  
  # Step C: Create the plot using the newly calculated plot_data
  p <- ggplot(plot_data,
              aes(species, repeat_bp / 1e6, fill = species)) +
    geom_boxplot(outlier.size = 0.4) +
    facet_wrap(~ group, scales = "free_y") +
    labs(y = "Repeat bp per Mb",
         title = paste("Repeat Distribution Per Chromosome (bp per Mb):", fam)) +
    scale_fill_manual(values = colors)
  
  # 4. Explicitly print
  print(p)
}

# 5. Close device
dev.off()


# We can run a t-test or wilcoxon (mann whitney U) tests depending on normality. BH correction because of 
# multiple testing. 
# Can make violin plots if we want to vizualize the data sligthly differently. I know Dr. Moore likes that. 
# 

###################################################################################################
# Looking at the correlation between position and quantity of repeats. 
# There is not many strong correlation values. Means repeats are pretty distributed along the chromosome
# and not specifically enriched toward the end of chromosomes. (Hypothesis)
repeat_density_corr <- repeat_density %>%
  group_by(species, group) %>%
  summarize(
    spearman_r = cor(pos_mb, repeat_bp, method = "spearman")
  )

# Plot the correlation 
ggplot(repeat_density,
       aes(pos_mb, repeat_bp, color = species)) +
  geom_point() +
  facet_wrap(~ group, scales = "free_x") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Genomic position (Mb)",
    y = "Repeat bp",
    title = "Relatonship of between repeats and chromosome positon"
  ) + 
  scale_color_manual(values = colors)

# could add R values, line, and or lm if thought chromsome size was to be included. 
# Could also subset by class_family to see if a specific family is trending toward the end of chromosomes. 



###############################################################################################

# This is specifically for the binned data
# Where the mass is concentrated
ggplot(repeat_density,
       aes(pos_mb, group = interaction(group, species),
           fill = species, weight = repeat_bp)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ group, scales = "free_x")


# 1. Define the unique families to iterate over
values <- unique(rm_df$class_family)

# 2. Open the PDF device
# You can set the width and height in inches
pdf("Repeat_DistributionsBinned_by_Family.pdf", width = 10, height = 8)

# 3. Loop through each family
for (fam in values) {
  
  # Filter data for the specific family
  plot_data <- subset(rm_df, class_family == fam)
  
  # Create the plot
  p <- ggplot(repeat_density,
              aes(pos_mb, group = interaction(group, species),
                  fill = species, weight = repeat_bp)) +
    geom_density(alpha = 0.6) +
    facet_wrap(~ group, scales = "free_x") +
    theme_minimal(base_size = 14) +
    labs(
      x = "Genomic Position (Mb)",
      y = "density",
      # Dynamically update the title for each page
      title = paste("Distribution for family:", fam)
    ) + 
    scale_fill_manual(values = colors)
  
  # 4. Explicitly print the plot (required inside loops)
  print(p)
}

# 5. Close the PDF device
dev.off()

# This is not useful because we see that the localization are so similar that there mass is also very similar. 


####################################################################################################
# The cv value quantifies how uneven repeat density is along a chromosome, independent of total 
# repeat content, by measuring window-to-window variability relative to the mean.

# CV is best visualized as a chromosome-level summary, using point/line plots or heatmaps—not spatial plots
# because it measures how uneven repeats are, not where they occur.

#------------------------------------------------------------
# 1. Compute CV (coefficient of variation)
#------------------------------------------------------------

cv_df <- repeat_density %>%
  group_by(species, group) %>%
  summarize(
    cv = sd(repeat_bp) / mean(repeat_bp),
    .groups = "drop"
  )

#------------------------------------------------------------
# 2. Point + line plot by chromosome (BEST default)
#------------------------------------------------------------

ggplot(cv_df,
       aes(group, cv, color = species)) +
  geom_point(
    size = 3,
    position = position_dodge(width = 0.4)
  ) +
  geom_line(
    aes(group = species),
    position = position_dodge(width = 0.4),
    linewidth = 1.2
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Chromosome",
    y = "Coefficient of variation (CV)",
    title = "Repeat density heterogeneity by chromosome"
  ) + 
  scale_color_manual(values = colors)

#------------------------------------------------------------
# 3. Heatmap of CV values
#------------------------------------------------------------

ggplot(cv_df,
       aes(group, species, fill = cv)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(
    option = "magma",
    name = "CV"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Chromosome",
    y = "Species",
    title = "Chromosome-level repeat heterogeneity"
  )

#------------------------------------------------------------
# 4. Boxplot of CV distributions (chromosomes as replicates)
#------------------------------------------------------------

ggplot(cv_df,
       aes(species, cv, fill = species)) +
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
    y = "Coefficient of variation (CV)",
    title = "Distribution of repeat heterogeneity across chromosomes"
  )


# Test for normality
ggplot(cv_df, aes(x = cv, fill = species)) + geom_histogram() +
  facet_wrap(~ species) +
  scale_fill_manual(values = colors)
# Can't tell from just looking at the histogram. 
# shaprio-wilks test for normality 
# Standard subsetting: [rows, columns]
shapiro.test(cv_df$cv[cv_df$species == "Gpenn"])  # Suggests normality
shapiro.test(cv_df$cv[cv_df$species == "Gfirm"])  # Suggests normality

# Run a t-test to compare distributions
t.test(cv ~ species, data = cv_df)

# No significant difference


#------------------------------------------------------------
# 5. Ranked (lollipop) plot of CV values
#------------------------------------------------------------

ggplot(cv_df,
       aes(reorder(group, cv), cv, color = species)) +
  geom_segment(
    aes(
      xend = group,
      y = 0,
      yend = cv
    ),
    position = position_dodge(width = 0.5)
  ) +
  geom_point(
    size = 3,
    position = position_dodge(width = 0.5)
  ) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    x = "Chromosome",
    y = "Coefficient of variation (CV)",
    title = "Ranked repeat heterogeneity by chromosome"
  )

#------------------------------------------------------------
# 6. CV vs chromosome length (OPTIONAL, if length available)
#------------------------------------------------------------
# Requires a dataframe with chromosome lengths:
# chrom_lengths: group, chrom_length_mb

# cv_df <- cv_df %>%
#   left_join(chrom_lengths, by = "group")
#
# ggplot(cv_df,
#        aes(chrom_length_mb, cv, color = species)) +
#   geom_point(size = 3) +
#   geom_smooth(
#     method = "lm",
#     se = FALSE
#   ) +
#   theme_minimal(base_size = 14) +
#   labs(
#     x = "Chromosome length (Mb)",
#     y = "Coefficient of variation (CV)",
#     title = "Repeat heterogeneity vs chromosome length"
#   )


####################################################################################

# Calculating CV over 5mb windowns along each chromosome. Then comparing the CV distributions. 
cv_df <- repeat_density %>%
  mutate(
    block_5mb = ceiling(window / 5)
  ) %>%
  group_by(species, group, block_5mb) %>%
  summarize(
    mean_bp = mean(repeat_bp, na.rm = TRUE),
    sd_bp   = sd(repeat_bp, na.rm = TRUE),
    CV      = sd_bp / mean_bp,
    .groups = "drop"
  ) %>%
  filter(
    is.finite(CV),
    mean_bp > 0
  )

repeat_density %>%
  count(species, group, block_5mb = ceiling(window / 5)) %>%
  filter(n != 5)

cv_df %>%
  count(species, group)

ggplot(cv_df,
       aes(species, CV, fill = species)) +
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
    y = "Coefficient of variation (CV)",
    title = "Distribution of repeat heterogeneity across chromosomes (5 mb)"
  ) + 
  facet_wrap(~ group) + scale_fill_manual(values = colors)

# Overall the repetitive landscape does not appear that different at this resolutions. 

###################
###################
# This now subsets with the 5mb resolutions to measure heterogeneity (i.e patchy or not) across each chromosome
# by each specifc repeat class_family 
# Lets try this subseting by each class_family feature
# 1. Define the unique families to iterate over
values <- unique(rm_df$class_family)

# 2. Open the PDF device
# You can set the width and height in inches
pdf("Repeat_CV_Distributions5mb_Family.pdf", width = 10, height = 8)

# 3. Loop through each family
for (fam in values) {
  
  # 1. Calculate density for ONLY this family
  # We need this to get the 'repeat_bp' per window for this family
  fam_density <- rm_df %>%
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
  
  # 2. Calculate CV based on the family-specific density
  cv_df <- fam_density %>%
    mutate(
      block_5mb = ceiling(window / 5)
    ) %>%
    group_by(species, group, block_5mb) %>%
    summarize(
      mean_bp = mean(repeat_bp, na.rm = TRUE),
      sd_bp   = sd(repeat_bp, na.rm = TRUE),
      CV      = sd_bp / mean_bp,
      .groups = "drop"
    ) %>%
    filter(
      is.finite(CV),
      mean_bp > 0
    )
  
  # 3. Skip if no data for this family (prevents ggplot errors)
  if (nrow(cv_df) == 0) next
  
  # 4. Create the plot using cv_df (NOT plot_data)
  p <- ggplot(cv_df, aes(x = species, y = CV, fill = species)) +
    geom_boxplot(
      width = 0.6,
      outlier.size = 0.8,
      alpha = 0.7
    ) +
    # Note: Jitter is better with color=species if you want to see points clearly
    geom_jitter(
      width = 0.1,
      size = 1.5,
      alpha = 0.4 
    ) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Species",
      y = "Coefficient of variation (CV)",
      title = paste("Repeat heterogeneity (5Mb blocks):", fam)
    ) + 
    facet_wrap(~ group) + 
    scale_fill_manual(values = colors)
  
  # 5. Print
  print(p)
}

# 5. Close the PDF device
dev.off()


# So just glancing over the pdf it looks like at 5mb resolution the heterogenity is pretty much the same. 
# It was worth trying but I probably won't put to much more effort into this unless I get specific feedback
# from Dr. Moore to investigate further. 
# 
# Understanding the Output
# This plot is measuring spatial clustering.
# 
# High CV: The repeat family is "clumpy"—it's very dense in some parts of the 5Mb block and absent in others.
# 
# Low CV: The repeat family is distributed evenly across that section of the chromosome.




############################################################
## 3️⃣ QUESTION:
## "Where along the chromosome do repeats accumulate
##  in each species?"
## (Element density per Mb)
############################################################
# Removed the unplaced scaffolds because they were messing with the scales. This is because the group 
# "unplaced" is a combination of many scaffolds. This analysis treats each as a one sequences so having 
# all the combined sequnces made it have many more elements per region. 


# N repeats per Mb where above it was N repeats per Bp
# I think this this is a little more intutitive to understand. 
df_elem_density <- rm_df %>% 
  filter(group %in% chr) %>%
  mutate(window = floor(midpoint / window_size)) %>%
  group_by(species, window, group) %>%
  summarize(
    pos_mb = mean(midpoint_mb),
    density_per_mb = n() / (window_size / 1e6),
    .groups = "drop"
  )

ggplot(df_elem_density, aes(pos_mb, density_per_mb, color = species)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ group, scales = "free_x") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Genomic position (Mb)",
    y = "Element density (per Mb)",
    title = "Repeat density trajectories by species"
  ) +
  scale_color_manual(values = colors)
