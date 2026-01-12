################################################################################################################
##
##                                Analysis for Divergence of repeat elements
##
################################################################################################################
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
############################################################
## Shared parameters
############################################################

window_size <- 1e6  # 1 Mb windows
chr <- c("X_chr", str_c("chr_", 1:14))

# Plot colors 
colors <- c("Gpenn" = "#618B4A", "Gfirm" = "#AFBC88")

############################################################
## 1️⃣ QUESTION:
## "Do the two species show different divergence trajectories
##  along the genome?"
## (Overlayed sliding-window median divergence)
############################################################

df_div_traj <- rm_df %>%
  mutate(window = floor(midpoint / window_size)) %>%
  group_by(species, window) %>%
  summarize(
    pos_mb = mean(midpoint_mb),
    median_div = median(perc_div),
    .groups = "drop"
  )

ggplot(df_div_traj, aes(pos_mb, median_div, color = species)) +
  geom_line(linewidth = 1) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Genomic position (Mb)",
    y = "Median percent divergence",
    title = "Comparative divergence trajectories"
  ) +
  scale_color_manual(values = colors)

# Notes: This plot X axis is not accurate because the genomic position is not continuous (RepeatMasker has 
# each entry in the fasta have its own coordinators.)

############################################################
##  QUESTION:
## "Do the two species show different divergence trajectories
##  along each chromosome?"
## (Overlayed sliding-window median divergence)
############################################################
df_div_traj_all <- rm_df %>%
  filter(group %in% chr) %>%
  mutate(window = floor(midpoint / window_size)) %>%
  group_by(species, window, group) %>%
  summarize(
    pos_mb = mean(midpoint_mb),
    median_div = median(perc_div),
    .groups = "drop"
  )

ggplot(df_div_traj_all, aes(pos_mb, median_div, color = species)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ group, scales = "free_x") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Genomic position (Mb)",
    y = "Median percent divergence",
    title = "Comparative divergence trajectories across chromosomes"
  ) +
  scale_color_manual(values = colors)

# This plot shows the divergence trajectories along each chromosome

##==================================================================================##
# Now this will separate them by class_family 
values <- unique(rm_df$class_family)

# 2. Open the PDF device
pdf("Repeat_MedianDivergence_by_Family.pdf", width = 12, height = 8)

# 3. Loop through each family
for (fam in values) {
  
  # Step A: Filter and calculate density FOR THIS SPECIFIC FAMILY
  plot_data <- rm_df %>%
    filter(class_family == fam, group %in% chr) %>%
    mutate(window = floor(midpoint / window_size)) %>%
    group_by(species, window, group) %>%
    summarize(
      pos_mb = mean(midpoint_mb),
      median_div = median(perc_div),
      .groups = "drop"
    )
  
  # Step B: Check if there is data to plot (prevents empty page errors)
  if (nrow(plot_data) == 0) next
  
  # Step C: Create the plot using the newly calculated plot_data
  p <- ggplot(plot_data, aes(pos_mb, median_div, color = species)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ group, scales = "free_x") +
    theme_minimal(base_size = 14) +
    labs(
      x = "Genomic position (Mb)",
      y = "Median percent divergence",
      title = paste("Comparative divergence trajectories across chromosomes:", fam)
    ) +
    scale_color_manual(values = colors)
  
  # 4. Explicitly print
  print(p)
}

# 5. Close device
dev.off()

############################################################
## 4️⃣ QUESTION:
## "How structured is the divergence distribution locally
##  along the chromosome in each species?"
## (Peak divergence density per window)
############################################################

df_peak_density <- rm_df %>%
  filter(group %in% chr) %>%
  mutate(window = floor(midpoint / window_size)) %>%
  group_by(species, window, group) %>%
  summarize(
    pos_mb = mean(midpoint_mb),
    peak_div_density = {
      d <- density(perc_div, from = 0, to = 50)
      max(d$y)
    },
    .groups = "drop"
  )

ggplot(df_peak_density, aes(pos_mb, peak_div_density, color = species)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ group, scales = "free_x") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Genomic position (Mb)",
    y = "Peak divergence density",
    title = "Local divergence density trajectories"
  ) + 
  scale_color_manual(values = colors)

# I do not completely understand what we are plotting. I get that it is computing the density(KDE). 

# This code collapses all repeat divergences in a genomic window into a single number that measures 
# how tightly clustered their divergence values are, using the height of a smoothed density peak.

# This is not useful at this point in time. Yet....


############################################################
## 5️⃣ QUESTION:
## "How does the FULL divergence density shift along the
##  chromosome for each species?"
## (Ridgeline density across genomic windows)
############################################################
# Plots it along a sliding window
rm_df2 <- rm_df %>%
  group_by(group) %>%  # chromosome
  mutate(
    window_mb = cut(
      midpoint_mb,
      breaks = seq(0, max(midpoint_mb, na.rm = TRUE), by = 5),
      include.lowest = TRUE
    )
  ) %>%
  ungroup()

plot_ridge_chr <- function(chrom) {
  ggplot(
    filter(rm_df2, group == chrom),
    aes(x = perc_div, y = window_mb, fill = species)
  ) +
    geom_density_ridges(scale = 2, alpha = 0.7) +
    facet_wrap(~ species) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Percent divergence",
      y = "Genomic window (Mb)",
      title = paste("Divergence density along", chrom)
    )
}

chroms <- unique(rm_df$group)

plots <- map(chroms, plot_ridge_chr)

plots[[1]]   # inspect



# Plots it along the chromosome
rm_df_chr <- rm_df %>%
  filter(group %in% chr)

ggplot(
  rm_df_chr,
  aes(
    x = perc_div,
    y = group,
    fill = species
  )
) +
  geom_density_ridges(
    alpha = 0.5,
    scale = 1,
    rel_min_height = 0.01, 
    aes(color = species)
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Percent divergence",
    y = "Chromosome",
    title = "Chromosome-wide divergence distributions (species overlaid)"
  )

# Integrates to 1 to make heights comparable
ggplot(
  rm_df_chr,
  aes(
    x = perc_div,
    y = group,
    fill = species,
    height = after_stat(density)
  )
) +
  geom_density_ridges(
    stat = "density",
    scale = 1,
    alpha = 0.5,
    rel_min_height = 0.01
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Percent divergence",
    y = "Chromosome",
    title = "Chromosome-wide divergence distributions (species overlaid)"
  )

#####################################################


ggplot(rm_df, aes(midpoint_mb, perc_div)) +
  geom_bin2d(bins = 75) +
  scale_fill_viridis_c() +
  facet_wrap(~ interaction(species, group)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Genomic position (Mb)",
    y = "Percent divergence",
    fill = "Density",
    title = "Joint density of divergence and position"
  )
# Need to revise this plot to get meaningful vizualizations.But it could be useful. 


