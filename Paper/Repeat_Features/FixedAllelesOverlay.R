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


# Fixed Alleles data
DM_f <- read_tsv("firmus.DMs.txt", col_names = FALSE)
colnames(DM_f) <- c("seqid", "location")

DM_g <- read_tsv("pennsylvanicus.DMs.txt", col_names = FALSE)
colnames(DM_g) <- c("seqid", "location")

DM_g$species <- "Gpenn" 
DM_f$species <- "Gfirm" 

DM <- bind_rows(DM_g, DM_f)

group_levels <- c("X_chr", str_c("chr_", 1:14), "Unplaced")

DM <- DM %>%
  mutate(
    group = case_when(
      seqid == "chr_1" ~ "X_chr",
      seqid %in% paste0("chr_", 2:15) ~ paste0("chr_", as.numeric(str_remove(seqid, "chr_")) - 1),
      TRUE ~ "Unplaced"
    ),
    group = factor(group, levels = group_levels))
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

# Getting Feature Density 
feature_density <- rm_df %>%
  filter(group %in% chr) %>%
  mutate(bin = floor(midpoint / window_size) * window_size) %>%
  group_by(species, group, bin,class_family) %>%
  summarise(
    density_per_mb = n() / (window_size / 1e6),
    .groups = "drop"
  ) %>%
  mutate(pos_mb = bin / 1e6)

# Allele density
allele_density <- DM %>%
  mutate(bin = floor(location / window_size) * window_size) %>%
  group_by(species, group, bin) %>%
  summarise(
    density_per_mb = n() / (window_size / 1e6),
    .groups = "drop"
  ) %>%
  mutate(pos_mb = bin / 1e6)


# Now create the difference 
shared_units <- feature_density %>%
  distinct(species, group, bin, class_family, pos_mb) %>%
  count(group, bin, class_family, pos_mb) %>%
  filter(n == 2) %>%          # present in both species
  select(group, bin, class_family, pos_mb)


df_diff <- feature_density %>%
  semi_join(shared_units, by = c("group", "bin", "class_family", "pos_mb")) %>%
  pivot_wider(
    id_cols    = c(group, bin, class_family, pos_mb),
    names_from = species,
    values_from = density_per_mb
  ) %>%
  mutate(density_diff = Gfirm - Gpenn)

df_diff <- df_diff %>% mutate(
  class = sub("/.*", "", class_family)
)



# Get it in a continious order
# Define chromosome order
chrom_order <- c("X_chr", paste0("chr_", 1:14))  

# Compute max window per chromosome
chr_lengths <- df_diff %>%
  group_by(group) %>%
  summarize(max_pos_mb = max(pos_mb), .groups = "drop") %>%
  mutate(group = factor(group, levels = chrom_order)) %>%
  arrange(group)

chr_lengths <- chr_lengths %>%
  mutate(
    chr_start = lag(cumsum(max_pos_mb + 1), default = 0)  # +1 because windows start at 0
  )

df_diff <- df_diff %>%
  left_join(chr_lengths %>% select(group, chr_start), by = "group") %>%
  mutate(
    genome_pos_mb = chr_start + pos_mb  # continuous genome-wide coordinate
  )

allele_density <- allele_density %>%
  left_join(chr_lengths %>% select(group, chr_start), by = "group") %>%
  mutate(
    genome_pos_mb = chr_start + pos_mb  # continuous genome-wide coordinate
  )

# Plot as multiple tracks
ggplot(df_diff, aes(x = genome_pos_mb, y = density_diff, color = class)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linewidth = 0.6) +
  scale_x_continuous(
    breaks = chr_lengths$chr_start + chr_lengths$max_pos_mb / 2,
    labels = chr_lengths$group
  ) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Genome (MB)", y = "Density difference", 
       title = "Net Density Variation", 
       subtitle = ("Net = Gfirm.elementDensity - Gpenn.elementDensity")) + facet_wrap(~ class) +
  # Fixed allele density
  geom_point(
    data = allele_density,
    aes(x = genome_pos_mb, y = density_per_mb, color = "Fixed alleles", fill = species),
    linewidth = 0.8
  ) 


# Plots every combination at the Class level
groups   <- unique(df_diff$group)
features <- unique(df_diff$class)

pdf("Fixed_Alleles_density_diff_by_group_and_class.pdf", width = 10, height = 4)

for (g in groups) {
  for (f in features) {
    
    df_sub <- df_diff %>%
      filter(group == g, class == f)
    
    df_allele <- allele_density %>%
      filter(group == g)
    
    if (nrow(df_sub) == 0) next   # skip empty combinations
    
    p <- ggplot(df_sub, aes(x = pos_mb, y = density_diff)) +
      geom_point() +
      geom_hline(yintercept = 0, color = "red", linewidth = 0.6)  +
      theme_bw() +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      labs(
        x = "Chromosome Spans (MB)",
        y = "Density difference",
        title = paste("Repeat Elements:", f, "along", g)
      ) + 
      # Fixed allele density
      geom_point(
        data = df_allele,
        aes(x = pos_mb, y = density_per_mb, color = "Fixed alleles")
      ) 
    
    print(p)   # one page per plot
  }
}
dev.off()


feature_density_bin <- feature_density %>%
  group_by(species, group, bin, pos_mb) %>%
  summarise(
    density_per_mb = sum(density_per_mb),
    .groups = "drop"
  )

big_df <- full_join(
  feature_density_bin,
  allele_density,
  by = c("species", "group", "bin", "pos_mb"),
  suffix = c("_feature", "_allele")
)



fd <- feature_density_bin
ad <- allele_density

key_fd <- paste(fd$species, fd$group, fd$bin)
key_ad <- paste(ad$species, ad$group, ad$bin)

keep <- key_fd %in% key_ad

fd_m <- fd[keep, ]
ad_m <- ad[match(key_fd[keep], key_ad), ]

nrow(fd_m) == nrow(ad_m)   # should be TRUE

stopifnot(
  all(paste(fd_m$species, fd_m$group, fd_m$bin) ==
        paste(ad_m$species, ad_m$group, ad_m$bin))
)

cor.test(fd_m$density_per_mb, ad_m$density_per_mb, method = "spearman")


ggplot() +
  geom_point(aes(x = fd_m$density_per_mb,
                 y = ad_m$density_per_mb),
             alpha = 0.6) +
  geom_smooth(aes(x = fd_m$density_per_mb,
                  y = ad_m$density_per_mb),
              method = "lm", se = FALSE) +
  theme_bw()

plot_df <- data.frame(
  species = fd_m$species,
  group   = fd_m$group,
  bin     = fd_m$bin,
  pos_mb  = fd_m$pos_mb,
  feature_density = fd_m$density_per_mb,
  allele_density  = ad_m$density_per_mb
)

ggplot(plot_df, aes(x = feature_density, y = allele_density)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ group + species, scales = "free") +
  theme_bw() +
  labs(x = "Feature density / Mb",
       y = "Allele density / Mb")



