library(tidyverse)
library(ggridges)

DM_f <- read_tsv("firmus.DMs.txt", col_names = FALSE)
colnames(DM_f) <- c("seqid", "location")

DM_g <- read_tsv("pennsylvanicus.DMs.txt", col_names = FALSE)
colnames(DM_g) <- c("seqid", "location")

common <- intersect(unique(DM_f$seqid), unique(DM_g$seqid))

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
    group = factor(group, levels = group_levels)
  )

colors <- c("Gpenn" = "#618B4A", "Gfirm" = "#AFBC88")

# Frequncy of Fixed Allele Along Chromsomes
ggplot(DM, aes(x = location, fill = species)) + 
  geom_histogram(binwidth = 1e6) + facet_wrap(~ seqid) +
  labs(x = "1 Mb Bins", 
       y = "Frequency")


# Location of Fixed Alleles along chromsomes 
ggplot(DM, aes(x = location/1e6, y = species, color = species)) +
  geom_point(size = 1.5) +
  facet_wrap(~ group, scales = "free_x") +
  labs(x = "Genomic position (Mb)", y = "Species", 
       title = "Fix Alleles Genomic Location Along Chromsomes") +
  scale_color_manual(values = colors)


# Add the chromosome lengths so we know the locations context
make_chr_lengths <- function(prefix,
                             chr_lengths_tsv,
                             asm_key_tsv,
                             group_levels) {
  
  read_tsv(chr_lengths_tsv, col_names = FALSE) %>%
    rename(Renamed = X1, Length_BP = X2) %>%
    left_join(
      read_tsv(asm_key_tsv),
      by = "Renamed"
    ) %>%
    select(Orginal, Length_BP) %>%
    slice_head(n = 15) %>%
    rename(seqid = Orginal) %>%
    mutate(
      seqid = case_when(
        seqid == "chr_1" ~ "X_chr",
        seqid %in% paste0("chr_", 2:15) ~
          paste0("chr_", as.numeric(str_remove(seqid, "chr_")) - 1),
        TRUE ~ "Unplaced"
      ),
      seqid = factor(seqid, levels = group_levels),
      species = prefix
    ) %>%
    arrange(seqid)
}

Gpenn_lengths <- make_chr_lengths(
  prefix = "Gpenn",
  chr_lengths_tsv = "Gpenn.chr.lengths.tsv",
  asm_key_tsv = "Gpenn.asm.key.tsv",
  group_levels = group_levels
)

Gfirm_lengths <- make_chr_lengths(
  prefix = "Gfirm",
  chr_lengths_tsv = "Gfirm.chr.lengths.tsv",
  asm_key_tsv = "Gfirm.asm.key.tsv",
  group_levels = group_levels
)
# Df with chromosome length information
chr_lengths <- bind_rows(Gpenn_lengths, Gfirm_lengths)



ggplot(DM, aes(x = location / 1e6, y = species, color = species)) +
  geom_blank(
    data = chr_lengths,
    aes(x = 0)
  ) +
  geom_blank(
    data = chr_lengths,
    aes(x = Length_BP / 1e6)
  ) +
  geom_point(size = 1) +
  facet_wrap(~ group, scales = "free_x") +
  labs(x = "Genomic position (Mb)", y = "Species")



# density style 
ggplot(
  DM,
  aes(
    x = location / 1e6,
    y = species,
    fill = species
  )
) +
  geom_density_ridges(
    scale = 1,
    rel_min_height = 0.01,
    alpha = 0.7,
    color = NA
  ) +
  facet_wrap(~ group, scales = "free_x") +
  labs(
    x = "Genomic position (Mb)",
    y = "Species"
  )


# Histogram Style
ggplot(
  DM,
  aes(
    x = location / 1e6,
    y = species,
    fill = species
  )
) +
  geom_blank(data = chr_lengths, aes(x = 0)) +
  geom_blank(data = chr_lengths, aes(x = Length_BP / 1e6)) +
  
  geom_density_ridges(
    stat = "binline",
    binwidth = 1,
    scale = 1,
    alpha = 0.8,
    color = NA
  ) +
  
  facet_wrap(~ group, scales = "free_x") +
  labs(
    x = "Genomic position (Mb)",
    y = "Species"
  )


############################################################################################################
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


# Allele density normalized per Mb
allele_density <- DM %>%
  mutate(bin = floor(location / window_size) * window_size) %>%
  group_by(species, seqid, bin) %>%
  summarise(
    density_per_mb = n() / (window_size / 1e6),
    .groups = "drop"
  ) %>%
  mutate(pos_mb = bin / 1e6)

# Feature density summed across all class_family
feature_density <- rm_df %>%
  filter(group %in% chr) %>%
  mutate(bin = floor(midpoint / window_size) * window_size) %>%
  group_by(species, group, bin) %>%
  summarise(
    density_per_mb = n() / (window_size / 1e6),
    .groups = "drop"
  ) %>%
  mutate(pos_mb = bin / 1e6)


# -----------------------
# 1️⃣ Allele density per Mb
# -----------------------
allele_density <- DM %>%
  mutate(bin = floor(location / window_size) * window_size) %>%
  group_by(species, seqid, bin) %>%
  summarise(
    density_per_mb = n() / (window_size / 1e6),
    .groups = "drop"
  ) %>%
  mutate(pos_mb = bin / 1e6)

# -----------------------
# 2️⃣ Feature density per Mb
# -----------------------
feature_density <- rm_df %>%
  filter(group %in% chr) %>%
  mutate(bin = floor(midpoint / window_size) * window_size) %>%
  group_by(species, group, bin) %>%
  summarise(
    density_per_mb = n() / (window_size / 1e6),
    .groups = "drop"
  ) %>%
  mutate(pos_mb = bin / 1e6,
         seqid = group)  # rename for consistency

# -----------------------
# 3️⃣ Overlay plot
# -----------------------
ggplot() +
  # Fixed allele density
  geom_line(
    data = allele_density,
    aes(x = pos_mb, y = density_per_mb, color = "Fixed alleles"),
    size = 0.8
  ) +
  # Feature density
  geom_line(
    data = feature_density,
    aes(x = pos_mb, y = density_per_mb, color = "Features"),
    size = 0.8
  ) +
  facet_wrap(~ species + group, scales = "free_x") +
  scale_color_manual(values = c("Fixed alleles" = "red", "Features" = "blue")) +
  labs(
    x = "Genomic position (Mb)",
    y = "Density per Mb",
    color = "Track"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold")
  )


allele_density <- allele_density %>%
  mutate(
    group = case_when(
      seqid == "chr_1" ~ "X_chr",
      seqid %in% paste0("chr_", 2:15) ~ paste0("chr_", as.numeric(str_remove(seqid, "chr_")) - 1),
      TRUE ~ "Unplaced"
    ),
    group = factor(group, levels = group_levels)
  )

pdf("allele_feature_density_by_chr.pdf", width = 10, height = 6)

for (chr in unique(allele_density$group)) {
  
  ad <- allele_density %>% filter(group == chr)
  fd <- feature_density %>% filter(group == chr)
  
  if (nrow(ad) == 0 || nrow(fd) == 0) next
  
  p <- ggplot() +
    
    geom_line(
      data = ad,
      aes(x = pos_mb, y = density_per_mb, color = "Fixed alleles"),
      size = 0.8
    ) +
    
    geom_line(
      data = fd,
      aes(x = pos_mb, y = density_per_mb, color = "Features"),
      size = 0.8
    ) +
    
    facet_wrap(~ species, scales = "free_x") +
    
    scale_color_manual(values = c(
      "Fixed alleles" = "red",
      "Features" = "blue"
    )) +
    
    labs(
      title = paste("Chromosome", chr),
      x = "Genomic position (Mb)",
      y = "Density per Mb",
      color = "Track"
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      strip.text = element_text(face = "bold")
    )
  
  print(p)
}

dev.off()



