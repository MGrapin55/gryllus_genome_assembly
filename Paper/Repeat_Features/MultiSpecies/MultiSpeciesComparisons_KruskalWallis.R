# Load some handy functions in
source(file = "../functions.R")

window <- 5e6  # 5 Mb windows
group_levels <- c("X_chr", paste0("chr_", 1:14), "Unplaced")

plot_fill <- c("#618B4A", "#AFBC88", "#7AA095", "#49306B")





# Read in all RepeatMasker Files
######################
# Reading in Data
#####################

# Reading in the data with a function read_rm_out() from functions.R 
p <- read_rm_out("../Gpenn.clean.out") %>% mutate(species = "Gpenn")
f <- read_rm_out("../Gfirm.clean.out") %>% mutate(species = "Gfirm")


# Read in chromsome lengths
Gpenn_lengths <- make_chr_lengths("../Gpenn.chr.lengths.tsv", "../Gpenn.asm.key.tsv", "Gpenn", group_levels = group_levels)
Gfirm_lengths <- make_chr_lengths("../Gfirm.chr.lengths.tsv", "../Gfirm.asm.key.tsv", "Gfirm", group_levels = group_levels)


# For Gbimac and Gassimils
Gassim_file = "../RepeatMaasker_out_files-selected/Gass_v1.0.genome.fasta.clean.out"
Gbimac_file = "../RepeatMaasker_out_files-selected/Gbim_v2.2.genome.fasta.clean.out"

# assimilis and bimac specific handling 
col <- c("score", "perc_div", "perc_del", "perc_ins", 
         "seq_id", "begin", "end", "left", 
         "strand", "repeat_name", "class_family", 
         "r_begin", "r_end", "r_left", "id", "star")

a <- read.table(Gassim_file,
                skip = 3,
                header = FALSE,
                sep = "",
                stringsAsFactors = FALSE,
                fill = TRUE,
                quote = "")
b <- read.table(Gbimac_file,
                skip = 3,
                header = FALSE,
                sep = "",
                stringsAsFactors = FALSE,
                fill = TRUE,
                quote = "")

# Assign column names and remove missing rows
colnames(a) <- col
a <- a %>% filter(!is.na(seq_id))
colnames(b) <- col
b <- b %>% filter(!is.na(seq_id))


# Add a species column 
a$species <- "Gassimilis"
b$species <- "Gbimaculatus"

# Bimac 
b <- b %>%
  mutate(
    group = recode(
      seq_id,
      chrX = "X_chr",
      !!!setNames(paste0("chr_", 1:14), paste0("chr", 1:14)),
      .default = "Unplaced"
    ),
    group = factor(group, levels = group_levels)
  )

# assimlis 
a <- a %>%
  mutate(
    group = case_when(
      seq_id == "Super-Scaffold_1" ~ "X_chr",
      str_detect(seq_id, "^Super-Scaffold_\\d+$") ~ {
        n <- as.numeric(str_remove(seq_id, "Super-Scaffold_"))
        ifelse(n >= 1 & n <= 14,
               paste0("chr_", n),
               "Unplaced")
      },
      TRUE ~ "Unplaced"
    ),
    group = factor(group, levels = group_levels)
  )

# Read in the chromsome length infomation
Gassim_lengths <- read_tsv("chr_lengths/Gassim_chr.lengths.tsv", col_names = c("NCBI_ID", "Length_BP")) %>% mutate(
  group = factor(c("X_chr", paste0("chr_", 1:14))), 
  species = "Gassimilis"
)


Gbimac_lengths <- read_tsv("chr_lengths/Gbiac2.2_chr.lengths.tsv", col_names = c("NCBI_ID", "Length_BP")) %>% mutate(
  group = factor(c(paste0("chr_", 1:14), "X_chr"), levels = c("X_chr", paste0("chr_", 1:14))), 
  species = "Gbimaculatus"
)



#######################
# Data Preprocessing 
######################
Gpenn_feature_density <- make_feature_density(
  chr_length_df = Gpenn_lengths,
  rm_df = p,
  window = window
)

Gfirm_feature_density <- make_feature_density(
  chr_length_df = Gfirm_lengths,
  rm_df = f,
  window = window
)

Gassim_feature_density <- make_feature_density(
  chr_length_df = Gassim_lengths,
  rm_df = a,
  window = window
)

Gbimac_feature_density <- make_feature_density(
  chr_length_df = Gbimac_lengths,
  rm_df = b,
  window = window
)

# Rbind here
feature_density <- rbind(Gpenn_feature_density, Gfirm_feature_density, Gassim_feature_density, Gbimac_feature_density )

# Remove class values that aren't in both 
feature_density <- feature_density %>%
  mutate(class = sub("/.*", "", class_family)) %>%
  group_by(class) %>%
  filter(n_distinct(species) == 4) %>%
  ungroup()

#######################
# Analysis 
######################

# Summary by broad class
summary_class <- feature_density %>%
  mutate(class = sub("/.*", "", class_family)) %>%
  group_by(species, group, class) %>%
  summarise(n_repeats = sum(n_repeats), .groups = "drop") %>%
  group_by(species, group) %>%
  mutate(pct = 100 * n_repeats / sum(n_repeats)) %>%
  ungroup()

# will have to update from here
summary_class <- summary_class %>%
  mutate(species = factor(species, levels = c("Gpenn", "Gfirm", "Gassimilis", "Gbimaculatus")))


pdf("plots/Mutlispecies_class_counts_chromosome_bar.pdf", width = 8, height = 6)

broad_classes <- unique(summary_class$class)

for (c in broad_classes) {
  
  df <- summary_class %>%
    filter(class == c)
  
  if (nrow(df) == 0) next
  
  p <- ggplot(df, aes(x = group, y = n_repeats, fill = species)) +
    geom_col(
      position = position_dodge(width = 0.9),
      color = "black",
      linewidth = 0.4
    ) +
    scale_fill_manual(values = plot_fill) +
    theme_classic() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      y = "Number of Identified Repeats",
      x = "Chromosome",
      title = paste("Repeat Class:", c)
    )
  
  print(p)
}
dev.off()

pdf("plots/Mutlispecies_class_percent_chromosome_bar.pdf", width = 8, height = 6)

broad_classes <- unique(summary_class$class)

for (c in broad_classes) {
  
  df <- summary_class %>%
    filter(class == c)
  
  if (nrow(df) == 0) next
  
  p <- ggplot(df, aes(x = group, y = pct, fill = species)) +
    geom_col(
      position = position_dodge(width = 0.9),
      color = "black",
      linewidth = 0.4
    ) +
    scale_fill_manual(values = plot_fill) +
    theme_classic() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      y = "% of Repeat Feature on Per Chromosome",
      x = "Chromosome",
      title = paste("Repeat Class:", c)
    )
  
  print(p)
}
dev.off()

library(rstatix)

# Test for normal distribution 
ggplot(summary_class,
       aes(x = n_repeats, fill = species)) +
  geom_histogram(bins = 100) +
  theme_classic() + facet_wrap(~class, scales = "free_y") +
  scale_fill_manual(values = plot_fill)
# I am just gonna use a non parametric test

library(dplyr)
library(rstatix)

pvals <- data.frame(class = character(),
                    p_value = numeric(),
                    stringsAsFactors = FALSE)

for (i in unique(summary_class$class)) {
  
  filter_df <- summary_class %>%
    filter(class == i)
  
  # Kruskal-Wallis test
  res <- kruskal_test(
    filter_df,
    n_repeats ~ species
  )
  
  pvals <- rbind(
    pvals,
    data.frame(
      class = i,
      p_value = res$p
    )
  )
}

pvals <- pvals %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    Significant = p_adj < 0.05
  )

pvals


class_order <- c("DNA","LTR","MITE","LINE","RC",
                 "Satellite","Penelope","SINE","PLE","Unknown")
species_order <- c("Gpenn", "Gfirm")

summary_class$class <- factor(summary_class$class,
                              levels = class_order)

pvals_pos$class <- factor(pvals_pos$class,
                          levels = class_order)


# after computing pvals
pvals_pos <- summary_class %>%
  group_by(class) %>%
  summarise(y = max(n_repeats, na.rm = TRUE) * 1.08) %>%
  left_join(pvals, by = "class") %>%
  mutate(label = ifelse(p_adj < 0.001, "***",
                        ifelse(p_adj < 0.01, "**",
                               ifelse(p_adj < 0.05, "*", " "))))



P <- ggplot(summary_class,
            aes(x = species, y = n_repeats, fill = species)) +
  geom_boxplot(width = 0.7, outlier.size = 0.7) +
  facet_wrap(~ class, scales = "free_y", axes = "all_x", ncol = 5) +
  scale_fill_manual(values = plot_fill) +
  
  labs(y = "Number of Repeats Features", x = NULL) +
  
  theme_classic(base_size = 14) +
  theme(
    text = element_text(face = "bold"),
    strip.text = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.spacing = unit(1.2, "lines"),
    
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    legend.position = "top",
    legend.justification = "center",
    legend.title = element_blank(),
    
    # Add a thick white border to the strip to simulate spacing
    strip.background = element_rect(color = "white", linewidth = 0)) +
  geom_text(data = pvals_pos,
            aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE,
            size = 6,
            fontface = "bold")

P

ggsave(filename = "MutliSpeciesREClass_CountsByChromosomeSpecies.png", plot = P, device = "png", path = "./Plot", 
       dpi = 300, width = 11.69 , height = 8.27, units = "in")




ggplot(summary_class,
       aes(x = species, y = n_repeats, fill = species)) +
  geom_boxplot(width = 0.7, outlier.size = 0.7) +
  facet_wrap(~ class, scales = "free_y", axes = "all_x", ncol = 5) +
  scale_fill_manual(values = plot_fill) +
  
  labs(y = "Number of Repeats Features", x = NULL) +
  
  theme_classic(base_size = 14) +
  theme(
    text = element_text(face = "bold"),
    strip.text = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.spacing = unit(1.2, "lines"),
    
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    legend.position = "top",
    legend.justification = "center",
    legend.title = element_blank(),
    
    # Add a thick white border to the strip to simulate spacing
    strip.background = element_rect(color = "white", linewidth = 0)) 


#################
Gassim_lengths <- select(Gassim_lengths, -NCBI_ID)
Gbimac_lengths <- select(Gbimac_lengths, -NCBI_ID)

lengths <- rbind(Gpenn_lengths, Gfirm_lengths, Gassim_lengths, Gbimac_lengths)

summary_overall <- feature_density %>%
  left_join(feature_density,lengths, by = c("species", "group")) %>%
  mutate(class = sub("/.*", "", class_family)) %>%
  
  group_by(species, group) %>%
  summarise(n_repeats_bp = sum(n_repeats_bp), .groups = "drop") %>%
  
  group_by(species) %>%
  mutate(
    species_total = sum(n_repeats_bp),
    pct = 100 * n_repeats_bp / length_BP
  ) %>%
  ungroup()

summary_overall <- feature_density %>%
  left_join(lengths, by = c("species", "group")) %>%
  group_by(species, group) %>%
  summarise(
    chr_repeat_bp = sum(n_repeats_bp),
    chr_length_bp = unique(Length_BP),
    .groups = "drop"
  ) %>%
  group_by(species) %>%
  summarise(
    genome_repeat_bp = sum(chr_repeat_bp),
    genome_length_bp = sum(chr_length_bp),
    pct = 100 * genome_repeat_bp / genome_length_bp,
    .groups = "drop"
  )

summary_by_class <- feature_density %>%
  left_join(lengths, by = c("species", "group")) %>%
  mutate(class = sub("/.*", "", class_family)) %>%
  group_by(species, class) %>%
  summarise(
    repeat_bp = sum(n_repeats_bp),
    .groups = "drop"
  ) %>%
  group_by(species) %>%
  mutate(
    species_total = sum(repeat_bp),
    pct = 100 * repeat_bp / species_total
  )



ggplot(summary_overall,
       aes(x = species, y = pct, fill = species)) +
  geom_boxplot(width = 0.7, outlier.size = 0.7) +
  scale_fill_manual(values = plot_fill) +
  
  labs(y = "Number of Repeats Features", x = NULL) +
  
  theme_classic(base_size = 14) +
  theme(
    text = element_text(face = "bold"),
    strip.text = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.spacing = unit(1.2, "lines"),
    
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    legend.position = "top",
    legend.justification = "center",
    legend.title = element_blank(),
    
    # Add a thick white border to the strip to simulate spacing
    strip.background = element_rect(color = "white", linewidth = 0)) 
