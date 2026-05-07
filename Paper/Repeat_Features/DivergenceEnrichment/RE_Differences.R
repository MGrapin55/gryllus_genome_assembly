################################################################################################################
##
##                              Gpenn vs Gfirm differences in Repeat Elements
##                                        Author: Michael Grapin 
################################################################################################################
# The Goal is to produce a figure that shows on a chromosome level Repeat Features by Class differ.


setwd("~/Downloads/MOORE_LAB_UNL/GRYLLUS_GENOME_ASSEMBLY/GIT_REPO/gryllus_genome_assembly/Paper/Repeat_Features/DivergenceEnrichment")

# Load some handy functions in
source(file = "../functions.R")

window <- 5e6  # 5 Mb windows
group_levels <- c("X_chr", paste0("chr_", 1:14), "Unplaced")
species_levels <- c("Gpenn", "Gfirm")
color_species <- c("#446739", "#e78e24")


######################
# Reading in Data
#####################

# Reading in the data with a function read_rm_out() from functions.R 
p <- read_rm_out("../Gpenn.clean.out")
f <- read_rm_out("../Gfirm.clean.out")

# Add a species column 
p$species <- "Gpenn"
f$species <- "Gfirm"


# Get the midpoint value
p <- p %>%
  mutate(midpoint = abs((begin + end) / 2))

f <- f %>%
  mutate(midpoint = abs((begin + end) / 2))


# Read in chromsome lengths
Gpenn_lengths <- make_chr_lengths("../Gpenn.chr.lengths.tsv", "../Gpenn.asm.key.tsv", "Gpenn", group_levels = group_levels)
Gfirm_lengths <- make_chr_lengths("../Gfirm.chr.lengths.tsv", "../Gfirm.asm.key.tsv", "Gfirm", group_levels = group_levels)

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

# Rbind here
feature_density <- rbind(Gpenn_feature_density, Gfirm_feature_density)

# Remove class_family values that aren't in both 
feature_density <- feature_density %>%
  group_by(class_family) %>%
  filter(n_distinct(species) == 2) %>%
  ungroup()




# Summary Class_family 
summary_class_family <- feature_density %>%
  mutate(class = sub("/.*", "", class_family)) %>%
  group_by(species, group, class_family) %>%
  summarise(n_repeats = sum(n_repeats), .groups = "drop") %>%
  group_by(species, group) %>%
  mutate(pct = 100 * n_repeats / sum(n_repeats)) %>%
  ungroup()


# bar plot 
# Summary by broad class
summary_class <- feature_density %>%
  mutate(class = sub("/.*", "", class_family)) %>%
  group_by(species, group, class) %>%
  summarise(n_repeats = sum(n_repeats), .groups = "drop") %>%
  group_by(species, group) %>%
  mutate(pct = 100 * n_repeats / sum(n_repeats)) %>%
  ungroup()

summary_class <- summary_class %>%
  mutate(species = factor(species, levels = c("Gpenn", "Gfirm")))

pdf("class_counts_chromosome_bar.pdf", width = 8, height = 6)

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
    scale_fill_manual(values = color_species) +
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










# Test for normal distribution 
ggplot(summary_class,
       aes(x = n_repeats, fill = species)) +
  geom_histogram(bins = 100) +
  theme_classic() + facet_wrap(~class, scales = "free_y") +
  scale_fill_manual(values = color_species)
# I am just gonna use a non parametric test

pvals <- data.frame(class = character(),
                    p_value = numeric(),
                    stringsAsFactors = FALSE)

for (i in unique(summary_class$class)) {
  
  filter_df <- summary_class %>%
    filter(class == i)
  
  # run Wilcoxon test
  res <- wilcox.test(n_repeats ~ species,
                     data = filter_df,
                     exact = TRUE,       # safer if ties/large n
                     alternative = "two.sided")
  
  # store result
  pvals <- rbind(pvals,
                 data.frame(class = i,
                            p_value = res$p.value))
  
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
  scale_fill_manual(values = color_species) +
  
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

ggsave(filename = "REClass_CountsByChromosomeSpecies.png", plot = P, device = "png", path = "./Plot", 
      dpi = 300, width = 11.69 , height = 8.27, units = "in")
