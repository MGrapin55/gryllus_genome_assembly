setwd("~/Downloads/MOORE_LAB_UNL/GRYLLUS_GENOME_ASSEMBLY/GIT_REPO/gryllus_genome_assembly/Paper/Repeat_Features/DivergenceEnrichment")

# Load some handy functions in
source(file = "../functions.R")

window <- 5e6  # 5 Mb windows
group_levels <- c("X_chr", paste0("chr_", 1:14), "Unplaced")
plot_fill <- c("#446739", "#e78e24")
species_levels <- c("Gpenn", "Gfirm")


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

lengths <- rbind(Gpenn_lengths, Gfirm_lengths)


## Analysis 
summary_by_class <- feature_density %>%
  left_join(lengths, by = c("species", "group")) %>%
  mutate(class = sub("/.*", "", class_family)) %>%
  
  # keep group here
  group_by(species, group, class) %>%
  summarise(
    repeat_bp = sum(n_repeats_bp),
    .groups = "drop"
  ) %>%
  
  # totals should still be per species
  group_by(species) %>%
  mutate(
    species_total = sum(repeat_bp),
    pct = 100 * repeat_bp / species_total
  ) %>%
  ungroup()


#Stats 
# Test for normal distribution 
ggplot(summary_by_class,
       aes(x = repeat_bp, fill = species)) +
  geom_histogram(bins = 100) +
  theme_classic() + facet_wrap(~class, scales = "free_y") +
  scale_fill_manual(values = plot_fill)
# I am just gonna use a non parametric test

pvals <- data.frame(class = character(),
                    p_value = numeric(),
                    stringsAsFactors = FALSE)

for (i in unique(summary_by_class$class)) {
  
  filter_df <- summary_by_class %>%
    filter(class == i)
  
  # run Mann-Whitney U test
  res <- wilcox.test(repeat_bp ~ species,
                     data = filter_df,
                     exact = TRUE,       # safer if ties/large n
                     alternative = "two.sided")
  
  # store result
  pvals <- rbind(pvals,
                 data.frame(class = i,
                            p_value = res$p.value,
                            test_statistic= res$statistic))
  
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

summary_by_class$class <- factor(summary_by_class$class,
                              levels = class_order)




# after computing pvals
pvals_pos <- summary_by_class %>%
  group_by(class) %>%
  summarise(y = max(repeat_bp, na.rm = TRUE) * 1.08) %>%
  left_join(pvals, by = "class") %>%
  mutate(label = ifelse(p_adj < 0.001, "***",
                        ifelse(p_adj < 0.01, "**",
                               ifelse(p_adj < 0.05, "*", " "))))


pvals_pos$class <- factor(pvals_pos$class,
                                 levels = class_order)

pvals_pos$class <- factor(pvals_pos$class,
                          levels = class_order)
## Rename and Relevel
summary_by_class <- summary_by_class %>%
  mutate(
    species_full = case_when(
      species == "Gpenn" ~ "G.pennsylvanicus",
      species == "Gfirm" ~ "G.firmus",
      TRUE ~ NA_character_
    ),
    species_full = factor(species_full,
                          levels = c("G.pennsylvanicus", "G.firmus"))
  )

P <- ggplot(summary_by_class,
            aes(x = species_full, y = repeat_bp, fill = species_full)) +
  geom_boxplot(width = 0.7, outlier.size = 0.7) +
  facet_wrap(~ class, scales = "free_y", axes = "all_x", ncol = 5) +
  scale_fill_manual(values = plot_fill) +
  
  labs(y = "Repeat Base Pairs", x = NULL) +
  
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
    strip.background = element_rect(color = "white", linewidth = 0))  +
  geom_text(data = pvals_pos,
            aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE,
            size = 6,
            fontface = "bold")

  P
ggsave(filename = "REClass_BPByChromosomeSpecies.png", plot = P, device = "png", path = "./Plot", 
       dpi = 300, width = 11.69 , height = 8.27, units = "in")
ggsave(filename = "REClass_BPByChromosomeSpecies.jpeg", plot = P, device = "png", path = "./Plot", 
       dpi = 300, width = 11.69 , height = 8.27, units = "in")
ggsave(filename = "REClass_BPByChromosomeSpecies.pdf", plot = P, device = "pdf", path = "./Plot", 
       dpi = 300, width = 11.69 , height = 8.27, units = "in")



###########################################################################################################

# LOGIC: 
# We tested for differences between each class of repeat between Gpenn and Gfirm. We found LTR and Satillites 
# are signifigantly more abundant in Gfirm. (Result 1)
# Barplot

# It would we interesting to see which chromosomes have the most repeats of this class 
# (i.e are they associated with speciation dynamics) but want to make sure chromosome size is not accounting
# for a large variation (R^2 result 2).

# Null: there is no correlation between chromsome size (BP) and Repeats (BP) (Biologically Interesting)
# Alternative: there is a correlation between chromsome size (BP) and Repeats (BP)
# A weak but significant correlation suggests:
#   
#   Repeat content may be influenced by:
#   
#   Chromosome size (weak effect)
# 
# PLUS stronger factors like:
#   
# • recombination rate
# • gene density
# • chromatin state
# • centromere size
# • transposable element activity
# • local selection pressure

# So firms has a bigger X chromosome but lets see if LTR/SATILLITEs are explain by this well.
firm_LTR <- left_join(summary_by_class, Gfirm_lengths, by = c("species", "group")) %>% drop_na() %>% filter(class == "LTR")
ggplot(firm_LTR, aes(x = Length_BP, y = repeat_bp)) + geom_point()
cor.test(firm_LTR$repeat_bp, firm_LTR$Length_BP, method = "spearman")

# Spearman's rank correlation rho
# 
# data:  firm_LTR$repeat_bp and firm_LTR$Length_BP
# S = 148, p-value = 0.00255
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7357143 

# R squared 
(0.7357143)^2 * 100
# Explains 54% of variation for LTR's

firm_Satellite <- left_join(summary_by_class, Gfirm_lengths, by = c("species", "group")) %>% drop_na() %>% filter(class == "Satellite")
ggplot(firm_Satellite, aes(x = Length_BP, y = repeat_bp)) + geom_point()
cor.test(firm_Satellite$repeat_bp, firm_Satellite$Length_BP, method = "spearman")

# Spearman's rank correlation rho
# 
# data:  firm_Satellite$repeat_bp and firm_Satellite$Length_BP
# S = 222, p-value = 0.01954
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.6035714 
(0.6035714)^2 * 100
# Explains 36% of variation for Satellite



## How much does is chromsome size explained by repeat content?
summary_by_chr <- feature_density %>%
  left_join(lengths, by = c("species", "group")) %>%
  mutate(class = sub("/.*", "", class_family)) %>%
  
  # keep group here
  group_by(species, group) %>%
  summarise(
    repeat_bp = sum(n_repeats_bp),
    .groups = "drop"
  ) %>%
  
  # totals should still be per species
  group_by(species) %>%
  mutate(
    species_total = sum(repeat_bp),
    pct = 100 * repeat_bp / species_total
  ) %>%
  ungroup()

firm_chr <- left_join(summary_by_chr, Gfirm_lengths, by = c("species", "group")) %>% drop_na()
ggplot(firm_chr, aes(x = Length_BP, y = repeat_bp)) + geom_point()
cor.test(firm_chr$repeat_bp, firm_chr$Length_BP, method = "spearman")
#0.8 strong
# ~70% 

penn_chr <- left_join(summary_by_chr, Gpenn_lengths, by = c("species", "group")) %>% drop_na()
ggplot(penn_chr, aes(x = Length_BP, y = repeat_bp)) + geom_point()
cor.test(penn_chr$repeat_bp, penn_chr$Length_BP, method = "spearman")
#0.8
# ~73% 

# So the conclusion is that repeats are explaining a good portion of the variation between Gpenn and Gfirm.
# And LTR's and Satellite repeats are more adbundant in Firmus and still have a strong correlation with chromosome 
# size but explain less of the variation. (I.e other factors are also playing a role). And Gfirm LTR is a outlier 
# on the X chromosome suggesting a potenial dynamic to speciation rearrangements. 

# For all classes interesting to see if outliers are on same chromosomes or differnt. (Insection Result 3)

o <- outlier_data %>% group_by(species,group) %>% filter(is.extreme == "TRUE") %>% summarise(Count_Outliers = n())
filter(o, species == "Gpenn") %>% nrow()
filter(o, species == "Gfirm") %>% nrow()

t <- as.data.frame(table(o$group))

# Nothing that interesting here We see DNA has more on Gpenn everything else is pretty similar. 
