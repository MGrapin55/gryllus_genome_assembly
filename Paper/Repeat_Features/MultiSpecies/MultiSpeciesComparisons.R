################################################################################################################
##
##                              Multi Gryllus Species RepeatMasker Comparison
##                                        Author: Michael Grapin 
################################################################################################################
out_dir = "~/Downloads/MOORE_LAB_UNL/GRYLLUS_GENOME_ASSEMBLY/GIT_REPO/gryllus_genome_assembly/Paper/Repeat_Features/MultiSpecies"

dir.create(out_dir)
setwd(out_dir)

# Repeat Masker "Cleaned" Files
Gpenn_file = "../Gpenn.clean.out"
Gfirm_file = "../Gfirm.clean.out"
Gassim_file = "../RepeatMaasker_out_files-selected/Gass_v1.0.genome.fasta.clean.out"
Gbimac_file = "../RepeatMaasker_out_files-selected/Gbim_v2.2.genome.fasta.clean.out"



# Load some handy functions in
source(file = "../functions.R")

# Reading in the data with a function read_rm_out() from functions.R 
p <- read_rm_out(Gpenn_file)
f <- read_rm_out(Gfirm_file)



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
p$species <- "Gpenn"
f$species <- "Gfirm"
a$species <- "Gassimilis"
b$species <- "Gbimaculatus"


## Figure out chr_nameing
group_levels <- c("X_chr", str_c("chr_", 1:14), "Unplaced")

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




# Combine dataframes
rm_df <- bind_rows(p, f, a, b)

# Remove class_family that aren't in all 4
rm_df <- rm_df %>%
  group_by(class_family) %>%
  filter(n_distinct(species) == 4) %>%
  ungroup()

# Set factor levels and get new columns
rm_df <- rm_df %>%
  mutate(species = factor(species, levels = c("Gpenn", "Gfirm", "Gassimilis", "Gbimaculatus")),
         midpoint = abs((begin + end) / 2),
         midpoint_mb = midpoint / 1e6 , 
           class = sub("/.*", "", class_family))

rm_df$class <- as.character(rm_df$class)


p %>%
  mutate( class = sub("/.*", "", class_family)) %>% 
  distinct(class)
###########################################################################################################
# Summary Statistics By Class_Family
summary <- rm_df %>%
  filter(group != "Unplaced") %>%
  group_by(species, group, class) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(percentage = count / sum(count) * 100)

summary_specific <- rm_df %>%
  filter(group != "Unplaced") %>%
  group_by(species, group, class_family) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(percentage = count / sum(count) * 100)


wide_df <- summary %>%
  select(species, group, class, count) %>%
  pivot_wider(names_from = class, values_from = count, values_fill = 0) %>% 
  select(1:12) %>%
  group_by(species) 

tab <- wide_df %>%   # your wide table
  group_by(species) %>%
  summarise(across(DNA:RC, sum), .groups = "drop")
tab <- tab %>% select(-MITE, -Penelope)

chisq.test(as.matrix(tab[,-1]))



tab <- wide_df %>%
  select(species, group, DNA) %>%
  pivot_wider(names_from = group, values_from = DNA, values_fill = 0)

chisq.test(as.matrix(tab[,-1]))


summary_species <- rm_df %>%
  filter(group != "Unplaced", class != "Unknown") %>%
  group_by(species, class) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(percentage = count / sum(count) * 100)



# The accompanying stacked barplot
ggplot(summary, aes(x = group, y = percentage, fill = class)) +
  geom_col(color = "black", width = 0.8) +
  facet_wrap(~ species, scales = "free_x") +
  labs(x = "Chromosome group",
       y = "Percentage",
       fill = "Class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Too much going on
ggplot(summary_specific, aes(x = group, y = percentage, fill = class_family)) +
  geom_col(color = "black", width = 0.8) +
  facet_wrap(~ species, scales = "free_x") +
  labs(x = "Chromosome group",
       y = "Percentage",
       fill = "Class_family") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Just by the species level 
ggplot(summary_species, aes(x = species, y = percentage, fill = class)) +
  geom_col(color = "black", width = 0.8) +
  labs(x = "Species",
       y = "Percentage",
       fill = "Class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

summary_species <- summary_species %>%
  group_by(species) %>%
  arrange(species, class) %>%
  mutate(label_pos = cumsum(percentage) - percentage/2)

ggplot(summary_species, aes(x = species, y = percentage, fill = class)) +
  geom_col(color = "black", width = 0.8) +
  geom_text(aes(y = label_pos,
                label = paste0(round(percentage, 1), "%")),
            color = "black", size = 3) +
  labs(x = "Species",
       y = "Percentage",
       fill = "Class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Pie chart
ggplot(summary_species, aes(x = "", y = percentage, fill = class)) +
  geom_col(width = 1, color = "black") +
  coord_polar(theta = "y") +
  facet_wrap(~ species) +
  labs(fill = "Class") +
  theme_void() +
  theme(strip.text = element_text(size = 12, face = "bold"))

## Format for latex table 
latex_cells <- summary_species %>% 
  select(species, class, percentage) %>%
  mutate(percentage = sprintf("%.2f", percentage), 
         percentage = paste0(percentage, "\\%")) %>%  
  pivot_wider(
    names_from = species,
    values_from = percentage,
    values_fill = "0.00"
  ) %>%
  arrange(class)

latex_rows <- apply(latex_cells, 1, function(x) {
  paste(x, collapse = " & ")
})

writeLines(latex_rows, "RE_Class_cells.txt")



##########################################################################################################
# Shared parameters
window_size <- 1e6  # 1 Mb windows
chr <- c("X_chr", str_c("chr_", 1:14))

# lets look at feature density along the chromosomes
feature_density <- rm_df %>%
  filter(group %in% chr) %>%
  mutate(bin = floor(midpoint / window_size) * window_size) %>%
  group_by(species, group, bin, class) %>%
  summarise(
    density_per_mb = n() / (window_size / 1e6),
    .groups = "drop"
  ) %>%
  mutate(pos_mb = bin / 1e6)


# The plot 
ggplot(feature_density, aes(pos_mb, density_per_mb, color = species)) +
  geom_line() +
  theme_minimal(base_size = 14) +
  labs(
    x = "Genomic position (Mb)",
    y = "Element density (per Mb)",
    title = "Repeat element density along chromosome"
  ) + facet_wrap(~ species)
#########################################################################################

feature_divergence <- rm_df %>%
  filter(group %in% chr) %>%
  mutate(bin = floor(midpoint / window_size) * window_size) %>%
  group_by(species, group, bin, class) %>%
  summarise(
    Avg_Percent_Divergence = mean(perc_div),
    .groups = "drop"
  ) %>%
  mutate(pos_mb = bin / 1e6)


ggplot(feature_divergence, aes(pos_mb, Avg_Percent_Divergence, color = species)) +
  geom_line() +
  theme_minimal(base_size = 14) +
  labs(
    x = "Genomic position (Mb)",
    y = "Average Percent Divergence (per Mb)",
    title = "Repeat element Divergence"
  ) + facet_wrap(~group)

young_divergence <- feature_divergence %>% filter(Avg_Percent_Divergence <= 5)

ggplot(young_divergence, aes(pos_mb, Avg_Percent_Divergence, color = class)) +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    x = "Genomic position (Mb)",
    y = "Average Percent Divergence (per Mb)",
    title = "Young Repeat element Divergence"
  ) + facet_wrap(~ species)

ggplot(feature_divergence, aes(x = Avg_Percent_Divergence, fill = species)) + geom_histogram() + 
  facet_wrap(~ species)


library(ggplot2)
library(ggridges)

ggplot(feature_divergence, aes(x = Avg_Percent_Divergence, y = species, fill = species)) +
  geom_density_ridges(alpha = 0.7, scale = 1) +
  labs(
    x = "Average Percent Divergence",
    y = "Species"
  ) +
  theme_ridges() +
  theme(legend.position = "none")


