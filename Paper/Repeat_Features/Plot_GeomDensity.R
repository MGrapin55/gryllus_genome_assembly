# Plot Densitys on the same plot

# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggridges)
#########################################
##          Preprocessing              ##
#########################################
# Data cleaned (Make * a new columb)
Gpenn_rm <- read.table("Gpenn.clean.out",
                       skip = 3,
                       header = FALSE,
                       sep = "",
                       stringsAsFactors = FALSE,
                       fill = TRUE,
                       quote = "")

Gfirm_rm <- read.table("Gfirm.clean.out",
                       skip = 3,
                       header = FALSE,
                       sep = "",
                       stringsAsFactors = FALSE,
                       fill = TRUE,
                       quote = "")

# Column Names
columns <- c(
  "SW_score", "perc_div", "perc_del", "perc_ins",
  "query_sequence", "query_start", "query_end", "query_left",
  "strand",
  "repeat_name", "repeat_classfam",
  "repeat_start", "repeat_end", "repeat_left",
  "ID", "star"
)

colnames(Gpenn_rm) <- columns
colnames(Gfirm_rm) <- columns


# Add a species column 
Gpenn_rm$species <- "Gpenn"
Gfirm_rm$species <- "Gfirm"

# Combine Rows
combined <- bind_rows(Gpenn_rm, Gfirm_rm)

# Want to calculated the center of each repeat feature
# feature_center = abs(query_start + query_end / 2) (its the arthimetic mean)
combined <- combined %>% mutate(
  repeat_center = abs(query_start + query_end / 2)
)
# Now we want to split up repeat_classfam
combined <- combined %>%
  mutate(
    repeat_class = sub("/.*", "", repeat_classfam)
  )

# Subset chromosomes
chromosomes <- str_c("chr_", 1:15)

# Define desired order: chr_X first, then chr_1:14, then Unplaced
group_levels <- c("chr_X", str_c("chr_", 1:14), "Unplaced")

combined <- combined %>%
  mutate(
    group = case_when(
      query_sequence == "chr_1" ~ "chr_X",
      query_sequence %in% paste0("chr_", 2:15) ~ paste0("chr_", as.numeric(str_remove(query_sequence, "chr_")) - 1),
      TRUE ~ "Unplaced"
    ),
    group = factor(group, levels = group_levels)
  )

# Open PDF
pdf("group_class_density.pdf", width = 8, height = 6)

groups <- unique(combined$group)
classes <- unique(combined$repeat_classfam)

for (g in groups) {
  for (c in classes) {
    
    df <- combined %>%
      filter(group == g, repeat_classfam == c)
    
    if (nrow(df) < 2) next  # skip empty or single-point datasets
    
    # Precompute densities for each species
    dens <- df %>%
      group_by(species) %>%
      group_map(~ {
        if(nrow(.x) < 2) return(NULL)
        d <- density(.x$repeat_center / 1e6, n = 512, adjust = 1)
        data.frame(species = .y$species, x = d$x, y = d$y)
      }) %>%
      bind_rows()
    
    # Plot
    p <- ggplot(dens, aes(x = x, y = y, fill = species)) +
      geom_ribbon(aes(ymin = 0, ymax = y), alpha = 0.4) +
      geom_line(aes(color = species), linewidth = 1) +
      scale_fill_manual(values = c("Gfirm" = "#AFBC88", "Gpenn" = "#618B4A")) +
      scale_color_manual(values = c("Gfirm" = "#AFBC88", "Gpenn" = "#618B4A")) +
      scale_x_continuous(
        name = "Genomic position (Mb)",
        breaks = seq(0, max(dens$x), by = 50),
        labels = comma
      ) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      labs(title = paste(g, ":", c))
    
    print(p)
  }
}

dev.off()
