library(tidyverse)

# Load the TSV file
df <- read_tsv("softmask_stats.tsv")

df <- df %>%
  # Extract chromosome number
  mutate(chr_num = as.numeric(str_extract(Region, "\\d+"))) %>%
  
  # Shift chromosome numbering
  mutate(new_chr = chr_num - 1) %>%
  
  # Build new Region column
  mutate(
    Region = case_when(
      chr_num == 1 ~ "X Chromosome",                         # special rename
      chr_num > 1  ~ paste0("Chromosome ", new_chr),         # shift numbers
      Region == "Unplaced Scaffolds" ~ "Unplaced Scaffolds", # keep
      TRUE ~ Region
    )
  ) %>%
  
  # Remove helper columns
  select(-chr_num, -new_chr)

# Set the levels for plotting 
df$Region <- factor(
  df$Region,
  levels = c(
    "X Chromosome",
    paste0("Chromosome ", 1:14),
    "Unplaced Scaffolds"
  ),
  ordered = TRUE
)


# Species name mapping
mapping <- c(
  "Gpenn.chr.final.masked.numt.sorted.renamed" = "G.pennsylvanicus",
  "Gfirm.chr.final.masked.numt.sorted.renamed" = "G.firmus",
  "GCA_046254815.1_ASM4625481v1_genomic" = "G.assimilis",
  "GCA_965638035.1_iqGryBima1.hap1.1_genomic.fna" = "G.bimaculatus"
)

df <- df %>%
  mutate(Species = recode(Species, !!!mapping))

# Custom color palette
colors <- c(
  "G.pennsylvanicus" = "#618B4A",
  "G.firmus" = "#AFBC88",
  "G.assimilis" = "#7AA095",
  "G.bimaculatus" = "#49306B"
)

# ---- Plot ----
custom_theme <- theme_classic(base_size = 14) +  # Use larger font sizes for readability
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Bold and centered title
    axis.title = element_text(size = 24, face = "bold"),  # Bold axis labels
    axis.text = element_text(size = 24, face = "bold"),  # Adjust axis tick labels size
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.margin = margin(10, 10, 10, 10),  # Adjust plot margins for better spacing
  )

p <- ggplot(df, aes(x = Region, y = Softmask_Percent, fill = Species)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 24, face = "bold"),          # tick labels
    axis.title = element_text(size = 24, face = "bold"),         # axis labels
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),                          # remove minor gridlines
    plot.margin = margin(10, 10, 10, 50),
    axis.ticks.length = unit(5, "pt"),
    axis.ticks = element_line(size = 1.5)# margins
  ) +
  labs(
    x = "",
    y = "Softmask (%)",
    title = "Softmask Percentage by Chromosome Across Species"
  ) 

# Save PNG
ggsave("softmask_boxplot.png", p, width = 14, height = 6, dpi = 300)

# Display plot
p
