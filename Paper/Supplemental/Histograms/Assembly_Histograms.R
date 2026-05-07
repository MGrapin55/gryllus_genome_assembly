library(readr)
library(ggplot2)
library(dplyr)
library(scales)

# Colors 
# Penn: #446739 Firm: #e78e24

# Making a histograms for Gpenn and Gfirm Assembly lengths

Gpenn <- read_tsv("Gpenn.final.tsv", col_names = c("seqid", "length.BP"))


A <- ggplot(Gpenn, aes(x = length.BP)) +
  geom_histogram(bins = 100, color = "black", fill = "#446739") +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = function(x) paste0(x / 1e6, " Mb")
  ) + 
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    #axis.text.x = element_blank(),
    panel.grid = element_blank(),
    
    # Border around each heatmap
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    
    # Horizontal legend on top
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    
    axis.title = element_text(size = 16, face = "bold", , color = "black"),
    axis.text  = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(face = "bold", size = 14, color = "black"),
    legend.title = element_text(face = "bold", size = 16, color = "black"),
    axis.ticks.y = element_line(color = "black"), 
    strip.text.x = element_text(size = 16, face = "bold")
  ) +
  labs(
    x = expression(Log[10]("Basepair Scaffold Lengths"))
  )



Gfirm <- read_tsv("Gfirm.final.tsv", col_names = c("seqid", "length.BP"))


B <- ggplot(Gfirm, aes(x = length.BP)) +
  geom_histogram(bins = 100, color = "black", fill = "#e78e24") +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = function(x) paste0(x / 1e6, " Mb")
  ) + 
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    #axis.text.x = element_blank(),
    panel.grid = element_blank(),
    
    # Border around each heatmap
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    
    # Horizontal legend on top
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    
    axis.title = element_text(size = 16, face = "bold", , color = "black"),
    axis.text  = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(face = "bold", size = 14, color = "black"),
    legend.title = element_text(face = "bold", size = 16, color = "black"),
    axis.ticks.y = element_line(color = "black"), 
    strip.text.x = element_text(size = 16, face = "bold")
  ) +
  labs(
    x = expression(Log[10]("Basepair Scaffold Lengths"))
  )

p <- A / B

ggsave(filename = "ScaffoldHistograms.pdf", plot = p, width = 11.69 , height = 8.27, units = "in", dpi = 300)

# Caption: 
# Distribution of scaffold lengths. Basepair Lengths were log10 transformed to visualize. 
# G.pennsylvanius on top and G. firmus on the bottom

# Could Plot the chromosome if needed
