################################################################################################################
##
##                                        Softmasked Comparison 
##                                        Author: Michael Grapin 
################################################################################################################
library(tidyverse)
library(ggplot2)

data <- read_tsv("softmasked_metrics_summary.tsv")

# Filter out the Unplaced scaffolds
data <- data %>% filter(Region != "Unplaced")

# Relevel species 
data$Species <- factor(data$Species, levels = c("G.pennsylvanicus", "G.firmus", "G.assimilis", "G.bimaculatus"))

## SUMMARY STATS
sink("Softmask_SummaryStats.txt")
# summary 
cat("Whole Data Set Summary Stats:")
summary(data)

cat("Gpenn Summary Stats:")
data %>% filter(Species == "G.pennsylvanicus") %>% summary()

cat("Gfirm Summary Stats:")
data %>% filter(Species == "G.firmus") %>% summary()

cat("Gassimilis Summary Stats:")
data %>% filter(Species == "G.assimilis") %>% summary()

cat("Gbimac Summary Stats:")
data %>% filter(Species == "G.bimaculatus") %>% summary()
sink()

## STATS
# Make it wide format 
wide_data <- data %>% 
  filter(Region != "Unplaced") %>%
  select(Region, Softmask_Prop, Species) %>%
  pivot_wider(
    names_from  = Species,
    values_from = Softmask_Prop)

# Drop Region and convert to numeric matrix
mat <- wide_data %>%
  select(-Region)%>%
  as.matrix()

# Make sure it's numeric
mode(mat) <- "numeric"

# Stats test
chisq.test(mat)

# Should I use this =? 
#kruskal.test(Softmask_Prop ~ Species, data = data)


## PLOTTING

# Custom color palette
colors <- c(
  "G.pennsylvanicus" = "#618B4A",
  "G.firmus" = "#AFBC88",
  "G.assimilis" = "#7AA095",
  "G.bimaculatus" = "#49306B"
)

custom_theme <- theme_classic(base_size = 14) +  # Use larger font sizes for readability
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Bold and centered title
    axis.title = element_text(size = 14, face = "bold"),  # Bold axis labels
    axis.text = element_text(size = 14, face = "bold"),  # Adjust axis tick labels size
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.margin = margin(10, 10, 10, 10),  # Adjust plot margins for better spacing
    legend.position = "none"
  )

# Reorder Species, add new color pallette, add breaks in y axis, and labels
p <- ggplot(data, aes(x = Species, y = Softmask_Prop, fill = Species)) + geom_boxplot(linewidth = 1) + 
  custom_theme + 
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    breaks = seq(30, 55, by = 5),
    limits = c(30, 57),
    expand = c(0.1, 0))+
  labs(x = "Species", 
       y = "% Softmasked")


# Display plot
p

# Save PNG
ggsave("Percent_softmask_boxplot.png", p, width = 11.69, height = 8.27, dpi = 300)





# Test specific chromsomes 
chr_10 <- wide_data %>% filter(Region == "10") %>% select(-Region) %>% as.matrix()
mode(chr_10) <- "numeric"
chisq.test(chr_10)


chr_13 <- wide_data %>% filter(Region == "13") %>% select(-Region) %>% as.matrix()
mode(chr_13) <- "numeric"
chisq.test(chr_13)
