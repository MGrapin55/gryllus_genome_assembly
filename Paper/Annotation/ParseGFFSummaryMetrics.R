# Gff Summary Stats 

library(dplyr)
library(readr)
library(tidyr)

# Bash cleaning for the formatting 
# grep -E "^[A-Za-z]" <features>.txt \
# | sed 's/[[:space:]]\{2,\}/\t/g' \
# > <features>.tsv

Gpenn <- read_tsv("Gpenn_stat_features.tsv", col_names = c("Feature", "Stat")) %>%
  mutate(Species = "Gpenn")
  
Gfirm <- read_tsv("Gfirm_stat_features.tsv", col_names = c("Feature", "Stat")) %>%
  mutate(Species = "Gfirm")

keep_rows  <-  c(
    "Number of gene",
    "Number of mrna",
    "Number of exon",
    "mean mrnas per gene",
    "mean exons per mrna",
    "mean introns in exons per mrna",
    "Total gene length (bp)",
    "Total mrna length (bp)",
    "Total cds length (bp)",
    "Total exon length (bp)",
    "mean gene length (bp)",
    "mean mrna length (bp)",
    "mean cds length (bp)",
    "mean exon length (bp)",
    "mean intron in exon length (bp)"
    )

summary <- rbind(Gpenn, Gfirm)

summary <- summary %>% filter(Feature %in% keep_rows)

latex_cells <- summary %>%
  filter(Feature %in% keep_rows) %>%
  mutate(
    Feature = factor(Feature, levels = keep_rows),
    Stat = as.numeric(Stat)
  ) %>%
  pivot_wider(
    names_from  = Species,
    values_from = Stat,
    values_fill = 0
  ) %>%
  arrange(Feature)


latex_rows <- apply(latex_cells, 1, function(x) {
  paste(paste(x, collapse = " & "), "\\\\")
})


writeLines(latex_rows, "Gff_summaryMetrics.txt")


