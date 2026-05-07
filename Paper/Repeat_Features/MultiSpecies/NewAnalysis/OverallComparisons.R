library(tidyverse)
source("functions.R")

Gbimac_tbl <- parse_repeat_summary("Gbim_v2.2.whole_summary.tbl", "G.bimaculatus") %>% distinct() %>% drop_na() %>% 
  group_by(class, species) %>% summarise(total_bp = sum(BpMasked_per_class_family), 
                                total_percent = sum(PercentMasked_per_class_family))

Gasim_tbl <- parse_repeat_summary("Gassim_whole.summary.tbl", "G.assimilis") %>% distinct() %>% drop_na() %>% 
  group_by(class, species) %>% summarise(total_bp = sum(BpMasked_per_class_family), 
                                total_percent = sum(PercentMasked_per_class_family))


Gpenn_tbl <- parse_repeat_summary("Gpenn.chr.final.fasta.summary.tbl", "G.pennsylvanicus")%>% distinct() %>% drop_na() %>% 
  group_by(class, species) %>% summarise(total_bp = sum(BpMasked_per_class_family), 
                                total_percent = sum(PercentMasked_per_class_family))

Gfirm_tbl <- parse_repeat_summary("Gfirm.chr.final.fasta.summary.tbl", "G.firmus")  %>% distinct() %>% drop_na() %>% 
  group_by(class, species) %>% summarise(total_bp = sum(BpMasked_per_class_family), 
                                total_percent = sum(PercentMasked_per_class_family))

species_binded <- rbind(Gasim_tbl, Gbimac_tbl, Gpenn_tbl, Gfirm_tbl)

# Remove class_family that aren't in all 4
unique_class <- species_binded %>%
  group_by(class) %>%
  filter(n_distinct(species) == 4) %>%
  ungroup()

## Make a table
# Make a table by BP
latex_cells <- unique_class %>% 
  select(species, class, total_bp) %>%
  mutate(total_bp = format(total_bp, big.mark = ",", scientific = FALSE)) %>%  
  pivot_wider(
    names_from = species,
    values_from = total_bp,
    values_fill = "0.00"
  ) %>%
  arrange(class)

latex_rows <- apply(latex_cells, 1, function(x) {
  paste(x, collapse = " & ")
})

writeLines(latex_rows, "RE_Class_cells_BP.txt")


# if you want a table of percentages
latex_cells <- unique_class %>% 
  select(species, class, total_percent) %>%
  pivot_wider(
    names_from = species,
    values_from = total_percent
  ) %>%
  arrange(class)

latex_rows <- apply(latex_cells, 1, function(x) {
  paste(x, collapse = " & ")
})

writeLines(latex_rows, "RE_Class_cells_Percentages.txt")



## If we want to compare these better should we do BP repeats per say MB of Genome. 