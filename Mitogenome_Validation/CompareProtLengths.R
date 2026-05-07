# Rscript to compare Protein Lengths from the 13 mitochrondrial genes


library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)



s <- c("pennsylvanicus", "firmus", "lineaticeps", "veletis") # exclude "bimac" (nameing convention and out group)

# Read and combine all TSVs
df <- map_dfr(s, function(i) {
  read_tsv(
    file = paste0("lengths/formated/", i, ".tsv"),
    col_names = c("gene", "length"),
    show_col_types = FALSE
  ) %>%
    mutate(
      species = i,
      gene = str_to_upper(gene)
    )
})

# Pivot to wide format (compare lengths across species)
df_wide <- df %>%
  pivot_wider(
    names_from = species,
    values_from = length
  )

df_wide

write_tsv(df_wide, file = "CDS_Mitogenomes_results.tsv")




# code for tRNA's

s <- c("Pennsylvanicus", "Firmus", "Lineaticeps", "Veletis", "Bimac") # exclude "bimac" (nameing convention and out group)

# Read and combine all TSVs
df <- map_dfr(s, function(i) {
  read_tsv(
    file = paste0("tRNA/", i, ".tRNA.tsv"),
    col_names = c("locus_tag", "product", "length"),
    show_col_types = FALSE, 
    skip = 1
  ) %>%
    mutate(
      species = i,
      product = str_to_upper(product)
    )
})

# Pivot to wide format (compare lengths across species)
df_wide <- df %>% select(-1) %>%
  pivot_wider(
    names_from = species,
    values_from = length
  )

df_wide

write_tsv(df_wide, file = "tRNA_Mitogenomes_results.tsv")