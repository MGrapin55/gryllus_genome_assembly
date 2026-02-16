library(tidyverse)
# Load the functions that are are used in the loop
source("../functions.R")
window <- 5e6  # 5 Mb windows
group_levels <- c("X_chr", paste0("chr_", 1:14), "Unplaced")

# Species
species <- c("Gpenn", "Gfirm")

# List of RE Class to Test for enrichment
# classes <- list(
#   Gpenn = c("DNA", "MITE", "Penelope"), 
#   Gfirm = c("LTR", "RC", "Satellite", "SINE")
# )

classes <- list(
  Gpenn = c(), 
  Gfirm = c("LTR","Satellite")
)
# RepeatMasker output files
RM_Files <- list(
  Gpenn = "../Gpenn.clean.out",
  Gfirm = "../Gfirm.clean.out"
)

# Allele files
Allele_Files <- list(
  Gpenn = "../pennsylvanicus.DMs.txt",
  Gfirm = "../firmus.DMs.txt"
)

# Chromosome length files
CHR_Files <- list(
  Gpenn = "../Gpenn.chr.lengths.tsv",
  Gfirm = "../Gfirm.chr.lengths.tsv"
)

# Key files
KEY_Files <- list(
  Gpenn = "../Gpenn.asm.key.tsv",
  Gfirm = "../Gfirm.asm.key.tsv"
)



results_list <- list()

for (sp in species) {
  
  message("\nProcessing species: ", sp)
  
  rm_file     <- RM_Files[[sp]]
  allele_file <- Allele_Files[[sp]]
  chr_file    <- CHR_Files[[sp]]
  key_file    <- KEY_Files[[sp]]
  species_classes <- classes[[sp]]
  
  # -----------------------------
  # RM data
  # -----------------------------
  rm_df <- read_rm_out(rm_file) %>%
    mutate(
      species = sp,
      midpoint = abs((begin + end) / 2),
      length_bp = abs(end - begin),
      class = sub("/.*", "", class_family)
    )
  
  # -----------------------------
  # Chromosome lengths
  # -----------------------------
  chr_lengths <- make_chr_lengths(
    chr_file, key_file, sp,
    group_levels = group_levels
  )
  
  # -----------------------------
  # Feature density (ALL classes)
  # -----------------------------
  feature_density <- make_feature_density(
    chr_length_df = chr_lengths,
    rm_df = rm_df,
    window = window
  ) %>%
    mutate(class = sub("/.*", "", class_family))
  
  # -----------------------------
  # Speciation analysis (RUN ONCE)
  # -----------------------------
  analysis <- Speciation_Analysis(
    species = sp,
    window = window,
    allele_file = allele_file,
    rm_file = rm_file,
    key = key_file,
    chr_length = chr_file,
    cutoff = 0.90,
    filter_class = NULL
  )
  
  AD <- analysis$allele_density
  FD <- analysis$feature_density
  
  full_df <- left_join(
    AD, FD,
    by = c("species", "group", "bin", "pos")
  )
  
  regions <- full_df %>%
    filter(Region_Speciation == 1) %>%
    select(group, bin)
  
  species_interest <- feature_density %>%
    semi_join(regions, by = c("group", "bin"))
  
  # -----------------------------
  # Hypergeometric components
  # -----------------------------
  # class_totals <- rm_df %>%
  #   group_by(group, class) %>%
  #   summarise(K = n(), .groups = "drop")
  # 
  # N_totals <- rm_df %>%
  #   count(group, name = "N")
  # 
  # region_totals <- species_interest %>%
  #   filter(class %in% classes[[sp]]) %>%
  #   group_by(group, class) %>%
  #   summarise(
  #     k = sum(n_repeats),
  #     .groups = "drop"
  #   )
  # 
  # n_totals <- species_interest %>%
  #   group_by(group) %>%
  #   summarise(n = sum(n_repeats), .groups = "drop")
  
  class_totals <- rm_df %>%
    group_by(group, class) %>%
    summarise(K = sum(length_bp), .groups = "drop")
  
  # N_totals <- rm_df %>%
  #   summarise(repeat_bp = sum(length_bp))
  
  N_totals <- rm_df %>%
    group_by(species, group) %>%
    summarise(group_bp = sum(length_bp), .groups = "drop") %>%
    mutate(N = sum(group_bp))
  
  
  region_totals <- species_interest %>%
    filter(class %in% classes[[sp]]) %>%
    group_by(group, class) %>%
    summarise(
      k = sum(n_repeats_bp),
      .groups = "drop"
    )
  
  n_totals <- species_interest %>%
    group_by(group) %>%
    summarise(n = sum(n_repeats_bp), .groups = "drop")
  
  
  # -----------------------------
  # Hypergeometric test
  # -----------------------------
  results_sp <- region_totals %>%
    left_join(class_totals, by = c("group", "class")) %>%
    left_join(N_totals,      by = "group") %>%
    left_join(n_totals,      by = "group") %>%
    replace_na(list(K = 0)) %>%   # class absent in rm_df on that chr
    mutate(
      species = sp
    ) %>%
    filter(class %in% species_classes) %>%
    mutate(
      species = sp,
      p_value = phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    ) %>%
    select(species, group, class, N, K, n, k, p_value)
  
  
  results_list[[sp]] <- results_sp
}

final_results <- bind_rows(results_list)
final_results

final_results <- final_results %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    Significant = p_adj < 0.05
  )

final_res_signif <- final_results %>% filter(p_adj < 0.05)

write_tsv(final_res_signif, file = "SpeciationRegions_Results/HypergeometricComparisons_BP_Class_Chromosome_results.tsv",)


filter(final_res_signif, species == "Gpenn") %>% count(class)

filter(final_res_signif, species == "Gfirm") %>% count(class)
