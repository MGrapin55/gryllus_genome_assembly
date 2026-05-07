# Getting our vector of interest from other raw df's 
p_ClassFamily <- p %>% filter(class == "DNA" | class == "MITE" | class == "Penelope") %>% 
  distinct(class_family) %>% pull()

f_ClassFamily <- f %>% filter(class == "LTR" | class == "RC" | class == "Satellite" | class == "SINE") %>%
  distinct(class_family) %>% pull()


# Load the functions that are are used in the loop
source("../functions.R")
window <- 5e6  # 5 Mb windows
group_levels <- c("X_chr", paste0("chr_", 1:14), "Unplaced")

# Species
species <- c("Gpenn", "Gfirm")

# List of RE Class to Test for enrichment
classes <- list(
  Gpenn = p_ClassFamily, 
  Gfirm = f_ClassFamily
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
    filter_class_family = NULL
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
  N <- nrow(rm_df)
  
  class_totals <- rm_df %>%
    count(class_family, name = "K")
  
  n_total <- sum(species_interest$n_repeats)
  
  region_totals <- species_interest %>%
    filter(class_family %in% classes[[sp]]) %>%
    group_by(class_family) %>%
    summarise(
      k = sum(n_repeats),
      .groups = "drop"
    ) %>%
    mutate(n = n_total)
  
  
  # -----------------------------
  # Hypergeometric test
  # -----------------------------
  results_sp <- region_totals %>%
    left_join(class_totals, by = "class_family") %>%
    filter(class_family %in% species_classes) %>%
    mutate(
      species = sp,
      N = N,
      p_value = phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    ) %>%
    select(species, class_family, N, K, n, k, p_value)
  
  results_list[[sp]] <- results_sp
}

final_results <- bind_rows(results_list)
final_results

final_res_signif <- final_results %>% filter(p_value < 0.05)

write_tsv(final_res_signif, file = "SpeciationRegions_Results/HypergeometricComparisons_ClassFamily_results.tsv",)


# Could also break this down by chromosome specific enriched regions (no just speciation regions overall)

# Rational: Most chromsomes have the same class but can vary widely on class/family 
# goal would be to retest the classes by chromosome 
# goal would be to retest class_family by chromosome