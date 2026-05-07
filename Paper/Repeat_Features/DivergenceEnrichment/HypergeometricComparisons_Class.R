# Running Hypergeometric Distribution Comparisons

# Outline of Hypergeometric distrubution 
# N = n total repeats (across all classes)
# K = n repeats (of class of interest)
# n = n total repeats in region of interest (specific mb region we categorized as speciation)
# k = n repeat (of class of interest) in (specific mb region we categorized as speciation)

# Load the functions that are are used in the loop
source("../functions.R")
window <- 5e6  # 5 Mb windows
group_levels <- c("X_chr", paste0("chr_", 1:14), "Unplaced")

# Species
species <- c("Gpenn", "Gfirm")

# List of RE Class to Test for enrichment
classes <- list(
  Gpenn = c("DNA", "MITE", "Penelope"), 
  Gfirm = c("LTR", "RC", "Satellite", "SINE")
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

## Run the main analysis loop 

# Store results
results_list <- list()

# Loop through each species
for (sp in species) {
  
  message("\nProcessing species: ", sp)
  
  # Files for this species
  rm_file <- RM_Files[[sp]]
  allele_file <- Allele_Files[[sp]]
  chr_file <- CHR_Files[[sp]]
  key_file <- KEY_Files[[sp]]
  
  # Get the repeat classes to test for this species
  species_classes <- classes[[sp]]
  
  # -----------------------------
  # Read and preprocess RM data
  # -----------------------------
  rm_df <- read_rm_out(rm_file) %>%
    mutate(
      species = sp,
      midpoint = abs((begin + end)/2),
      class = sub("/.*", "", class_family)
    )
  
  for (filter_class in species_classes) {
    message("Filtering by repeat class: ", filter_class)
    
    # -----------------------------
    # Read chromosome lengths
    # -----------------------------
    chr_lengths <- make_chr_lengths(chr_file, key_file, sp, group_levels = group_levels)
    
    # -----------------------------
    # Calculate feature density
    # -----------------------------
    feature_density <- make_feature_density(
      chr_length_df = chr_lengths,
      rm_df = rm_df,
      window = 5e6
    ) %>% 
      mutate(
      class = sub("/.*", "", class_family)
    )
    
    # -----------------------------
    # Run speciation analysis
    # -----------------------------
    analysis <- Speciation_Analysis(
      species = sp,
      window = 5e6,
      allele_file = allele_file,
      rm_file = rm_file,      
      key = key_file,
      chr_length = chr_file,
      cutoff = 0.90,
      filter_class = filter_class
    )
    
    # Extract the density tables
    AD <- analysis$allele_density
    FD <- analysis$feature_density
    
    # Join allele & feature density
    full_df <- left_join(AD, FD, by = c("species", "group", "bin", "pos"))
    
    # Hypergeometric test variables
    N <- nrow(rm_df)
    K <- nrow(rm_df %>% filter(class == filter_class))
    
    regions <- full_df %>% filter(Region_Speciation == 1) %>% select(group, bin)
    
    species_interest <- feature_density %>%
      semi_join(regions, by = c("group", "bin"))
    
    n <- sum(species_interest$n_repeats)
    k <- sum(species_interest %>% filter(class == filter_class) %>% pull(n_repeats))
    
    p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    
    # Save results in a list
    results_list[[paste0(sp, "_", filter_class)]] <- tibble(
      species = sp,
      class = filter_class,
      N = N,
      K = K,
      n = n,
      k = k,
      p_value = p_value
    )
  }
}

# Combine all results into one data frame
final_results <- bind_rows(results_list)
final_results

final_results %>% filter(p_value < 0.05)

write_tsv(final_results, file = "SpeciationRegions_Results/HypergeometricComparisons_results.tsv",)
