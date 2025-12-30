

# Some questions we have to test 

# 1.) Type of Data we have to test 
# perc_div == numeric,
# raw_count == count
# nr_count == count
# overlap_fraction == numeric
# nr_length == numeric, continious

species <- c("Gpenn", "Gfirm")
chromosomes <- c("X_chr", paste0("chr_", 1:14), "unplaced")

# Repeat levels
class_groups <- paste0("class_", 1:11)          # replace with your actual class names
family_groups <- paste0("family_", 1:95)        # replace with actual family names

# -----------------------------
# 1. Global comparisons (all chromosomes)
# -----------------------------
global_class <- expand.grid(
  species1 = "Gpenn",
  species2 = "Gfirm",
  chromosome = "all",
  repeat_category = class_groups,
  repeat_level = "class"
)

global_family <- expand.grid(
  species1 = "Gpenn",
  species2 = "Gfirm",
  chromosome = "all",
  repeat_category = family_groups,
  repeat_level = "class/family"
)

# -----------------------------
# 2. Chromosome-level comparisons
# -----------------------------
chr_class <- expand.grid(
  species1 = "Gpenn",
  species2 = "Gfirm",
  chromosome = chromosomes,
  repeat_category = class_groups,
  repeat_level = "class"
)

chr_family <- expand.grid(
  species1 = "Gpenn",
  species2 = "Gfirm",
  chromosome = chromosomes,
  repeat_category = family_groups,
  repeat_level = "class/family"
)

# -----------------------------
# 3. Combine all
# -----------------------------
all_comparisons <- rbind(global_class, global_family, chr_class, chr_family)



results <- map_dfr(seq_len(nrow(all_comparisons)), function(i) {
  
  comp <- all_comparisons[i, ]
  
  # Subset data for this comparison
  df_sub <- master_df %>%
    filter(species %in% c(comp$species1, comp$species2)) %>%
    filter(chromosome == comp$chromosome | comp$chromosome == "all") %>%
    filter(
      (comp$repeat_level == "class" & class == comp$repeat_category) |
        (comp$repeat_level == "class/family" & family == comp$repeat_category)
    )
  
  # Skip if not enough data
  if(nrow(df_sub) < 2) return(NULL)
  
  # Select metric to test
  metric <- if(comp$repeat_level == "class") "nr_length_class" else "nr_length_family"
  
  # Run Wilcoxon test (non-parametric)
  test_res <- wilcox.test(df_sub[[metric]] ~ df_sub$species)
  
  # Return results as a tibble
  tibble(
    species1 = comp$species1,
    species2 = comp$species2,
    chromosome = comp$chromosome,
    repeat_level = comp$repeat_level,
    repeat_category = comp$repeat_category,
    p_value = test_res$p.value
  )
})
