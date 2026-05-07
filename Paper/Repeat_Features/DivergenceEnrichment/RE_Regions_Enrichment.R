################################################################################################################
##
##                              Gpenn vs Gfirm regions of divergence overlapping with RE regions
##                                        Author: Michael Grapin 
################################################################################################################

# Outline: 
# The goal is to find regions of speciation (opposite fixed alleles) and compare RE (repeat element) features
# if the are more or less overlapping with speciation windows

# Scales: (window)
# 10 kb 
# 100 kb 
# 1 mb 

# Calculating SNP (opposite fixed alleles per window) then plotting with kernel density estimator function. Obtaining 
# 5% cutoff for binary classification threshold as window is different (contribute to speciation) or not different 
# no contribution. 

# Overlay RE elements with windows of speciation. Compute Hypergeometric distribution to to see if features are 
# statistically enriched. 

# LIBRARIES: 
library(tidyverse)
library(ggplot2)

# SETTING PARAMETERS: 
setwd("~/Downloads/MOORE_LAB_UNL/GRYLLUS_GENOME_ASSEMBLY/GIT_REPO/gryllus_genome_assembly/Paper/Repeat_Features/DivergenceEnrichment")

out_dir <- "SpeciationRegions_Results/Correct"

# Scales (window size)
#window = 1e4 # 10 kb
#window = 1e5 # 100 kb 
#window = 1e6 # 1 mb 
#window_sizes <- c(1e4, 1e5, 1e6, 1e7, 1e8)
#window_sizes <- c(1e6, 5e6, 1e7, 2e7, 3e7, 4e7, 5e7, 7e7, 7e7, 8e7, 9e7, 1e7)
#window_sizes <- c(5e6, 1e7, 1.5e7, 2e7)
window_sizes <- c(5e6)

# Experiment with window size 5mb = 5e6, 50mb = 50e6

species <- c("Gpenn", "Gfirm")

RM_Files <- c(
  Gpenn = "../Gpenn.clean.out",
  Gfirm = "../Gfirm.clean.out"
)
Allele_Files <- c(
  Gpenn = "../pennsylvanicus.DMs.txt",
  Gfirm = "../firmus.DMs.txt"
)

CHR_Files <- c(
  Gpenn = "../Gpenn.chr.lengths.tsv",
  Gfirm = "../Gfirm.chr.lengths.tsv"
)

KEY_Files <- c(
  Gpenn = "../Gpenn.asm.key.tsv",
  Gfirm = "../Gfirm.asm.key.tsv"
)
cutoff_quantile = 0.90
  
# Load some handy functions in
source(file = "../functions.R")

# Add other parameters as necessary.... 
###########################################################################################################
dir.create(out_dir, showWarnings = FALSE)

# ------------------------------
# Loop over species
# ------------------------------
for (sp in species) {
  
  pdf_file <- file.path(
    out_dir,
    paste0(sp, "_", cutoff_quantile,"Quantile", "_SpeciationTracks.pdf")
  )
  
  pdf(pdf_file, width = 11.69, height = 8.27)
  
  # ------------------------------
  # Loop over window sizes
  # ------------------------------
  for (w in window_sizes) {
    
    analysis <- Speciation_Analysis(
      species = sp,
      window = w,
      allele_file = Allele_Files[[sp]],
      rm_file = RM_Files[[sp]],
      chr_length = CHR_Files[[sp]],
      key = KEY_Files[[sp]],
      cutoff = cutoff_quantile
    )
    
    p <- plot_speciation_tracks(
      analysis,
      speciation_colors = c("grey80", "red"),
      line_width = 0.8
    )
    
    print(p)
  }
  
  dev.off()
}



## STATISTICAL TESTING 
# Manuelly done for one
analysis <- Speciation_Analysis(
  species = "Gpenn",
  window = 5e6,
  allele_file = "../pennsylvanicus.DMs.txt",
  rm_file = "../Gpenn.clean.out",
  key = "../Gpenn.asm.key.tsv", 
  chr_length = "../Gpenn.chr.lengths.tsv",
  cutoff = 0.90
)

p <- plot_speciation_tracks(
  analysis,
  speciation_colors = c("grey80", "red"),
  line_width = 0.8
)

print(p)

AD <- analysis$allele_density
FD <- analysis$feature_density

full_df <- left_join(AD, FD, by = c("species", "group",   "bin",     "pos"))



AD %>% group_by(group) %>% distinct(cutoff)
# # A tibble: 15 × 2
# # Groups:   group [15]
# group   cutoff
# <fct>    <dbl>
#   1 X_chr  226.   
# 2 chr_1    0    
# 3 chr_2   36.4  
# 4 chr_3    0    
# 5 chr_4    0    
# 6 chr_5    0    
# 7 chr_6    0    
# 8 chr_7    0    
# 9 chr_8    0.800
# 10 chr_9    1.2  
# 11 chr_10   2.6  
# 12 chr_11 249.   
# 13 chr_12   1.8  
# 14 chr_13 262.   
# 15 chr_14   0    
# Zeros mean it is really not concentrated so it is probably worth while to just excluded them


library(MASS)

# 1. Fit a Negative Binomial Model (handles overdispersion)
# If window sizes vary, add + offset(log(window_size))
fit_nb <- glm.nb(repeat_count ~ Region_Speciation, data = df)

summary(fit_nb)

# 2. Check the coefficients
# Exp(coef) tells you the fold-change
exp(coef(fit_nb))



library(MASS)

# --- Step 1: Fit the Real Model ---
# We use glm.nb to account for overdispersion (variance > mean)
real_model <- glm.nb(repeat_count ~ Region_Speciation, data = df)
obs_coef   <- coef(real_model)["Region_Speciation"]

print(paste("Observed Coefficient:", round(obs_coef, 4)))

# --- Step 2: Run Permutations ---
set.seed(42)       # For reproducibility
n_perms <- 1000    # 1000 is a standard starting point
perm_coefs <- numeric(n_perms)

# Using a progress bar so you know it's working
pb <- txtProgressBar(min = 0, max = n_perms, style = 3)

for(i in 1:n_perms) {
  
  # A. Create a temporary dataframe
  df_perm <- df
  
  # B. Shuffle the Speciation labels ONLY. 
  # This breaks the link between Speciation and Repeats.
  df_perm$Region_Speciation <- sample(df$Region_Speciation)
  
  # C. Run the model on this shuffled "nonsense" data
  # We use try() because sometimes models fail to converge on random noise
  try({
    # use quiet=TRUE to suppress convergence warnings during the loop
    perm_fit <- glm.nb(repeat_count ~ Region_Speciation, data = df_perm)
    perm_coefs[i] <- coef(perm_fit)["Region_Speciation"]
  }, silent = TRUE)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# --- Step 3: Calculate the Empirical P-Value ---
# Remove any failed runs (NAs) if the model didn't converge
perm_coefs <- perm_coefs[!is.na(perm_coefs)]

# Two-tailed test: How often was the random coefficient 'more extreme' 
# (positive or negative) than our observed one?
p_value <- mean(abs(perm_coefs) >= abs(obs_coef))

# --- Step 4: Visualize Results ---
hist(perm_coefs, breaks = 30, col = "grey", border = "white",
     main = "Permutation Test Distribution",
     xlab = "Coefficient under Null Hypothesis")
abline(v = obs_coef, col = "red", lwd = 2, lty = 2) # Your real result

cat("\n------------------------------------------------\n")
cat("Empirical P-value:", p_value, "\n")
if(p_value < 0.05) {
  cat("Result: SIGNIFICANT. The relationship is likely real.\n")
} else {
  cat("Result: NOT SIGNIFICANT. This could be due to random chance.\n")
}
cat("------------------------------------------------\n")



library(MASS)
library(dplyr)

# --- Helper Function: Circular Shift ---
# This function rotates a vector 'x' by 'n' positions.
# If x = [1, 2, 3, 4, 5] and n = 2, result = [4, 5, 1, 2, 3]
circular_shift <- function(x, n) {
  if (n == 0) return(x)
  len <- length(x)
  n <- n %% len # Handle shifts larger than length
  c(tail(x, n), head(x, len - n))
}

# --- Step 1: Fit the Real Model ---
# Using Negative Binomial (glm.nb) as discussed
real_model <- glm.nb(repeat_count ~ Region_Speciation, data = df)
obs_coef   <- coef(real_model)["Region_Speciation"]

print(paste("Observed Coefficient:", round(obs_coef, 4)))

# --- Step 2: Run Circular Permutations ---
set.seed(123)
n_perms <- 1000
perm_coefs <- numeric(n_perms)

# Progress bar
pb <- txtProgressBar(min = 0, max = n_perms, style = 3)

for(i in 1:n_perms) {
  
  # A. Create a shifted version of the dataframe
  # We group by chromosome so the shift wraps around per-chromosome,
  # preventing data from 'leaking' between different chromosomes.
  df_perm <- df %>%
    mutate(
      # Calculate a random shift size for THIS chromosome
      shift_amount = sample(1:n(), 1),
      
      # Apply the shift to the Speciation column only
      Region_Speciation_Perm = circular_shift(Region_Speciation, first(shift_amount))
    ) %>%
    ungroup()
  
  # B. Run the model on the shifted data
  # Note: We regress against the NEW permuted column 'Region_Speciation_Perm'
  try({
    perm_fit <- glm.nb(repeat_count ~ Region_Speciation_Perm, data = df_perm)
    perm_coefs[i] <- coef(perm_fit)["Region_Speciation_Perm"]
  }, silent = TRUE)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# --- Step 3: Calculate P-Value ---
# Filter out NAs (failed convergences)
perm_coefs <- perm_coefs[!is.na(perm_coefs)]

# Calculate Empirical P-value (Two-tailed)
p_value <- mean(abs(perm_coefs) >= abs(obs_coef))

# --- Step 4: Visualize ---
hist(perm_coefs, breaks = 30, col = "lightblue", border = "white",
     main = "Circular Permutation Distribution",
     xlab = "Coefficient under Null (Shifted Data)",
     xlim = range(c(perm_coefs, obs_coef))) # Ensure real coef is visible

abline(v = obs_coef, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Null Distribution", "Observed"), 
       fill = c("lightblue", "red"), density = c(NA, 40))

cat("\n------------------------------------------------\n")
cat("Circular Permutation P-value:", p_value, "\n")
cat("------------------------------------------------\n")
###########################################################################################################
#                                   Manuel Coded Process
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
# # LOADING DATA: 
# # Reading in the data with a function read_rm_out() from functions.R 
# p <- read_rm_out("../Gpenn.clean.out")
# f <- read_rm_out("../Gfirm.clean.out")
# 
# # Add a species column 
# p$species <- "Gpenn"
# f$species <- "Gfirm"
# 
# # Combine dataframes
# rm_df <- bind_rows(p, f)
# ##=======================================##
# # SNP Data 
# DM_f <- read_tsv("../firmus.DMs.txt", col_names = FALSE)
# colnames(DM_f) <- c("seqid", "location")
# 
# DM_g <- read_tsv("../pennsylvanicus.DMs.txt", col_names = FALSE)
# colnames(DM_g) <- c("seqid", "location")
# 
# DM_g$species <- "Gpenn" 
# DM_f$species <- "Gfirm" 
# 
# DM <- bind_rows(DM_g, DM_f)
# 
# group_levels <- c("X_chr", str_c("chr_", 1:14), "Unplaced")
# 
# DM <- DM %>%
#   mutate(
#     group = case_when(
#       seqid == "chr_1" ~ "X_chr",
#       seqid %in% paste0("chr_", 2:15) ~ paste0("chr_", as.numeric(str_remove(seqid, "chr_")) - 1),
#       TRUE ~ "Unplaced"
#     ),
#     group = factor(group, levels = group_levels))
# 
# # SNP DENSITY
# # Allele density normalized per Mb
# allele_density <- DM %>%
#   mutate(bin = floor(location / window) * window) %>%
#   group_by(species, group, bin) %>%
#   summarise(
#     density = n() / (window / window),
#     .groups = "drop"
#   ) %>%
#   mutate(pos = bin / window)
# 
# 
# # PLOT SNP DENSITY 
# # Visualize it
# ggplot(
#   allele_density,
#   aes(
#     x = density,
#     y = species,
#     fill = species
#   )
# ) +
#   geom_density_ridges(
#     scale = 1,
#     rel_min_height = 0.01,
#     alpha = 0.7,
#     color = NA 
#     ) + 
#   stat_density_ridges(
#     quantile_lines = TRUE, 
#     quantiles = c(0.95),
#     alpha = 0.7
#   ) +
#   labs(
#     x = "Genomic position",
#     y = "Species"
#   ) + facet_wrap(~group)
# 
# 
# # CALCULATE CUTOFF
# # 5% cuffoff
# #cutoff <- quantile(x, 0.95, na.rm = TRUE)
# 
# cutoffs <- allele_density %>%
#   group_by(species, group) %>%
#   summarise(
#     cutoff_95 = quantile(density, 0.95, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# 
# # CLASSIFY WINDOW 
# # Speciation or Not
# allele_density <- allele_density %>%
#   left_join(cutoffs, by = c("species", "group"))
# 
# allele_density <- allele_density %>%
#   mutate(
#     Region_Speciation = if_else(density >= cutoff_95, 1L, 0L)
#   )
# 
# 
# # PLOT RE ELEMENTS
# # Remove class_family values that aren't in both 
# rm_df <- rm_df %>%
#   group_by(class_family) %>%
#   filter(n_distinct(species) == 2) %>%
#   ungroup()
# 
# # Set factor levels and get new columns
# rm_df <- rm_df %>%
#   mutate(species = factor(species, levels = c("Gpenn", "Gfirm")),
#          midpoint = abs((begin + end) / 2),
#          midpoint_mb = midpoint / 1e6)
# 
# # Feature density summed across all class_family
# feature_density <- rm_df %>%
#   filter(group %in% group_levels[1:15]) %>%
#   mutate(bin = floor(midpoint / window) * window) %>%
#   group_by(species, group, bin) %>%
#   summarise(
#     density = n() / (window / window),
#     .groups = "drop"
#   ) %>%
#   mutate(pos = bin / window)
# 
#   # 1. MATH SETUP
#   max_y <- max(feature_density$density, na.rm = TRUE)
#   track_height <- max_y * 0.1 
#   track_y_pos  <- -(track_height / 2)
#   
#   # 2. PLOT
#   ggplot() +
#     # --- Track 1: Speciation ---
#     geom_tile(
#       data = allele_density,
#       # FIX IS HERE: Use as.factor() to tell R these numbers are categories
#       aes(x = pos, y = track_y_pos, fill = as.factor(Region_Speciation)),
#       height = track_height, 
#       width = 1
#     ) +
#     
#     # --- Track 2: Feature Density ---
#     geom_line(
#       data = feature_density,
#       aes(x = pos, y = density),
#       linewidth = 0.8
#     ) +
#     
#     # --- Formatting ---
#     scale_fill_manual(
#       values = c("grey80", "red"), 
#       name = "Speciation",
#       # Optional: labels = c("No", "Yes") if you want to rename 0/1 in the legend
#     ) +
#     
#     facet_wrap(~ species + group, scales = "free_x") +
#     coord_cartesian(ylim = c(-track_height, NA)) +
#     labs(x = "Genomic position (Mb)", y = "Density") +
#     theme_bw() +
#     theme(legend.position = "top")
# 
# 
# # STATISTIC TEST
# # speciation windows
# # RE windows
# # equally distributed or not equally distributed
# 
# 
# 
# 
# 
