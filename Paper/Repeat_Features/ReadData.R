# Libraries
library(tidyverse)
library(GenomicRanges)



#########################################
# 1. Function to parse the RM .out file #
#########################################
read_rm_out <- function(file_path) {
  # Read the file, skipping the header lines (usually top 3)
  # RM .out is whitespace delimited but can be irregular. 
  # read_table usually handles it well.
  rm_data <- read.table(file_path,
                        skip = 3,
                        header = FALSE,
                        sep = "",
                        stringsAsFactors = FALSE,
                        fill = TRUE,
                        quote = "")
    #read_table(file_path, skip = 3, col_names = FALSE, show_col_types = FALSE)
  
  # Assign standard RM column names
  colnames(rm_data) <- c("score", "perc_div", "perc_del", "perc_ins", 
                         "seq_id", "begin", "end", "left", 
                         "strand", "repeat_name", "class_family", 
                         "r_begin", "r_end", "r_left", "id", "star")
  
  # Filter out rows that might be empty or malformed
  rm_data <- rm_data %>% filter(!is.na(seq_id))
  
  # Chromosomes Names
  chromosomes <- str_c("chr_", 1:15)
  # Define desired order: chr_X first, then chr_1:14, then Unplaced
  group_levels <- c("X_chr", str_c("chr_", 1:14), "Unplaced")
  
  rm_data <- rm_data %>%
    mutate(
      group = case_when(
        seq_id == "chr_1" ~ "X_chr",
        seq_id %in% paste0("chr_", 2:15) ~ paste0("chr_", as.numeric(str_remove(seq_id, "chr_")) - 1),
        TRUE ~ "Unplaced"
      ),
      group = factor(group, levels = group_levels)
    )
  
  return(rm_data)
}
#-----------------------------------------------------------------------------------#
#############################
# 2. Main processing script #
#############################
get_rm_summary <- function(rm_df) {
  
  # --- A. Clean and Split Class/Family ---
  clean_df <- rm_df %>%
    mutate(
      # Split "Class/Family" into two columns. 
      # If no slash, Family = Class.
      class = str_split_fixed(class_family, "/", 2)[,1],
      family = str_split_fixed(class_family, "/", 2)[,2],
      family = ifelse(family == "", class, family),
      
      # Ensure numeric coordinates
      begin = as.numeric(begin),
      end = as.numeric(end), 
      
      # Convert RepeatMasker 'C' strand to '-' for GenomicRanges compatibility
      strand = ifelse(strand == "C", "-", "+")
    )
  
  
  # --- B. Calculate Non-Redundant (NR) Lengths using GenomicRanges ---
  # This is crucial. We cannot just sum(end - begin) because of overlaps.
  
  # Convert to GRanges object
  gr <- makeGRangesFromDataFrame(clean_df, 
                                 seqnames.field = "seq_id", 
                                 start.field = "begin", 
                                 end.field = "end",
                                 keep.extra.columns = TRUE)
  
  # 1. Total Genome-wide Masked Length (Merge ALL overlaps)
  gr_reduced <- reduce(gr)
  total_masked_genome <- sum(width(gr_reduced))
  
  # 2. Summary by Class (Merge overlaps ONLY within the same class)
  nr_by_class <- split(gr, mcols(gr)$class) %>% 
    map(~ sum(width(reduce(.)))) %>% 
    enframe(name = "class", value = "nr_length")
  
  # 3. Summary by Family
  nr_by_family <- split(gr, mcols(gr)$family) %>% 
    map(~ sum(width(reduce(.)))) %>% 
    enframe(name = "family", value = "nr_length")
  
  # --- C. Calculate Counts and Averages ---
  # While 'length' requires overlap resolution, 'count' usually refers to 
  # the number of fragments found in the file.
  
  stats_summary <- clean_df %>%
    group_by(class, family, group, species) %>%
    summarise(
      count_fragments = n(),
      avg_divergence = mean(perc_div, na.rm = TRUE),
      median_divergence = median(perc_div, na.rm = TRUE),
      stdev_divergence = stats::sd(perc_div, na.rm = TRUE),
      raw_length_bp = sum(end - begin + 1), # Length WITH overlaps (for comparison)
      .groups = "drop"
    ) %>%
    # Join with the Non-Redundant lengths calculated above
    left_join(nr_by_family, by = "family") %>%
    dplyr::rename(nr_length_bp = nr_length)
  
  return(list(
    summary_table = stats_summary,
    total_masked_bp = total_masked_genome
  ))
}
################################################################################
################################################################################
################################################################################
file_paths <- c("Gpenn.clean.out", "Gfirm.clean.out")

rm_data <- file_paths %>%
  purrr::map(function(file_path) {
    
    species <- sub("\\..*$", "", basename(file_path))
    
    s <- read_rm_out(file_path)
    s$species <- species
    results <- get_rm_summary(s)
    
    data <- as.data.frame(results$summary_table)
    
    message(
      "Total masked genome size (", species, "): ",
      results$total_masked_bp / 1e6, " Mb"
    )
    
    
    data
  }) %>%
  purrr::list_rbind()

# Write to csv
write_csv(rm_data, file = "Gpenn_Gfirm_RM_Summary.csv")

# Looking at counts 
# Summary_family
family_sum_chr <- rm_data %>%
  group_by(class, family, group, species) %>%
  summarise(count_fragments = sum(count_fragments)) %>%
  mutate(species = factor(species, levels = c("Gpenn", "Gfirm")))

family_sum <- rm_data %>%
  group_by(class, family, species) %>%
  summarise(count_fragments = sum(count_fragments)) %>%
  mutate(species = factor(species, levels = c("Gpenn", "Gfirm")))


# Summary_class
class_sum_chr <- rm_data %>%
  group_by(class, group, species) %>%
  summarise(count_fragments = sum(count_fragments)) %>%
  mutate(species = factor(species, levels = c("Gpenn", "Gfirm")))

class_sum <- rm_data %>%
  group_by(class, species) %>%
  summarise(count_fragments = sum(count_fragments)) %>%
  mutate(species = factor(species, levels = c("Gpenn", "Gfirm")))

