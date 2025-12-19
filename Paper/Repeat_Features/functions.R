# Library 
library(tidyverse)
library(GenomicRanges)

#########################################
# 1. Function to parse the RM .out file #
#########################################
# usage: read_rm_out(<Repeat Masker TSV>)
# Accepts Repeat Masker *.out that has been parsed into TSV format. 
# Columns:
#  ("score", "perc_div", "perc_del", "perc_ins", 
#   "seq_id", "begin", "end", "left", 
#   "strand", "repeat_name", "class_family", 
#   "r_begin", "r_end", "r_left", "id", "star")

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