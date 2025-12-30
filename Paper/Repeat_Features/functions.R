# Library 
library(tidyverse)
library(GenomicRanges)
library(dplyr)
library(purrr)
library(tibble)

#########################################
# 1. Function to parse the RM .out file #
#########################################
#----------------------------------------------------------
# Function: read_rm_out
#----------------------------------------------------------
# usage: read_rm_out(<Repeat Masker TSV>)
# Accepts Repeat Masker *.out that has been parsed into TSV format. Returns a dataframe.
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

##=======================================================================================================##
#----------------------------------------
# Function: fragment_overlap_stats
#----------------------------------------
# Inputs:
#   gr          : GRanges object
#   group_cols  : character vector of metadata columns to group by (e.g., c("class") or c("family"))
# Returns:
#   tibble with raw_count, nr_count, overlap_fraction per seq_id × group
#----------------------------------------
fragment_overlap_stats <- function(gr, group_cols = c("class")) {
  
  if (!inherits(gr, "GRanges")) stop("Input must be a GRanges object")
  if (!all(group_cols %in% names(mcols(gr)))) stop("group_cols not found in GRanges metadata")
  
  # 1. Create a composite key for grouping
  # We combine seqnames and metadata into a unique string ID. 
  grp_cols_char <- as.data.frame(mcols(gr)[, group_cols, drop = FALSE])
  grouping_key <- do.call(paste, c(list(as.character(seqnames(gr))), 
                                   grp_cols_char, 
                                   sep = "|||"))
  
  # 2. Split GRanges by this key
  gr_split <- split(gr, grouping_key)
  
  # 3. Compute Counts
  nr_counts  <- map_int(gr_split, ~ length(reduce(.)))
  raw_counts <- map_int(gr_split, length)
  
  # 4. Recover Metadata Safely
  metadata_lookup <- gr %>%
    mcols() %>%
    as_tibble() %>%
    mutate(seq_id = as.character(seqnames(gr)), 
           key = grouping_key) %>%
    select(key, seq_id, all_of(group_cols)) %>%
    distinct(key, .keep_all = TRUE)
  
  # 5. Assemble Final Tibble (Overlap calculation removed)
  results <- tibble(
    key = names(gr_split),
    raw_count = raw_counts,
    nr_count  = nr_counts
  ) %>%
    left_join(metadata_lookup, by = "key") %>%
    select(-key) # Remove the helper key
  
  return(results)
}

##====================================================================================================##
#----------------------------------------------------------
# Function: compute_nr_by_group
#----------------------------------------------------------
# Inputs:
#   gr         : GRanges object
#   group_cols : character vector of metadata columns to group by (e.g., c("class") or c("family"))
# Returns:
#   tibble with seq_id, group_cols, raw_length, nr_length, overlap_fraction
#----------------------------------------------------------
compute_nr_by_group <- function(gr, group_cols) {
  
  # Check input
  if (!inherits(gr, "GRanges")) stop("Input must be a GRanges object")
  if (!all(group_cols %in% names(mcols(gr)))) stop("group_cols must exist in GRanges metadata")
  
  # Compute total genome-wide NR masked length
  total_masked_genome <- sum(unlist(
    lapply(coverage(gr), function(rle) sum(runLength(rle)[runValue(rle) > 0]))
  ))
  
  message("Total non-redundant masked genome (Mb): ", round(total_masked_genome / 1e6, 2))
  
  # -----------------------------
  # Compute NR lengths per group
  # -----------------------------
  
  # If multiple grouping columns, loop over each separately
  group_col <- group_cols[1]  # assume only one column for now
  gr_split <- split(gr, interaction(seqnames(gr), mcols(gr)[[group_col]], drop = TRUE))
  
  # NR lengths using coverage
  nr_lengths <- map_dbl(gr_split, function(x) {
    cov <- coverage(x)
    sum(unlist(lapply(cov, function(rle) sum(runLength(rle)[runValue(rle) > 0]))))
  })
  
  # Raw lengths (non-reduced)
  raw_lengths <- map_dbl(gr_split, ~ sum(width(.)))
  
  # Combine into tibble
  df <- tibble(
    grp = names(nr_lengths),
    nr_length  = nr_lengths,
    raw_length = raw_lengths
  ) %>%
    mutate(overlap_fraction = 1 - (nr_length / raw_length)) %>%
    separate(grp, into = c("seq_id", group_col), sep = "\\.")
  
  return(df)
}


##====================================================================================================##
#----------------------------------------------------------
# Function: granges_object
#----------------------------------------------------------
# Inputs:
#   rm_df : data.frame of RepeatMasker results
#           Must include: seq_id, begin, end, class_family, strand
# Returns:
#   GRanges object with cleaned columns and extra metadata
#----------------------------------------------------------
granges_object <- function(rm_df) {
  
  # Check required columns
  required_cols <- c("seq_id", "begin", "end", "class_family", "strand")
  if (!all(required_cols %in% colnames(rm_df))) {
    stop("rm_df must contain columns: seq_id, begin, end, class_family, strand")
  }
  
  # 1. Clean data
  clean_df <- rm_df %>%
    mutate(
      # Split "Class/Family" into separate columns
      class  = str_split_fixed(class_family, "/", 2)[, 1],
      family = str_split_fixed(class_family, "/", 2)[, 2],
      family = ifelse(family == "", class, family),
      
      # Ensure numeric coordinates
      begin = as.numeric(begin),
      end   = as.numeric(end),
      
      # Standardize strand
      strand = ifelse(strand == "C", "-", "+")
    )
  
  # 2. Convert to GRanges
  gr <- makeGRangesFromDataFrame(
    clean_df,
    seqnames.field = "seq_id",
    start.field    = "begin",
    end.field      = "end",
    strand.field   = "strand",
    keep.extra.columns = TRUE
  )
  
  return(gr)
}
##===============================================================================================##


