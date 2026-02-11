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
#' Speciation Analysis from SNP and RepeatMasker Data
#'
#' This function integrates SNP density and RepeatMasker feature density
#' across fixed genomic windows to identify putative speciation regions.
#' Genomic windows whose SNP density exceeds a user-defined quantile
#' (default = 95%) within each chromosome are classified as
#' \emph{speciation-enriched}.
#'
#' The function automatically:
#' \itemize{
#'   \item Standardizes chromosome naming conventions
#'   \item Bins SNPs and repeats into fixed genomic windows
#'   \item Computes per-window densities
#'   \item Calculates chromosome-specific quantile cutoffs
#'   \item Classifies windows as speciation or non-speciation
#' }
#'
#' @param species Character. Species name (e.g. "Gpenn", "Gfirm").
#' @param window Numeric. Genomic window size (e.g. 1e4, 1e5, 1e6).
#' @param allele_file Character. Path to SNP/allele TSV file
#'   with columns: \code{seqid, location}.
#' @param rm_file Character. Path to RepeatMasker .out file
#'   parsed using \code{read_rm_out()}.
#' @param chr_length Character. Path to chromosome lengths TSV file 
#'    with columns: \code{seqid, Length_BP}
#'    @param key Character. Path to chromosome names key file 
#'    with columns: \code{Renamed, Orginal}

#' @param cutoff Numeric. Quantile cutoff for SNP density
#'   (default = 0.95).
#'
#' @return A named list of class \code{"SpeciationAnalysis"} with:
#' \describe{
#'   \item{allele_density}{SNP density per genomic window}
#'   \item{feature_density}{RepeatMasker feature density per window}
#'   \item{cutoffs}{Chromosome-specific SNP density thresholds}
#'   \item{rm_df}{Annotated RepeatMasker data}
#'   \item{snps}{Annotated SNP data}
#' }
#'
#' @examples
#' \dontrun{
#' analysis <- Speciation_Analysis(
#'   species = "Gpenn",
#'   window = 1e6,
#'   allele_file = "alleles.tsv",
#'   rm_file = "repeatmasker.out", 
#'   chr_length = "Gpenn.lengths.tsv",
#'   key = "Gpenn.key.tsv"
#' )
#' }

Speciation_Analysis <- function(
    species,
    window = 1e6,
    allele_file,
    rm_file,
    chr_length,
    key,
    cutoff = 0.95,
    filter_class = NULL,
    filter_class_family = NULL
) 
  {
  # ------------------------------
  # Input validation
  # ------------------------------
  if (!is.character(species) || length(species) != 1) {
    stop("`species` must be a single character string (e.g. 'Gpenn' or 'Gfirm').")
  }
  
  if (!is.numeric(window) || window <= 0) {
    stop("`window` must be a positive numeric value (e.g. 1e4, 1e5, 1e6).")
  }
  
  if (!is.numeric(cutoff) || cutoff <= 0 || cutoff >= 1) {
    stop("`cutoff` must be a numeric quantile between 0 and 1 (e.g. 0.95).")
  }
  
  if (!is.character(allele_file) || !file.exists(allele_file)) {
    stop("`allele_file` must be a valid path to a TSV file with columns: seqid, location.")
  }
  
  if (!is.character(rm_file) || !file.exists(rm_file)) {
    stop("`rm_file` must be a valid path to a RepeatMasker .out file readable by read_rm_out().")
  }
  
  if (!is.character(chr_length) || !file.exists(chr_length)) {
    stop("`chr_length` must be a valid path to a TSV file with columns: seqid, Length_BP")
  }
  
  if (!is.character(key) || !file.exists(key)) {
    stop("`key` must be a valid path to a TSV file with columns: Renamed, Orginal")
  }
  
  
  message("Running Speciation Analysis for: ", species)
  message("Window size: ", window)
  message("Quantile cutoff: ", cutoff)
  
  
  # ------------------------------
  # Load required libraries
  # ------------------------------
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
  })
  
  # ------------------------------
  # Load RepeatMasker data
  # ------------------------------
  rm_df <- read_rm_out(rm_file)
  
  rm_df <- rm_df %>%
    mutate(
      species = species,
      midpoint = abs((begin + end) / 2),
      class = sub("/.*", "", class_family)
    )
  
  # if flag filter_class == value Run this block else skips it filter_class flag is off by default 
  if (!is.null(filter_class)) {
    message("Filtering by repeat class: ", filter_class)
    rm_df <- rm_df %>%
      filter(class == filter_class)
  }
  
  if (!is.null(filter_class_family)) {
    message("Filtering by repeat class/family: ", filter_class_family)
    rm_df <- rm_df %>%
      filter(class == filter_class_family)
  }
  
  # ------------------------------
  # Load SNP data
  # ------------------------------
  DM <- read_tsv(allele_file, col_names = FALSE, show_col_types = FALSE)
  colnames(DM) <- c("seqid", "location")
  
  DM <- DM %>%
    mutate(species = species)
  
  # ------------------------------
  # Load Chromosome Length Data
  # ------------------------------
  group_levels <- c("X_chr", paste0("chr_", 1:14), "Unplaced")

  # Read inputs
  chr_df <- read_tsv(chr_length,
                     col_names = c("Renamed", "Length_BP"),
                     show_col_types = FALSE)
  
  key_df <- read_tsv(key, show_col_types = FALSE)
  
  # Join and keep only what we need
  chr_df <- left_join(chr_df, key_df, by = "Renamed")
  chr_df <- chr_df[, c("Orginal", "Length_BP")]
  
  # Keep first 15 rows
  chr_df <- chr_df[1:15, ]
  
  # Rename for consistency
  names(chr_df)[names(chr_df) == "Orginal"] <- "group"
  
  # Remap chromosome names
  chr_df$group <- ifelse(chr_df$group == "chr_1", "X_chr", chr_df$group)
  
  idx <- chr_df$group %in% paste0("chr_", 2:15)
  chr_df$group[idx] <- paste0(
    "chr_",
    as.numeric(str_remove(chr_df$group[idx], "chr_")) - 1
  )
  
  chr_df$group[!chr_df$group %in% group_levels] <- "Unplaced"
  
  # Apply factor order
  chr_df$group <- factor(chr_df$group, levels = group_levels)
  
  # Order by chromosome
  chr_length_df <- chr_df[order(chr_df$group), ]
  
  
  # ------------------------------
  # Chromosome re-mapping
  # ------------------------------
  DM <- DM %>%
    mutate(
      group = case_when(
        seqid == "chr_1" ~ "X_chr",
        seqid %in% paste0("chr_", 2:15) ~
          paste0("chr_", as.numeric(str_remove(seqid, "chr_")) - 1),
        TRUE ~ "Unplaced"
      ),
      group = factor(group, levels = group_levels)
    )
  
  rm_df <- rm_df %>%
    mutate(
      group = factor(group, levels = group_levels)
    )
  
  # ------------------------------
  # SNP density per window
  # ------------------------------
  all_windows <- chr_length_df %>%
    mutate(
      bin = map(Length_BP, ~ seq(0, .x - 1, by = window))  #Fixes the last window (i.e 10 -> 15) not being a 0 because no genomic content
    ) %>%
    unnest(bin) %>%
    select(group, bin)
  
  snp_bins <- DM %>%
    mutate(bin = floor(location / window) * window) %>%
    count(species, group, bin, name = "n_snps")
  
  allele_density <- expand_grid(
    species = unique(DM$species),
    group   = unique(chr_length_df$group)
  ) %>%
    left_join(all_windows, by = "group") %>%
    left_join(snp_bins, by = c("species", "group", "bin")) %>%
    mutate(
      n_snps = replace_na(n_snps, 0),
      pos = bin / window
    )
  
  # ------------------------------
  # Quantile cutoffs per chromosome
  # ------------------------------
  cutoffs <- allele_density %>%
    group_by(group) %>%
    summarise(
      cutoff = quantile(n_snps, cutoff, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ------------------------------
  # Speciation classification
  # ------------------------------
  allele_density <- allele_density %>%
    left_join(cutoffs, by = "group") %>%
    mutate(
      Region_Speciation = if_else(n_snps > cutoff, 1L, 0L)        # was >= vs just > (Matters for the visual)
    )
  
  # ------------------------------
  # Feature density (RepeatMasker)
  # ------------------------------
  RM_bins <- rm_df %>%
    filter(group %in% group_levels[1:15]) %>%
    mutate(bin = floor(midpoint / window) * window) %>%
    count(species, group, bin, name = "n_repeats")
  
  feature_density <- expand_grid(
    species = unique(rm_df$species),
    group   = unique(chr_length_df$group)
  ) %>%
    left_join(all_windows, by = "group") %>%
    left_join(RM_bins, by = c("species", "group", "bin")) %>%
    mutate(
      n_repeats = replace_na(n_repeats, 0),
      pos = bin / window
    )
  
  
  # feature_density <- rm_df %>%
  #   filter(group %in% group_levels[1:15]) %>%
  #   mutate(bin = floor(midpoint / window) * window) %>%
  #   group_by(group, bin) %>%
  #   summarise(
  #     density = n() / (window / window),
  #     .groups = "drop"
  #   ) %>%
  #   mutate(pos = bin / window)
  
  # ------------------------------
  # Output object
  # ------------------------------
  analysis <- list(
    species = species,
    window = window,
    cutoff = cutoff,
    allele_density = allele_density,
    feature_density = feature_density,
    cutoffs = cutoffs,
    rm_df = rm_df,
    snps = DM,
    lengths = chr_length_df
  )
  
  class(analysis) <- "SpeciationAnalysis"
  
  return(analysis)
}

##===============================================================================================##
#' Plot Genome-wide Speciation and Feature Density Tracks
#'
#' This function visualizes the output of \code{Speciation_Analysis()} as
#' a two-track genome plot. The lower track encodes speciation-enriched
#' windows (based on SNP density cutoffs) as colored tiles, while the
#' upper track displays RepeatMasker feature density as a continuous line.
#'
#' Each chromosome (or scaffold group) is shown in a separate facet,
#' allowing direct comparison of genomic structure and speciation signal
#' across the genome.
#'
#' @param analysis A \code{SpeciationAnalysis} object returned by
#'   \code{Speciation_Analysis()}.
#' @param speciation_colors Character vector of length 2 specifying colors
#'   for non-speciation (0) and speciation (1) windows.
#' @param line_width Numeric. Line width for the feature density track
#'   (default = 0.8).
#'
#' @return A \code{ggplot} object containing the genome track visualization.
#'
#' @examples
#' \dontrun{
#' analysis <- Speciation_Analysis(
#'   species = "Gpenn",
#'   window = 1e6,
#'   allele_file = "alleles.tsv",
#'   rm_file = "repeatmasker.out"
#' )
#'
#' p <- plot_speciation_tracks(analysis)
#' print(p)
#' }



plot_speciation_tracks <- function(analysis,
                                   speciation_colors = c("grey80", "red"),
                                   line_width = 0.8) {
  
  # ------------------------------
  # Input checks
  # ------------------------------
  if (!inherits(analysis, "SpeciationAnalysis")) {
    stop("Input must be an object returned by Speciation_Analysis().")
  }
  
  if (!all(c("allele_density", "feature_density") %in% names(analysis))) {
    stop("SpeciationAnalysis object is missing required slots.")
  }
  
  allele_density  <- analysis$allele_density
  feature_density <- analysis$feature_density
  
  # ------------------------------
  # Math setup for tracks
  # ------------------------------
  max_y <- max(feature_density$n_repeats, na.rm = TRUE)
  
  if (!is.finite(max_y)) {
    stop("Feature density contains no finite values.")
  }
  
  track_height <- max_y * 0.1
  track_y_pos  <- -(track_height / 2)
  
  # ------------------------------
  # Plot
  # ------------------------------
  p <- ggplot() +
    
    # --- Track 1: Speciation windows ---
    geom_tile(
      data = allele_density,
      aes(x = pos,
          y = track_y_pos,
          fill = as.factor(Region_Speciation)),
      height = track_height,
      width  = 1
    ) +
    
    # --- Track 2: Feature density ---
    geom_line(
      data = feature_density,
      aes(x = pos, y = n_repeats),
      linewidth = line_width
    ) +
    
    # --- Formatting ---
    scale_fill_manual(
      values = speciation_colors,
      name = "Speciation"
    ) +
    
    facet_wrap(~group, scales = "free_x") +
    coord_cartesian(ylim = c(-track_height, NA)) +
    labs(
      x = "Genomic position (windows)",
      y = "Number of Repeats Per Window",
      title = paste0("Speciation tracks: ", analysis$species),
      subtitle = paste0("Window = ", format(analysis$window, scientific = TRUE),
                        " | Cutoff = ", analysis$cutoff)
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.spacing = unit(0.8, "lines"),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    )
  
  return(p)
}
##===============================================================================================##
# Function: make_chr_lengths
# Reads in a tsv of seqid, Length_BP and formats it to a df used in Speciation Analysis() function.
make_chr_lengths <- function(chr_length, key, prefix, group_levels) {
  
  library(dplyr)
  library(readr)
  library(stringr)
  
  # Read inputs
  chr_df <- read_tsv(chr_length,
                     col_names = c("Renamed", "Length_BP"),
                     show_col_types = FALSE)
  
  key_df <- read_tsv(key, show_col_types = FALSE)
  
  # Join and keep only what we need
  chr_df <- left_join(chr_df, key_df, by = "Renamed")
  chr_df <- chr_df[, c("Orginal", "Length_BP")]
  
  # Keep first 15 chromosomes
  chr_df <- chr_df[1:15, ]
  
  # Rename for consistency
  names(chr_df)[names(chr_df) == "Orginal"] <- "group"
  
  # Remap chromosome names
  chr_df$group <- ifelse(chr_df$group == "chr_1", "X_chr", chr_df$group)
  
  idx <- chr_df$group %in% paste0("chr_", 2:15)
  chr_df$group[idx] <- paste0(
    "chr_",
    as.numeric(str_remove(chr_df$group[idx], "chr_")) - 1
  )
  
  chr_df$group[!chr_df$group %in% group_levels] <- "Unplaced"
  
  # Apply factor order
  chr_df$group <- factor(chr_df$group, levels = group_levels)
  
  # Add species prefix
  chr_df$species <- prefix
  
  # Order by chromosome
  chr_df <- chr_df[order(chr_df$group), ]
  
  return(chr_df)
}

# This function could be cleaned up it it just fits more module rather than writing all this out. 

##===============================================================================================##
# Function to make repeat feature densities per species along the whole chromosome
make_feature_density <- function(chr_length_df, rm_df, window) {
  
  # ------------------------------
  # 1) Build windows per chromosome (and species if present)
  # ------------------------------
  all_windows <- chr_length_df %>%
    mutate(
      bin = purrr::map(Length_BP, ~ seq(0, .x - 1, by = window))
    ) %>%
    tidyr::unnest(bin) %>%
    dplyr::select(species, group, bin)
  
  # ------------------------------
  # 2) Bin repeat midpoints
  # ------------------------------
  RM_bins <- rm_df %>%
    mutate(
      midpoint  = (begin + end) / 2,
      length_bp = abs(end - begin),
      bin       = floor(midpoint / window) * window
    ) %>%
    group_by(species, group, bin, class_family) %>%
    summarise(
      n_repeats    = n(),
      n_repeats_bp = sum(length_bp),
      .groups = "drop"
    )
  
  # ------------------------------
  # 3) Join repeats into windows
  # ------------------------------
  feature_density <- all_windows %>%
    left_join(RM_bins,
              by = c("species", "group", "bin")) %>%
    tidyr::replace_na(list(n_repeats = 0)) %>%
    mutate(pos = bin / window)
  
  return(feature_density)
}
##===============================================================================================##
