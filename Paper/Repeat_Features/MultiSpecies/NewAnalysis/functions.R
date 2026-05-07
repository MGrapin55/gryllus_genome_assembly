## Functions 

##==== parse_repeat_summary ====## 
# parses summary file from Repeat Masker buildSummary.pl
# file_path = path to summary file 
# species_name = species group name
parse_repeat_summary <- function(file_path, species_name) {
  # Read the contents of the file
  lines <- readLines(file_path)
  
  # 1. Extract total length for the 'total_length_bp' column
  total_length_line <- grep("^Total Length:", lines, value = TRUE)[1]
  total_length_bp <- as.numeric(gsub("[^0-9]", "", total_length_line))
  
  # 2. Find the boundaries of the "Repeat Classes" table
  header_idx <- grep("^Class\\s+Count\\s+bpMasked", lines)
  start_idx <- header_idx + 1 # Target the separator under the header
  
  # The table ends just before the line starting with "Total "
  total_indices <- grep("^Total\\s+[0-9]+", lines)
  end_idx <- total_indices[total_indices > start_idx][1]
  
  # Extract relevant lines
  table_lines <- lines[(start_idx + 1):(end_idx - 1)]
  
  # Clean up lines: remove separators, the 'total interspersed' summary, and blanks
  table_lines <- table_lines[!grepl("^-+$", table_lines)]
  table_lines <- table_lines[!grepl("total interspersed", table_lines, ignore.case = TRUE)]
  table_lines <- table_lines[trimws(table_lines) != ""]
  
  # 3. Initialize vectors to build the dataframe
  class_vec <- character()
  class_family_vec <- character()
  count_fam_vec <- character()
  bp_fam_vec <- character()
  pct_fam_vec <- character()
  
  # Track state while looping
  curr_class <- NA
  
  # Iterate over lines to extract classes and families
  for (i in seq_along(table_lines)) {
    line <- table_lines[i]
    
    # Check if line is a subclass/family (indicated by leading whitespace)
    if (grepl("^\\s+", line)) {
      parts <- strsplit(trimws(line), "\\s+")[[1]]
      fam_name <- parts[1]
      
      class_vec <- c(class_vec, curr_class)
      class_family_vec <- c(class_family_vec, paste0(curr_class, "/", fam_name))
      
      count_fam_vec <- c(count_fam_vec, parts[2])
      bp_fam_vec <- c(bp_fam_vec, parts[3])
      pct_fam_vec <- c(pct_fam_vec, parts[4])
      
    } else {
      # Line is a main Class (left-aligned with no leading space)
      parts <- strsplit(trimws(line), "\\s+")[[1]]
      curr_class <- parts[1]
      
      # Add the parent class as its own row with "none" as the family
      class_vec <- c(class_vec, curr_class)
      class_family_vec <- c(class_family_vec, paste0(curr_class, "/none")) 
      
      # Map the parent's stats directly into the family columns
      count_fam_vec <- c(count_fam_vec, parts[2])
      bp_fam_vec <- c(bp_fam_vec, parts[3])
      pct_fam_vec <- c(pct_fam_vec, parts[4])
    }
  }
  
  # Helper to clean missing values ("--") and percentages ("%"), then convert to numeric
  clean_num <- function(x) {
    x <- gsub("%", "", x)
    x[!is.na(x) & x == "--"] <- NA
    as.numeric(x)
  }
  
  # 4. Construct the final dataframe
  df <- data.frame(
    species = species_name,
    class = class_vec,
    class_family = class_family_vec,
    count_per_class_family = clean_num(count_fam_vec),
    BpMasked_per_class_family = clean_num(bp_fam_vec),
    PercentMasked_per_class_family = clean_num(pct_fam_vec),
    total_length_bp = total_length_bp,
    stringsAsFactors = FALSE
  )
  
  return(df)
}
