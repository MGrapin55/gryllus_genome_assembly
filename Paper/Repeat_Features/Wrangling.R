
### Purpose: Preprocesses data to calculate the difference between overlapping and non overlapping 
###          Repeat Masker repeative features
###

# set working dirctory 
setwd("~/Downloads/MOORE_LAB_UNL/GRYLLUS_GENOME_ASSEMBLY/GIT_REPO/gryllus_genome_assembly/Paper/Repeat_Features")

# Load the functions into the global environment 
source("functions.R")

# Set file to work on
file_path = "Gpenn.clean.out"
species <- "Gpenn"
############################################################################################

# Read in the RM data
data <- read_rm_out(file_path)


# Create a Granges object 
gr <- granges_object(data)


# -----------------------------
# Fragment overlap stats
# -----------------------------
# For class
frag_class_df <- fragment_overlap_stats(gr, group_cols = c("class"))

# For family
frag_family_df <- fragment_overlap_stats(gr, group_cols = c("family"))


# -----------------------------
# NR basepair overlap stats
# -----------------------------
# For class
nr_class_df <- compute_nr_by_group(gr, group_cols = c("class"))
nr_class_df <- rename(nr_class_df, "overlap_fraction" = "overlap_fraction_class")

# For family
nr_family_df <- compute_nr_by_group(gr, group_cols = c("family"))
nr_family_df <- rename(nr_family_df, "overlap_fraction" = "overlap_fraction_family")


# Join all this data together
class_df <- frag_class_df %>%
  left_join(
    nr_class_df,
    by = c("seq_id", "class")
  )

family_df <- frag_family_df %>%
  left_join(
    nr_family_df,
    by = c("seq_id", "family")
  )

order <- c("seq_id","class", "family", "raw_count", "nr_count", "raw_length", "nr_length", "overlap_fraction")


# -----------------------------
# 1. Split class / family in RM data
# -----------------------------
data <- data %>%
  mutate(
    class  = str_split_fixed(class_family, "/", 2)[, 1],
    family = str_split_fixed(class_family, "/", 2)[, 2],
    family = if_else(family == "", class, family)
  )

# -----------------------------
# 2. Rename summary columns to avoid collisions
# -----------------------------
# Prefix class-level columns
class_df <- class_df %>%
  rename_with(
    ~ paste0(.x, "_class"),
    c("raw_count", "nr_count", "raw_length", "nr_length")
  )

# Prefix family-level columns
family_df <- family_df %>%
  rename_with(
    ~ paste0(.x, "_family"),
    c("raw_count", "nr_count", "raw_length", "nr_length")
  )


# -----------------------------
# 3. Join class- and family-level stats
# -----------------------------
master_df <- data %>%
  left_join(class_df,  by = c("seq_id", "class")) %>%
  left_join(family_df, by = c("seq_id", "family"))

# Write to a csv
filename <- paste(species, "RM_Summary_Master.csv")
write_csv(master_df, filename = filename)


