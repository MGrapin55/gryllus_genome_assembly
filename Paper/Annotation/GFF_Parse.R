###############################################################################################################
## Author: Michael Grapin 

## Purpose: Read, edit, and write information about our Gryllus genome annotations (gff) to a figure for our
##          paper.
###############################################################################################################
# Set Working directory 
setwd("~/Downloads/MOORE_LAB_UNL/GRYLLUS_GENOME_ASSEMBLY/GIT_REPO/gryllus_genome_assembly/Paper/Annotation")

gff = "../../Gfirm/Annotation/Gryllus_firmus_blast.gff3"
key_file = "../../Paper/Repeat_Features/Gfirm.asm.key.tsv"
outdir = "./annotation_figures"

# Libraries 
library(rtracklayer)
library(readr)
library(tidyverse)
library(Gviz)
library(GenomicRanges)
library(S4Vectors)

# Read the GFF file into a data frame
gff_object <- readGFF("Gpenn_annotation.gff3") #Useful if we net to plot the track 
#########################################################
#########################################################
read_gff_with_seqid_key <- function(gff_file, key_file) {
  # ----------------------------
  # Read scaffold key
  # ----------------------------
  key <- read_tsv(
    key_file,
    col_names = c("seqid", "original_seqid"),
    show_col_types = FALSE
  )
  
  # ----------------------------
  # Read GFF (simple tab-delimited)
  # ----------------------------
  gff_df <- read.delim(
    gff_file,
    header = FALSE,
    comment.char = "#",
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  colnames(gff_df) <- c(
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes"
  )
  
  # ----------------------------
  # Join key and reorder columns
  # ----------------------------
  final_df <- gff_df %>%
    left_join(key, by = "seqid") %>%
    select(
      original_seqid,
      source,
      type,
      start,
      end,
      score,
      strand,
      phase,
      attributes
    ) %>%
    rename(seqid = original_seqid)
  
  return(final_df)
}
###############################################
## Shared Variables 

# Chromosomes Names
chromosomes <- str_c("chr_", 1:15)
# Define desired order: chr_X first, then chr_1:14, then Unplaced
group_levels <- c("X_chr", str_c("chr_", 1:14), "Unplaced")

###############################################
# Read in the gff to a table 
Gpenn <- read_gff_with_seqid_key(gff_file = gff, key_file = key_file)
#write_tsv(Gpenn, file = "Gfirm_annotation.gff3")

# Add a grouping variable by chromosome
Gpenn <- Gpenn %>%
  mutate(
    group = case_when(
      seqid == "chr_1" ~ "X_chr",
      seqid %in% paste0("chr_", 2:15) ~ paste0("chr_", as.numeric(str_remove(seqid, "chr_")) - 1),
      TRUE ~ "Unplaced"
    ),
    group = factor(group, levels = group_levels)
  )


# group by and summarize features 
Gpenn_summary <- Gpenn %>% group_by(group, type) %>% count()
####################################################################################################




######################################################################

# Fix strand
strand_fixed <- ifelse(Gpenn$strand %in% c("+", "-"), Gpenn$strand, "*")

# Convert metadata to S4 DataFrame
metadata_df <- DataFrame(
  source     = as.character(Gpenn$source),
  type       = as.character(Gpenn$type),
  score      = as.character(Gpenn$score),
  phase      = as.character(Gpenn$phase),
  attributes = as.character(Gpenn$attributes),
  group      = as.character(Gpenn$group)
)

# Create GRanges
gff_gr <- GRanges(
  seqnames = as.character(Gpenn$seqid),
  ranges   = IRanges(start = Gpenn$start, end = Gpenn$end),
  strand   = strand_fixed,
  mcols    = metadata_df
)

# Rename metadata columns if needed
colnames(mcols(gff_gr)) <- c("source", "type", "score", "phase", "attributes", "group")


# Subset to chromosome of interest
gff_chr14 <- gff_gr[seqnames(gff_gr) == "chr_14"]  
gff_chr14 <- dropSeqlevels(gff_chr14, setdiff(seqlevels(gff_chr14), "chr_14"))


# Optionally filter by type
gff_chr14 <- gff_chr14[mcols(gff_chr14)$type == "gene"]


# Genome axis
axisTrack <- GenomeAxisTrack()



# 1. Create the Track
# stacking = "dense" forces all genes into one horizontal strip
# col = NULL removes the borders around the bars for a cleaner "band" look
bandTrack <- AnnotationTrack(gff_chr14, 
                             name = "Gene Locs", 
                             stacking = "dense", 
                             fill = "darkblue",
                             col = NULL, 
                             alpha = 0.8)

# 2. Plot
plotTracks(list(axisTrack, bandTrack), 
           from = 1, 
           to = max(end(gff_chr14)),
           main = "Gene Locations (Bands)")

plotTracks(list(axisTrack, bandTrack), 
           from = 100000, 
           to = 2000000,
           main = "Gene Locations (Bands)")
################
library(GenomicRanges)
library(Gviz)

# 1. Define the chromosome length
# Ideally, use the actual seqlength if known. If not, approximate with the last gene position.
chr_max <- max(end(gff_chr1))
# If your object already has accurate seqlengths, you can skip the line above and use:
# chr_max <- seqlengths(gff_chr1)["chr_1"]

# 2. Create bins (windows) along the chromosome (e.g., 100kb window)
# Adjust 'tilewidth' to change resolution (larger = smoother, smaller = more detailed)
bins <- tileGenome(seqlengths = c("chr_1" = chr_max), 
                   tilewidth = 100000, 
                   cut.last.tile.in.chrom = TRUE)

# 3. Count the number of genes in each bin
# This adds a "score" column to the bins object
score(bins) <- countOverlaps(bins, gff_chr1)

# Create the DataTrack
# type = "mountain" fills the area; type = "hist" makes a bar chart; type = "l" is a line
densityTrack <- DataTrack(range = bins, 
                          name = "Gene Density", 
                          type = "mountain",  
                          col.mountain = "darkblue",
                          fill.mountain = c("darkblue", "lightblue"))

# If you also want to see the individual genes below the density:
geneTrack <- AnnotationTrack(gff_chr1, name = "Genes", fill = "salmon")

plotTracks(list(axisTrack, densityTrack, geneTrack), 
           from = 1, 
           to = chr_max,
           main = "Gene Density on Chromosome 1")
