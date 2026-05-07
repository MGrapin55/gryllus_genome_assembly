source("../functions.R")
# Reading in the data with a function read_rm_out() from functions.R 
p <- read_rm_out("../Gpenn.clean.out")

# Add a species column 
p$species <- "Gpenn"


# Get the midpoint value
p <- p %>%
  mutate(midpoint = abs((begin + end) / 2),
         class = sub("/.*", "", class_family)
  )

analysis_DNA <- Speciation_Analysis(
  species = "Gpenn",
  window = 5e6,
  allele_file = "../pennsylvanicus.DMs.txt",
  rm_file = "../Gpenn.clean.out",
  key = "../Gpenn.asm.key.tsv", 
  chr_length = "../Gpenn.chr.lengths.tsv",
  cutoff = 0.90,
  filter_class = "LTR"
)

p <- plot_speciation_tracks(
  analysis_DNA,
  speciation_colors = c("grey80", "red"),
  line_width = 0.8
)

print(p)

# Extract the orginal df's from the slots and join to make a full_df
AD <- analysis_DNA$allele_density
FD <- analysis_DNA$feature_density

full_df <- left_join(AD, FD, by = c("species", "group",   "bin",     "pos"))

# Total rows of RE features
N = nrow(p)
# Total rows of RE features of specific class 
K = p %>% filter(class == "LTR") %>% nrow

# Join the regions of speciation
Gpenn_regions <- full_df %>%
  filter(Region_Speciation == 1) %>%
  select(group, bin)

Gpenn_feature_density <- Gpenn_feature_density %>% mutate(
  class = sub("/.*", "", class_family)
)
Gpenn_interest <- Gpenn_feature_density %>%
  semi_join(Gpenn_regions, by = c("group", "bin"))

n = Gpenn_interest %>% pull(n_repeats) %>% sum()

k = Gpenn_interest %>% filter(class == "LTR") %>% pull(n_repeats) %>% sum()

p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
p_value
