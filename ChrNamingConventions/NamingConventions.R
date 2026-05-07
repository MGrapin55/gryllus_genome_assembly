
library(tidyr)
# Read in the columns of interest and save to csv 

a <- read_tsv("G.assimilis_SeqReport.tsv")

a <- a %>% select(`GenBank seq accession`, `Chromosome name`, `Sequence name`, `Seq length`) %>% mutate(species = "assimilis", 
                                                                                                        key = "WGS") %>% slice(1:15)



b <- read_tsv("G.bimaculatus_SeqReport.tsv")
b <- b %>% select(`GenBank seq accession`, `Chromosome name`, `Sequence name`, `Seq length`) %>% mutate(species = "bimaculatus", 
                                                                                                        key = "WGS") %>% slice(1:15)


# Gpenn and Gfirm Naming 
p <- read_tsv("../Paper/Supplemental/Histograms/Gpenn.final.tsv", 
              col_names = c("Sequence name", "Seq length" )) %>% slice(1:15) %>% mutate(species = "pennsylvanicus", 
                                                                                        key = "WGS", 
                                                                                        `GenBank seq accession` = NA, 
                                                                                        `Chromosome name` = NA)

f <- read_tsv("../Paper/Supplemental/Histograms/Gfirm.final.tsv", 
              col_names = c("Sequence name", "Seq length" )) %>% slice(1:15) %>% mutate(species = "firmus", 
                                                                                    key = "WGS", 
                                                                                    `GenBank seq accession` = NA, 
                                                                                    `Chromosome name` = NA)



combined <- bind_cols(a, b, f, p)

write_csv(combined, "combined_wide.csv")


