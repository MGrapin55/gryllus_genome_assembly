library(stringr)
library(dplyr)
library(readr)

txt <- read_lines("collinearity.txt")

# Keep only alignment blocks
blocks <- txt[grepl("^## Alignment|^\\s+\\d+-", txt)]

res <- list()
current <- NULL

for (l in blocks) {
  
  # Block header
  if (grepl("^## Alignment", l)) {
    
    m <- str_match(l,
                   "Alignment (\\d+):.* (\\S+)&(\\S+) (plus|minus)")
    
    current <- list(
      block = as.integer(m[2]),
      chr1  = m[3],
      chr2  = m[4],
      orient = m[5]
    )
  }
  
  # Gene pair lines
  if (grepl("^\\s+\\d+-", l)) {
    
    m <- str_match(l,
                   "\\d+-\\s*\\d+:\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)")
    
    res[[length(res) + 1]] <- data.frame(
      block = current$block,
      chr1  = current$chr1,
      chr2  = current$chr2,
      orient = current$orient,
      gene1 = m[2],
      gene2 = m[3],
      evalue = as.numeric(m[4]),
      stringsAsFactors = FALSE
    )
  }
}

syn <- bind_rows(res)
head(syn)


syn <- syn %>%
  mutate(
    gene1 = sub(".*\\|", "", gene1),
    gene2 = sub(".*\\|", "", gene2),
    chr1  = sub(".*\\|", "", chr1),
    chr2  = sub(".*\\|", "", chr2)
  )


library(GenomicRanges)
library(rtracklayer)

gff <- import("Gryllus_firmus.gff3")
genes <- gff[gff$type == "gene"]

pos <- data.frame(
  gene = genes$ID,
  chr = as.character(seqnames(genes)),
  start = start(genes),
  end = end(genes)
)



syn2 <- syn %>%
  left_join(pos, by = c("gene1" = "gene")) %>%
  rename(start1 = start, end1 = end) %>%
  left_join(pos, by = c("gene2" = "gene")) %>%
  rename(start2 = start, end2 = end)


library(gggenomes)

gggenomes(
  seqs = pos,
  links = syn2,
  seq_id = chr,
  start = start,
  end = end,
  link_id = block
) +
  geom_seq() +
  geom_link(alpha = 0.25)

