library("phyloseq")
library("dplyr")

### ASV table
load("cache/seqtabNoC20.rda")
### taxonomy table
load("cache/taxTab.rda")
### metadata
samdf <- read.csv("data/BG_sample_data.csv", header = TRUE)
rownames(samdf) <- samdf$Sample_ID ## use sampleID as rownames, otherwise ps won't accept it

#=======================================================================
#### make phyloseq object
#=======================================================================

ps <- phyloseq(otu_table(seqtabNoC20, taxa_are_rows=FALSE), 
  sample_data(samdf), 
  tax_table(taxTab))

### make padded numbers for asvs and rename asvs (numbers instead of sequences)
asvnum <- paste0("asv_", formatC(c(1:ncol(otu_table(ps))), width = 6, format = "d", flag = "0"))

### remember sequences for later
names(asvnum) <- taxa_names(ps)

save(asvnum, file = "data/sequences.rda")

## use asv names instead of sequences
taxa_names(ps) <- asvnum

save(ps, file = "cache/ps.rda") ## save phyloseq object


