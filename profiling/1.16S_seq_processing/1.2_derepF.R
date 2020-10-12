#================================================================================================
# SCRIPT 2 – derepF.R
#================================================================================================

### derep is fast but takes a lot of memory. used 96G for 576 samples

library("knitr")
library("BiocStyle")
library("ggplot2")
library("gridExtra")
library("dada2")
library("phyloseq")
library("DECIPHER")
library("phangorn")

load("cache/sampleNames.rda")
load("cache/filtFs.rda")

# dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
names(derepFs) <- sampleNames
save(derepFs, file = "largedata/derepFs.rda")
print("LOG ==> saved derep F")
