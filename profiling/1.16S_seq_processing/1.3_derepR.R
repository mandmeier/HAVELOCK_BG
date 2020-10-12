#================================================================================================
# SCRIPT 3 – derepR.R
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
load("cache/filtRs.rda")


derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepRs) <- sampleNames
save(derepRs, file = "largedata/derepRs.rda")
print("LOG ==> saved derep R")


