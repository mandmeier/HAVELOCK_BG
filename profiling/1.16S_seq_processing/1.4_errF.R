#================================================================================================
# SCRIPT 4 – errF.R
#================================================================================================

library("knitr")
library("BiocStyle")
library("ggplot2")
library("gridExtra")
library("dada2")
library("phyloseq")
library("DECIPHER")
library("phangorn")

load("cache/filtFs.rda")
# filter sequence variants that are probably PCR or sequencing errors and not biological
errF <- learnErrors(filtFs, multithread=TRUE)
save(errF, file = "cache/errF.rda")
print("LOG ==> saved errF")

