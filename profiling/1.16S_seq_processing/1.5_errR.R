#================================================================================================
# SCRIPT 5 – errR.R
#================================================================================================

library("knitr")
library("BiocStyle")
library("ggplot2")
library("gridExtra")
library("dada2")
library("phyloseq")
library("DECIPHER")
library("phangorn")

load("cache/filtRs.rda")
# filter sequence variants that are probably PCR or sequencing errors and not biological
errR <- learnErrors(filtRs, multithread=TRUE)
save(errR, file = "cache/errR.rda")
print("LOG ==> saved errR")

