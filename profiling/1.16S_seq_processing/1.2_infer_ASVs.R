#================================================================================================
# SCRIPT 2 – infer_ASVs.R
#================================================================================================


#================================================================================================
### infer sequence variants from forward and reverse fastq reads
### derep is fast but takes a lot of memory. used 96G for 576 samples
### mergePairs takes a few hrs probably
#================================================================================================

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
load("cache/filtRs.rda")

# dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

save(derepFs, file = "cache/derepFs.rda")
save(derepRs, file = "cache/derepRs.rda")
print("LOG ==> saved derep")

# filter sequence variants that are probably PCR or sequencing errors and not biological
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

save(errF, file = "cache/errF.rda")
save(errR, file = "cache/errR.rda")
print("LOG ==> saved err")

# independent inference by sample
dadaFs <- dada(derepFs, err=errF, pool = FALSE, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, pool = FALSE, multithread=TRUE)

save(dadaFs, file = "cache/dadaFs.rda")
save(dadaRs, file = "cache/dadaRs.rda")
print("LOG ==> saved dada")

# merge together the inferred forward and reverse sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
save(mergers, file = "largedata/mergers.rda")
print("LOG ==> saved mergers")

