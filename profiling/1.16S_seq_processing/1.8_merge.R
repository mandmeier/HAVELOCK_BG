library("dada2")

load("largedata/derepFs.rda")
load("largedata/dadaFs.rda")

load("largedata/derepRs.rda")
load("largedata/dadaRs.rda")


# merge together the inferred forward and reverse sequences
mergers <- mergePairs(dadaRs, derepRs, dadaRs, derepRs, verbose=TRUE)
save(mergers, file = "largedata/mergers.rda")
print("LOG ==> saved mergers")




