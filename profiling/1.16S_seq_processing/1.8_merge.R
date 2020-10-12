library("dada2")

load("largedata/derepFs.rda")
load("largedata/dadaFs.rda")

load("largedata/derepRs.rda")
load("largedata/dadaRs.rda")


# merge together the inferred forward and reverse sequences
mergers <- mergePairs(dadaRs, derepRs, dadaRs, derepRs)
save(mergers, file = "largedata/mergers.rda")
print("LOG ==> saved mergers")


seqtabAll <- makeSequenceTable(mergers)
save(seqtabAll, file = "largedata/seqtabAll.rda")
print("LOG ==> saved seqtabAll")

# remove chimeric sequences by comparing each inferred sequence to the others in the table
seqtabNoC <- removeBimeraDenovo(seqtabAll)
save(seqtabNoC, file = "largedata/seqtabNoC.rda")
print("LOG ==> saved seqtabNoC")


