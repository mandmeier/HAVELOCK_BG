library("dada2")

load("largedata/derepFs_3.rda")
print("derepFs_3.rda")
load("cache/sampleNames.rda")
names(derepFs_3) <- sampleNames[2401:3313]
print(names(derepFs_3))

load("largedata/dadaFs_3.rda")
dadaFs <- sv
print("dadaFs_3.rda")
print(names(dadaFs))

load("largedata/derepRs_3.rda")
print("derepRs_3.rda")
names(derepRs_3) <- sampleNames[2401:3313]
print(names(derepRs_3))

load("largedata/dadaRs_3.rda")
dadaRs <- sv
print("dadaRs_3.rda")
print(names(dadaRs))


# merge together the inferred forward and reverse sequences
mergers_3 <- mergePairs(dadaFs, derepFs_3, dadaRs, derepRs_3, verbose=TRUE)
save(mergers_3, file = "largedata/mergers_3.rda")
