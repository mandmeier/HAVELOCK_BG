library("dada2")

load("largedata/derepFs_2.rda")
print("derepFs_2.rda")
load("cache/sampleNames.rda")
names(derepFs_2) <- sampleNames[1201:2400]
print(names(derepFs_2))

load("largedata/dadaFs_2.rda")
dadaFs <- sv
print("dadaFs_2.rda")
print(names(dadaFs))

load("largedata/derepRs_2.rda")
print("derepRs_2.rda")
names(derepRs_2) <- sampleNames[1201:2400]
print(names(derepRs_2))

load("largedata/dadaRs_2.rda")
dadaRs <- sv
print("dadaRs_2.rda")
print(names(dadaRs))


# merge together the inferred forward and reverse sequences
mergers_2 <- mergePairs(dadaFs, derepFs_2, dadaRs, derepRs_2, verbose=TRUE)
save(mergers_2, file = "largedata/mergers_2.rda")
