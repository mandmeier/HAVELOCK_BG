library("dada2")

load("largedata/derepFs_1.rda")
print("derepFs_1.rda")
load("cache/sampleNames.rda")
names(derepFs_1) <- sampleNames[1:1200]
print(names(derepFs_1))

load("largedata/dadaFs_1.rda")
dadaFs <- sv
print("dadaFs_1.rda")
print(names(dadaFs))

load("largedata/derepRs_1.rda")
print("derepRs_1.rda")
names(derepRs_1) <- sampleNames[1:1200]
print(names(derepRs_1))

load("largedata/dadaRs_1.rda")
dadaRs <- sv
print("dadaRs_1.rda")
print(names(dadaRs))


# merge together the inferred forward and reverse sequences
mergers_1 <- mergePairs(dadaFs, derepFs_1, dadaRs, derepRs_1, verbose=TRUE)
save(mergers_1, file = "largedata/mergers_1.rda")
