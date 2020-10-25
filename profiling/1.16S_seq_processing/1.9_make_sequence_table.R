library("dada2")

## combine mergers

load("largedata/mergers_1.rda")
load("largedata/mergers_2.rda")
load("largedata/mergers_3.rda")

mergers <- c(mergers_1, mergers_2, mergers_3)

seqtabAll <- makeSequenceTable(mergers)
save(seqtabAll, file = "largedata/seqtabAll.rda")
print("LOG ==> saved seqtabAll")

# remove chimeric sequences by comparing each inferred sequence to the others in the table
seqtabNoC <- removeBimeraDenovo(seqtabAll, verbose=TRUE)
save(seqtabNoC, file = "largedata/seqtabNoC.rda")
print("LOG ==> saved seqtabNoC")

