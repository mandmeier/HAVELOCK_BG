library("dada2")

load("largedata/seqtabAll.rda")

## remove ASVs with fewer than 20 total observations
seqtabAll20 <- seqtabAll[,colSums(seqtabAll) >= 20]

save(seqtabAll20, file = "largedata/seqtabAll20.rda")
print("LOG ==> saved seqtabAll20")

# remove chimeric sequences by comparing each inferred sequence to the others in the table
seqtabNoC <- removeBimeraDenovo(seqtabAll20, verbose=TRUE)
save(seqtabNoC, file = "largedata/seqtabNoC.rda")
print("LOG ==> saved seqtabNoC")
