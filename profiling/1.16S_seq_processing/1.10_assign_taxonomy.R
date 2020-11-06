library("dada2")

# load("/work/jyanglab/mmeier/HAVELOCK_BG/cache/seqtabNoC.rda")

## remove ASVs with fewer than 20 total observations
# seqtabNoC20 <- seqtabNoC[,colSums(seqtabNoC) >= 20]

load("/work/jyanglab/mmeier/HAVELOCK_BG/cache/seqtabNoC20.rda")

# save(seqtabNoC20, file = "cache/seqtabNoC20.rda")
# print("LOG ==> saved seqtabNoC20")

### Assign taxonomy

## filepath to SILVA train set for DADA 2 https://zenodo.org/record/1172783#.XhN6UxdKh24
fastaRef <- "/work/jyanglab/mmeier/HAVELOCK_BG/largedata/silva_nr99_v138_wSpecies_train_set.fa.gz"
taxTab <- assignTaxonomy(seqtabNoC20, refFasta = fastaRef, multithread=TRUE, verbose=TRUE)
unname(head(taxTab))
## save taxTab
save(taxTab, file = "cache/taxTab.rda")
