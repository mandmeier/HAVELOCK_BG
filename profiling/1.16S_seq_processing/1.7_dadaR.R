library("dada2")

load("largedata/derepRs.rda")
load("cache/errR.rda")


derepRs_1 <-derepRs[1:1200]
save(derepRs_1, file = "largedata/derepRs_1.rda")
rm(derepRs_1)

derepRs_2 <-derepRs[1201:2400]
save(derepRs_2, file = "largedata/derepRs_2.rda")
rm(derepRs_2)

derepRs_3 <-derepRs[2401:length(derepRs)]
save(derepRs_3, file = "largedata/derepRs_3.rda")
rm(derepRs_3)


dadaRs <- dada(derepRs, err=errR, pool = FALSE, multithread=TRUE)


sv <- dadaRs[1:1200]
save(sv, file = "largedata/dadaRs_1.rda")
sv <- dadaRs[1201:2400]
save(sv, file = "largedata/dadaRs_2.rda")
sv <- dadaRs[2401:3313]
save(sv, file = "largedata/dadaRs_3.rda")


save(dadaRs, file = "largedata/dadaRs.rda")

print("LOG ==> saved dada R")
