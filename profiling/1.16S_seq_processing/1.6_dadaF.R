library("dada2")

load("largedata/derepFs.rda")
load("cache/errF.rda")


derepFs_1 <-derepFs[1:1200]
save(derepFs_1, file = "largedata/derepFs_1.rda")
rm(derepFs_1)

derepFs_2 <-derepFs[1201:2400]
save(derepFs_2, file = "largedata/derepFs_2.rda")
rm(derepFs_2)

derepFs_3 <-derepFs[2401:length(derepFs)]
save(derepFs_3, file = "largedata/derepFs_3.rda")
rm(derepFs_3)


dadaFs <- dada(derepFs, err=errF, pool = FALSE, multithread=TRUE)


sv <- dadaFs[1:1200]
save(sv, file = "largedata/dadaFs_1.rda")
sv <- dadaFs[1201:2400]
save(sv, file = "largedata/dadaFs_2.rda")
sv <- dadaFs[2401:3313]
save(sv, file = "largedata/dadaFs_3.rda")

save(dadaFs, file = "largedata/dadaFs.rda")

print("LOG ==> saved dada F")

