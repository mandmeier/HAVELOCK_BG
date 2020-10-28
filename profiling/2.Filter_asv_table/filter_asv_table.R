
library("phyloseq")
library("tidyverse")

load("cache/ps.rda")

#as.data.frame(tax_table(ps))

#table(burk$Species, useNA = c("always"))

#burk <- as.data.frame(tax_table(ps)) %>%
#  filter(Genus == "Burkholderia-Caballeronia-Paraburkholderia")

ps

#### 2) Kingdom bacteria or archaea ####
# sort(unique(as.character(tax_table(ps)[, 1]))) ### Kingdom
ps_ab <- subset_taxa(ps, Kingdom %in% c("Bacteria", "Archaea"))
#save(ps_ab, file = "cache/ps_ab.rda")
#sum(otu_table(ps_ab))


#### 3) remove Chloroplasts ####
# sort(unique(as.character(tax_table(ps)[, 4]))) ### Order "Chloroplast"
ps_noC <- subset_taxa(ps_ab, Order!="Chloroplast")
#save(ps_noC, file = "cache/ps_noC.rda")
#sum(otu_table(ps_noC))

#### 4) remove Mitochondria ####
# sort(unique(as.character(tax_table(ps)[, 5]))) ### family "Mitochondria"
ps_noMC <- subset_taxa(ps_noC, Family!="Mitochondria")
#save(ps_noMC, file = "cache/ps_noMC.rda")
#sum(otu_table(ps_noMC))

#### 5) retain taxa that are found in both years ####

# find common ASVs
y2018_asvs <- colSums(otu_table(subset_samples(ps_noMC, year == "Y2018")))
y2019_asvs <- colSums(otu_table(subset_samples(ps_noMC, year == "Y2019")))
y2018_asvs_lt0 <- y2018_asvs[y2018_asvs > 0]
y2019_asvs_lt0 <- y2019_asvs[y2019_asvs > 0]
common_asvs <- Reduce(intersect, list(names(y2018_asvs_lt0),names(y2019_asvs_lt0)))
length(names(y2019_asvs_lt0)[!names(y2019_asvs_lt0) %in% common_asvs])
length(names(y2018_asvs_lt0)[!names(y2018_asvs_lt0) %in% common_asvs])

# sort(unique(as.character(tax_table(ps)[, 5]))) ### family "Mitochondria"
ps_common <- prune_taxa(common_asvs, ps_noMC)
ps_common

#save(ps_common, file = "cache/ps_common.rda")
#sum(otu_table(ps_common))



load("cache/ps_common.rda")

#### 7) remove 7 samples with < 100 unique ASVs ####

## presence absence table
pa <- otu_table(ps_common)
pa[pa > 0] <- 1

unique_ASV_counts <- sort(rowSums(pa))

sparse_samples <- names(unique_ASV_counts[unique_ASV_counts < 100])
samples_to_keep <- names(unique_ASV_counts[unique_ASV_counts >= 100])

ps_short <- prune_samples(samples_to_keep, ps_common)

save(ps_short, file = "cache/ps_short.rda")
#sum(otu_table(ps_short))




#### 8) remove ASVs observed in less than 166 samples (5%, 3306*0.05) ####

pa <- otu_table(ps_short)
pa[pa > 0] <- 1

sample_counts <- sort(colSums(pa))

#length(sample_counts[sample_counts >= 166])

asvs_to_keep <- names(sample_counts[sample_counts >= 166])

ps_core <- prune_taxa(asvs_to_keep, ps_short)

save(ps_core, file = "cache/ps_core.rda")
#sum(otu_table(ps_core))




#### 9) remove taxa (genus or family) with less than 5 observed ASVs ####

# fill in family where genus unknown
taxtab <- as.data.frame(tax_table(ps_core))
taxtab$taxa <- ifelse(is.na(taxtab$Genus), paste0("f_",taxtab$Family), taxtab$Genus)
tax_table(ps_core) <- as.matrix(taxtab)

## get taxa frequency table
taxfreq <- plyr::count(taxtab, 'taxa') %>%
  arrange(desc(freq))

taxa_to_keep <- taxfreq %>%
  filter(freq >= 5)

asvs_to_keep <- taxtab %>%
  filter(taxa %in% taxa_to_keep$taxa)
asvs_to_keep <- rownames(asvs_to_keep)

ps_core5 <- prune_taxa(asvs_to_keep, ps_core)

#save(ps_core5, file = "cache/ps_core5.rda")
#sum(otu_table(ps_core5))


#load("cache/ps_core5.rda")


## get taxa counts


#### Construct phylogenetic tree from 3626 taxa ####


library("seqinr")
library("ape")

## get sequences from original table
#load("cache/seqtabNoC20.rda")
seq <- colnames(seqtabNoC20)

### make padded numbers for asvs and rename asvs (numbers instead of sequences)
asvs <- paste0("asv_", formatC(c(1:ncol(otu_table(ps))), width = 6, format = "d", flag = "0"))

### save ASVs ids and sequences for later
names(asvs) <- seq

#sequences <- asvs
#save(sequences, file = "cache/sequences.rda")

#find sequences of filtered ps object (the ones to keep)
headers <- asvs[asvs %in% taxa_names(ps_core5)]
seqs <- names(headers)

### write sequences to make tree from to fasta file:
write.fasta(as.list(seqs), headers, "cache/unaligned_ps_core5.fasta" , open = "w", nbchar = 60, as.string = FALSE)


### noe copy fasta file to HCC to run mafft, fasttree

" Shell script

#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=12:00:00
#SBATCH --mem=16gb
#SBATCH --job-name=mafttree
#SBATCH --error=./LOG/mafttree_%J.err
#SBATCH --output=./LOG/mafttree_%J.out

module load mafft/7.407
mafft cache/unaligned_ps_core5.fasta > cache/aligned_ps_core5.fasta

module load fasttree/2.1
FastTree -gtr -nt cache/aligned_ps_core5.fasta > cache/ps_core5.tre

"


# combine ps object with newly generated tree



# attach new tree file to ps object
tre <- read.tree("cache/ps_core5.tre")
ps_core5t <- merge_phyloseq(ps_core5, tre)


# fix N treatment data

sample_data(ps_core5t)[, 'nitrogen'] <- c(ifelse(sample_data(ps_core5t)[, 'row'] > 2000 & sample_data(ps_core5t)[, 'row'] < 3999, "-N", "+N"))

save(ps_core5t, file = "cache/ps_core5t.rda")











