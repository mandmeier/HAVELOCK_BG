# alpha and beta diversity


library("tidyverse")
library("phyloseq")
library("ggplot2")
library("officer")
library("rvg")
library("vegan")



### import final phyloseq object
load("cache/ps_noMC.rda")

# use ps object with shared ASVs (5% SAMPLES before) filtering 150 most abundant groups

#### 6) remove ASVs observed in less than 166 samples (5%, 3313*0.05) ####

pa <- otu_table(ps_noMC)
pa[pa > 0] <- 1

sample_counts <- sort(colSums(pa))

#length(sample_counts[sample_counts >= 166])

asvs_to_keep <- names(sample_counts[sample_counts >= 166])

ps_common <- prune_taxa(asvs_to_keep, ps_noMC)




## remove 6 samples with < 1000 ASVs

## HAV_5505 HAV_6038 HAV_2574 HAV_5774 HAV_2501 HAV_5773

unique_ASV_counts <- sort(rowSums(otu_table(ps_common)))

samples_to_exclude <- names(unique_ASV_counts[unique_ASV_counts < 1000])

samples_to_keep <- names(unique_ASV_counts[unique_ASV_counts >= 1000])

ps_common <- prune_samples(samples_to_keep, ps_common)


save(ps_common, file = "cache/ps_common")


### make phylogenetic tree



#### Construct phylogenetic tree from 4632 taxa ####


library("seqinr")
library("ape")


#sequences <- asvs
#save(sequences, file = "cache/sequences.rda")

load("cache/sequences.rda")

head(sequences)

#find sequences of filtered ps object (the ones to keep)
headers <- sequences[sequences %in% taxa_names(ps_common)]

seqs <- names(headers)

### write sequences to make tree from to fasta file:
write.fasta(as.list(seqs), headers, "cache/unaligned_ps_common.fasta" , open = "w", nbchar = 60, as.string = FALSE)


### now copy fasta file to HCC to run mafft, fasttree

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
mafft cache/unaligned_ps_common.fasta > cache/aligned_ps_common.fasta

module load fasttree/2.1
FastTree -gtr -nt cache/aligned_ps_common.fasta > cache/ps_common.tre

"


# combine ps object with newly generated tree



# attach new tree file to ps object
tre <- read.tree("cache/ps_common.tre")
ps_common <- merge_phyloseq(ps_common, tre)


# fix N treatment data

sample_data(ps_common)
sample_data(ps_common)[, 'nitrogen'] <- c(ifelse(sample_data(ps_common)[, 'row'] > 2000 & sample_data(ps_common)[, 'row'] < 3999, "-N", "+N"))

save(ps_common, file = "cache/ps_common.rda")



#### Venn Diagram show shared taxa 2018 vs 2019 ####

# find common ASVs
y2018_asvs <- colSums(otu_table(subset_samples(ps_common, year == "Y2018")))
y2019_asvs <- colSums(otu_table(subset_samples(ps_common, year == "Y2019")))

y2018_asvs_lt0 <- y2018_asvs[y2018_asvs > 0]
y2019_asvs_lt0 <- y2019_asvs[y2019_asvs > 0]

common_asvs <- Reduce(intersect, list(names(y2018_asvs_lt0),names(y2019_asvs_lt0)))


#### 3728 ASVs are common to both years, 904 only found in 2019


## plot shannon diversity between 4632 (common_asvs) and 3728 sets (ps_bothyears)

ps_bothyears <- prune_taxa(common_asvs, ps_common)



rch_4632 <- estimate_richness(ps_common, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "Fisher"))
rch_4632$asv_set <- "all_ASVs"
median(rch_4632$Shannon) # 6.434921

rch_3728 <- estimate_richness(ps_bothyears, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "Fisher"))
rch_3728$asv_set <- "ASVs_common_to_both_years"
median(rch_3728$Shannon) # 6.287561
diversity <- rbind(rch_4632, rch_3728)



library("ggpubr")

colors = c("#743282", "#dc9200")

shannon_plot <- ggplot(diversity, aes(x = asv_set, y = Shannon)) +
  geom_jitter(aes(color = asv_set), size=1, width=0.2) +
  geom_boxplot(outlier.shape = NA) +
  #stat_compare_means(aes(group = asv_set), label = "p.format") +
  scale_color_manual(values=colors) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

shannon_plot


(1-median(rch_3728$Shannon)/median(rch_4632$Shannon)) *100 ## 2.290005 % diversity lost

(1-median(rch_3728$Observed)/median(rch_4632$Observed)) *100 ## 12.65823 % diversity lost

3728/4632



ps_common_family <- tax_glom(ps_common, taxrank="Family")

ps_bothyears_family <- tax_glom(ps_bothyears, taxrank="Family")




rch_4632_fam <- estimate_richness(ps_common_family, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "Fisher"))
rch_4632_fam$asv_set <- "all_ASVs"
rch_3728_fam <- estimate_richness(ps_bothyears_family, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "Fisher"))
rch_3728_fam$asv_set <- "ASVs_common_to_both_years"
diversity_fam <- rbind(rch_4632_fam, rch_3728_fam)

# quartz()

shannon_plot_fam <- ggplot(diversity_fam, aes(x = asv_set, y = Simpson)) +
  geom_jitter(aes(color = asv_set), size=1, width=0.2) +
  geom_boxplot(outlier.shape = NA) +
  #scale_color_manual(values=colors) +
  theme_minimal()

shannon_plot_fam

(1-median(rch_3728_fam$Shannon)/median(rch_4632_fam$Shannon)) *100 ## 4.570855 % diversity lost



?geom_jitter


divplot <- ggplot(diversity, aes(x= asv_set, y = Shannon )) +
  geom_bo
  

divplot



length(names(y2019_asvs_lt0)[!names(y2019_asvs_lt0) %in% common_asvs])
length(names(y2018_asvs_lt0)[!names(y2018_asvs_lt0) %in% common_asvs])





load("cache/ps_common.rda")



#### Constrained ordination ####

## https://rdrr.io/rforge/vegan/man/adonis.html

# convert to relative abundance, aka Total Sum Normalization [TSS] and log transform ASV counts
ps <- transform_sample_counts(ps_common, function(x) x / sum(x))


sdat <- data.frame(sample_data(ps))


# exclude checks
ps <- subset_samples(ps, !genotype %in% c("CHECK", "CHECK2"))




## get distance matrix: this takes a long time!
#dm <- distance(ps, method = "wunifrac")


### replace missing values with median
#dm[is.na(dm)] <- median(dm, na.rm = TRUE)
#save(dm, file='largedata/distance_matrix_wunifrac.rda')


load("largedata/distance_matrix_wunifrac.rda")

## calculate bray-custis distance matrix
bray.cap.whole <- capscale(as.dist(dm) ~ year + genotype + nitrogen + block + sp + spb, data = sdat, add = T, na.action = na.exclude)


#save(bray.cap.whole, file='largedata/bray.cap.whole.rda')

#anova(bray.cap.whole, by = "terms")

bray.cap.whole.axes <- data.frame(cbind(sdat, scores(bray.cap.whole)$sites))

## Permutation ANOVA (this takes 24h or so...)
#permanova <- adonis(dm~year + genotype + nitrogen + block + sp + spb,data = sdat,add = T)

#save(permanova, file = "largedata/permanova.rda")

load("largedata/permanova.rda")

permanova

percent_explained <- bray.cap.whole$CCA$eig / sum(bray.cap.whole$CCA$eig) * 100
percent_explained

colors <- c("#ff0000", "#0000ff")



#find genotypes common to both years, without check
sdat <- data.frame(sample_data(ps_asv))
y2018 <- sdat %>%
  filter(year == "Y2018")
y2019 <- sdat %>%
  filter(year == "Y2019")
common_genotypes <- Reduce(intersect, list(unique(y2018$genotype),unique(y2019$genotype)))
common_genotypes <- common_genotypes[-1]



bray.cap.whole.axes$color_groups <- ifelse(bray.cap.whole.axes$genotype %in% common_genotypes, "common_HN", bray.cap.whole.axes$nitrogen)

bray.cap.whole.axes$color_groups <- ifelse(bray.cap.whole.axes$genotype %in% common_genotypes & bray.cap.whole.axes$nitrogen == "-N", "common_LN", bray.cap.whole.axes$color_groups)


bray.cap.whole.axes$shape_groups <- ifelse(bray.cap.whole.axes$genotype %in% common_genotypes, "common", "unique")

bray.cap.whole.axes$shape_groups <- factor(bray.cap.whole.axes$shape_groups, levels = c("unique", "common"))

colors <- c("-N"="#ff0000","common_LN"="#C00001", "+N"="#0000ff", "common_HN"="#000080")

colors <- c("#ff0000", "#0000ff")

colors = c("#ffc30b", "#000000")


ord_plot <- ggplot(bray.cap.whole.axes, aes(x = CAP1, y = CAP2, color = shape_groups, shape = nitrogen)) +
  geom_point(size=1, alpha = 0.5) +
  facet_wrap(~ year, nrow = 1) +
  scale_color_manual(values=colors) +
  labs(x = "Constrained PCo1 (31.80%)", y = "Constrained PCo2 (26.24%)") +
  theme_bw()

ord_plot    
