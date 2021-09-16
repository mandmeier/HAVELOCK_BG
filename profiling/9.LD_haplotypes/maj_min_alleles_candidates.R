#maj_min_alleles_candidates

setwd("~/Desktop/Labwork/Havelock_2019/HAVELOCK_BG")

library("tidyverse")
library("ggpubr")



hmp <- read_delim('largedata/hapmaps_v4/hapmap_chr2.txt', delim="\t", na = c("", "NA", "NN"))

hmp2 <- hmp %>%
  filter(pos > 236700000 & pos <= 237000000)


hmp <- read_delim('largedata/hapmaps_v4/hapmap_chr3.txt', delim="\t", na = c("", "NA", "NN"))

hmp3 <- hmp %>%
  filter(pos > 187210000 & pos <= 187510000)


hmp <- read_delim('largedata/hapmaps_v4/hapmap_chr5.txt', delim="\t", na = c("", "NA", "NN"))

hmp5 <- hmp %>%
  filter(pos > 188350000 & pos <= 188650000)


hmp <- read_delim('largedata/hapmaps_v4/hapmap_chr6.txt', delim="\t", na = c("", "NA", "NN"))

hmp6 <- hmp %>%
  filter(pos > 171410000 & pos <= 171710000)



hmp_candidates <- rbind(hmp2, hmp3, hmp5, hmp6)


hmp_candidates <- rbind(hmp2, hmp5)



#names <- read_csv("data/BG_MM_Gen_names.csv")

#test <- read_csv("data/microbe_counts_per_genotype.csv")

#test <- test %>%
  #rename(MM_name=genotype) %>%
  #left_join(dplyr::select(names, MM_name, GX_name))

#genotypes_w_count <- unique(test$GX_name)

#zz <- file("largedata/genotypes_w_count.txt", "wb")
#writeBin(paste(genotypes_w_count, collapse="\n"), zz ) 
#close(zz)


#genotypes_test <- read_delim("largedata/genotypes_test.txt", delim = ' ', col_names = FALSE)
#genotypes_test <- pull(genotypes_test, "X1")
#genotypes_test[genotypes_test %in% names$GX_name]
#genotypes_test[!(genotypes_test %in% names$GX_name)]
#genotypes_test[!(genotypes_test %in% names$GX_name)]
#sort(unique(names$GX_name))


save(hmp_candidates, file="cache/GWAS2/hmp_candidates.rda")


hmp_candidates <- hmp2

colnames(hmp_candidates)[1] <- "snp_id"

BG_panel <- colnames(hmp_candidates)[5: length(colnames(hmp_candidates))]


#load("largedata/GWAS2/plot_dat_v5_gt5.rda")

#chr2 <- gwas_dat %>%
  #filter(startsWith(rs, "2")) %>%
  #mutate(snp_id=paste0("S2_", ps))

#gwas_snps_chr2 <- chr2$snp_id

#colnames(hmp)[1] <- "snp_id"

#hmp_snps_chr2 <- unique(pull(hmp, snp_id))

# snps above 10 -5
snps_HN <- read_csv("largedata/GWAS2/merged_HNshort_5.csv")
snps_LN <- read_csv("largedata/GWAS2/merged_HNshort_5.csv")


### fix snp position
##snps_LN$ps <- as.numeric(sapply(strsplit(snps_LN$rs, split='-', fixed=TRUE), `[`, 2))

raw_snpids_chr2 <- snps_HN %>%
  mutate(snp_id=paste0("S2_", ps)) %>%
  filter(chr == "2")




#length(gwas_snps_chr2) # above 10 -2

length(raw_snpids_chr2$snp_id) # above 10 -5

length(hmp_snps_chr2) # all SNPs



##
#length(gwas_snps_chr2[gwas_snps_chr2 %in% hmp_snps_chr2]) # 0.0112326 % of hmp snps


length(raw_snpids_chr2$snp_id[raw_snpids_chr2$snp_id %in% hmp_snps_chr2])





names <- read_csv("data/BG_MM_Gen_names.csv")

hmp_data <- hmp_candidates %>%
  filter(snp_id %in% c("S2_236826599", "S5_188470211", "S6-171561709")) %>%
  pivot_longer(BG_panel, names_to = "genotype", values_to = "haplotype") %>%
  #dplyr::select(-starts_with("assembly"), -center, -protLSID, -assayLSID, -panelLSID, -QCcode) %>%
  rename(GX_name = genotype) %>%
  left_join(dplyr::select(names, GX_name, MM_name)) %>%
  mutate(haplotype_NN = ifelse(is.na(haplotype), "NN", haplotype)) %>%
  group_by(snp_id) %>%
  add_count(name = "total_allele_obs") %>%
  group_by(snp_id, haplotype_NN) %>%
  add_tally(name = "allele_freq")

## get NA counts to calculate allele frequency correctly
NAs <- hmp_data %>%
  dplyr::select(snp_id, GX_name, haplotype_NN, allele_freq) %>%
  filter(haplotype_NN == "NN") %>%
  ungroup() %>%
  dplyr::select(-GX_name, -haplotype_NN) %>%
  rename(NA_freq = allele_freq) %>%
  unique()


## add taxa names
load("data/group_data.rda")




## get abundance by genotype blup data

blup_HN <- read_delim("data/Blup_and_Heritability/data_blup_HN_result.txt", delim = "\t")
blup_LN <- read_delim("data/Blup_and_Heritability/data_blup_LN_result.txt", delim = "\t")


abundance_HN <- blup_HN %>%
  pivot_longer(starts_with("asv_"), names_to = "ASV", values_to = "abundance") %>%
  rename(GX_name = id) %>%
  mutate(nitrogen = "HN")


abundance_LN <- blup_LN %>%
  pivot_longer(starts_with("asv_"), names_to = "ASV", values_to = "abundance") %>%
  rename(GX_name = id) %>%
  mutate(nitrogen = "LN")


abundance_data <- rbind(abundance_HN, abundance_LN)



hmp_data_candidates <- hmp_data %>%
  left_join(NAs) %>%
  mutate(total_not_na = total_allele_obs-NA_freq) %>%
  mutate(allele = ifelse(allele_freq > total_not_na/2, "major allele", "minor allele")) %>%
  rename(rep_snp = snp_id) %>%
  rename(snp_id = rep_snp, genotype = MM_name) %>%
  left_join(abundance_data) %>%
  left_join(dplyr::select(group_data, ASV, tax_group)) %>%
  filter(haplotype %in% c("AA", "TT", "CC", "GG"))



save(hmp_data_candidates, file="cache/GWAS2/hmp_data_candidates.rda")



#load("cache/GWAS2/hmp_data_candidates.rda")



#load("data/yield_analysis/abundance_vs_phenotype.rda")


### get abundance data
#load("cache/BLUP/blup_data.rda")

#mean_blups <- blup_data %>%
  #group_by(nitrogen, GX_name) %>%
  #summarize(mean_blup = mean(blup_logrel))


### get correlation data
#load("data/yield_analysis/corr_data.rda")



#load("data/group_data.rda")


tgr <- "Sphingoaurantiacus"
snp <- "S2_236826599"
chr <- "2"
nitr <- "HN"


tgr <- "Steroidobacter"
snp <- "S5_188470211"
chr <- "5"
nitr <- "LN"





load("cache/GWAS2/hmp_data_candidates.rda")


colors <- c("major allele"="#ecb602", "minor allele"="#a45ee5")




plot_dat <- hmp_data_candidates %>%
  filter(snp_id == snp & tax_group == tgr & nitrogen == nitr)
  

p <- ggplot(plot_dat, aes(x=nitrogen, y=abundance, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  #facet_wrap(~snp_id, nrow = 1) +
  stat_compare_means(aes(group = haplotype), label = "p.signif") +
  ylab("log(microbe relative abundance)") +
  #ylab("relative microbe abundance") +
  #ylab("ASV count") +
  scale_fill_manual(values = colors) +
  #geom_point(data=NAM_tax, aes(x=nitrogen, y=log_relab, fill=allele), pch = 21, size = 4, position = position_dodge(width = 1.5)) +
  #geom_text(data=NAM_tax[NAM_tax$allele=='major allele',], aes(label=genotype), size=3, hjust=1, fontface="bold", nudge_x = -0.45) +
  #geom_text(data=NAM_tax[NAM_tax$allele=='minor allele',], aes(label=genotype), size=3, hjust=0, fontface="bold", nudge_x = 0.45) +
  #geom_point(data=NAM_tax, aes(x=nitrogen, y=log_relab, fill=allele), pch = 21, size = 4, position = position_dodge(width = 0.75)) +
  #geom_text(data=NAM_tax[NAM_tax$allele=='major allele',], aes(label=genotype), size=3, hjust=1, fontface="bold", nudge_x = -0.25) +
  #geom_text(data=NAM_tax[NAM_tax$allele=='minor allele',], aes(label=genotype), size=3, hjust=0, fontface="bold", nudge_x = 0.25) +
  
  #scale_color_manual(values = cols) +
  ggtitle(paste(tgr, "|", nitr)) +
  theme_bw() +
  #theme(axis.title.x = element_blank())
  theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())

p



### plot gene expression 7 tissue for maj and min allele

load("cache/RNAseq/MAPL_rna.rda")


"Sphingoaurantiacus"
genes <- c("Zm00001d007639", "Zm00001d007640", "Zm00001d007642", "Zm00001d007643", "Zm00001d007644", "Zm00001d007645", "Zm00001d007646", "Zm00001d007647", "Zm00001d007650", "Zm00001d007652")

"Steroidobacter"

genes <- c("Zm00001d017190", "Zm00001d017191", "Zm00001d017192", "Zm00001d017193", "Zm00001d017195", "Zm00001d017197")




gendat <- read_csv("data/RNA_seq/FPKM_Old_new_four_lines2.csv")




rnadat <- MAPL_rna %>%
  filter( gene %in% genes)

unique(rnadat$gene)



plot_dat_rna <- plot_dat %>%
  left_join(rnadat) %>%
  filter(gene %in% genes) %>%
  #filter(gene %in% c("Zm00001d007644", "Zm00001d007645", "Zm00001d007650")) %>%
  #filter(gene %in% c("Zm00001d017193")) %>%
  drop_na()





p <- ggplot(plot_dat_rna, aes(x=tissue, y=expression, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~gene, nrow = 7, scales = "free_y") +
  stat_compare_means(aes(group = haplotype), label = "p.signif") +
  ylab("gene expression") +
  #ylab("relative microbe abundance") +
  #ylab("ASV count") +
  scale_fill_manual(values = colors) +
  #geom_point(data=NAM_tax, aes(x=nitrogen, y=log_relab, fill=allele), pch = 21, size = 4, position = position_dodge(width = 1.5)) +
  #geom_text(data=NAM_tax[NAM_tax$allele=='major allele',], aes(label=genotype), size=3, hjust=1, fontface="bold", nudge_x = -0.45) +
  #geom_text(data=NAM_tax[NAM_tax$allele=='minor allele',], aes(label=genotype), size=3, hjust=0, fontface="bold", nudge_x = 0.45) +
  #geom_point(data=NAM_tax, aes(x=nitrogen, y=log_relab, fill=allele), pch = 21, size = 4, position = position_dodge(width = 0.75)) +
  #geom_text(data=NAM_tax[NAM_tax$allele=='major allele',], aes(label=genotype), size=3, hjust=1, fontface="bold", nudge_x = -0.25) +
  #geom_text(data=NAM_tax[NAM_tax$allele=='minor allele',], aes(label=genotype), size=3, hjust=0, fontface="bold", nudge_x = 0.25) +
  
  #scale_color_manual(values = cols) +
  ggtitle(paste(tgr, "|", nitr)) +
  theme_bw()
  theme(axis.title.x = element_blank())
  

p







#### plot gene expression gen



tgr <- "Sphingoaurantiacus"
snp <- "S2_236826599"
chr <- "2"
nitr <- "HN"
genes <- c("Zm00001d007639", "Zm00001d007640", "Zm00001d007642", "Zm00001d007643", "Zm00001d007644", "Zm00001d007645", "Zm00001d007646", "Zm00001d007647", "Zm00001d007650", "Zm00001d007652")



tgr <- "Steroidobacter"
snp <- "S5_188470211"
chr <- "5"
nitr <- "LN"
genes <- c("Zm00001d017190", "Zm00001d017191", "Zm00001d017192", "Zm00001d017193", "Zm00001d017195", "Zm00001d017197")




rna_dat_gen <- gendat %>%
  rename(gene = X1) %>%
  pivot_longer(-gene, names_to = "id", values_to = "expression") %>%
  separate(id, c("nitrogen", "genotype", "tissue", "rep"), "_") %>%
  filter(gene %in% genes)

unique(rna_dat_gen$gene)


gen_4genotypes <- c("B73", "SD40", "W153R", "K55")


plot_dat_gen <- hmp_data_candidates %>%
  filter(snp_id == snp & tax_group == tgr & nitrogen == nitr) %>%
  dplyr::select(-ASV,-abundance, -tax_group ) %>%
  unique() %>%
  filter(genotype %in% gen_4genotypes) %>%
  left_join(rna_dat_gen) %>%
  filter(tissue %in% c("Root", "leaf"))




p <- ggplot(plot_dat_gen, aes(x=tissue, y=expression, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~gene, nrow = 2, scales = "free_y") +
  stat_compare_means(aes(group = haplotype), label = "p.signif") +
  ylab("gene expression") +
  #ylab("relative microbe abundance") +
  #ylab("ASV count") +
  scale_fill_manual(values = colors) +
  #geom_point(data=NAM_tax, aes(x=nitrogen, y=log_relab, fill=allele), pch = 21, size = 4, position = position_dodge(width = 1.5)) +
  #geom_text(data=NAM_tax[NAM_tax$allele=='major allele',], aes(label=genotype), size=3, hjust=1, fontface="bold", nudge_x = -0.45) +
  #geom_text(data=NAM_tax[NAM_tax$allele=='minor allele',], aes(label=genotype), size=3, hjust=0, fontface="bold", nudge_x = 0.45) +
  #geom_point(data=NAM_tax, aes(x=nitrogen, y=log_relab, fill=allele), pch = 21, size = 4, position = position_dodge(width = 0.75)) +
  #geom_text(data=NAM_tax[NAM_tax$allele=='major allele',], aes(label=genotype), size=3, hjust=1, fontface="bold", nudge_x = -0.25) +
  #geom_text(data=NAM_tax[NAM_tax$allele=='minor allele',], aes(label=genotype), size=3, hjust=0, fontface="bold", nudge_x = 0.25) +
  
  #scale_color_manual(values = cols) +
  ggtitle(paste(tgr, "|", nitr)) +
  theme_bw()
theme(axis.title.x = element_blank())


p



