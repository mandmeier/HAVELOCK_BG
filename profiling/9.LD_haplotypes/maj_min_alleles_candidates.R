#maj_min_alleles_candidates


library("tidyverse")


hmp <- read_delim('largedata/Hapmaps/hapmap_chr2.txt', delim="\t", na = c("", "NA", "NN"))

hmp2 <- hmp %>%
  filter(pos > 236690000 & pos <= 236990000)


hmp <- read_delim('largedata/Hapmaps/hapmap_chr3.txt', delim="\t", na = c("", "NA", "NN"))

hmp3 <- hmp %>%
  filter(pos > 187210000 & pos <= 187510000)


hmp <- read_delim('largedata/Hapmaps/hapmap_chr5.txt', delim="\t", na = c("", "NA", "NN"))

hmp5 <- hmp %>%
  filter(pos > 188350000 & pos <= 188650000)


hmp <- read_delim('largedata/Hapmaps/hapmap_chr6.txt', delim="\t", na = c("", "NA", "NN"))

hmp6 <- hmp %>%
  filter(pos > 171410000 & pos <= 171710000)



hmp_candidates <- rbind(hmp2, hmp3, hmp5, hmp6)

unique(hmp_data_candidates$GX_name)


load()

names <- read_csv("data/BG_MM_Gen_names.csv")

test <- read_csv("data/microbe_counts_per_genotype.csv")

test <- test %>%
  rename(MM_name=genotype) %>%
  left_join(dplyr::select(names, MM_name, GX_name))

genotypes_w_count <- unique(test$GX_name)

zz <- file("largedata/genotypes_w_count.txt", "wb")
writeBin(paste(genotypes_w_count, collapse="\n"), zz ) 
close(zz)


genotypes_test <- read_delim("largedata/genotypes_test.txt", delim = ' ', col_names = FALSE)
genotypes_test <- pull(genotypes_test, "X1")
genotypes_test[genotypes_test %in% names$GX_name]
genotypes_test[!(genotypes_test %in% names$GX_name)]
genotypes_test[!(genotypes_test %in% names$GX_name)]
sort(unique(names$GX_name))


#save(hmp_candidates, file="cache/GWAS2/hmp_candidates.rda")


hmp_candidates <- read_delim('largedata/hmp_v4.txt', delim="\t", na = c("", "NA", "NN"))



names <- read_csv("data/BG_MM_Gen_names.csv")

colnames(hmp_candidates)[1] <- "snp_id"

BG_panel <- colnames(hmp_candidates)[12: length(colnames(hmp_candidates))]

hmp_data <- hmp_candidates %>%
  filter(snp_id %in% c("S2-230389923", "S3-184434608", "S5-183779458", "S6-167618482")) %>%
  pivot_longer(BG_panel, names_to = "genotype", values_to = "haplotype") %>%
  dplyr::select(-starts_with("assembly"), -center, -protLSID, -assayLSID, -panelLSID, -QCcode) %>%
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


hmp_data <- hmp_data %>%
  left_join(NAs) %>%
  mutate(total_not_na = total_allele_obs-NA_freq) %>%
  mutate(allele = ifelse(allele_freq > total_not_na/2, "major allele", "minor allele")) %>%
  rename(rep_snp = snp_id) %>%
  left_join(unique(dplyr::select(top_12_MAPLs, rep_snp, tax_group, GWAS_nitrogen))) %>%
  rename(snp_id = rep_snp, genotype = MM_name)



hmp_data_candidates <- hmp_data

#save(hmp_data_candidates, file="cache/GWAS2/hmp_data_candidates.rda")



load("cache/GWAS2/hmp_data_candidates.rda")



load("data/yield_analysis/abundance_vs_phenotype.rda")


### get abundance data
load("cache/BLUP/blup_data.rda")

mean_blups <- blup_data %>%
  group_by(nitrogen, GX_name) %>%
  summarize(mean_blup = mean(blup_logrel))


### get correlation data
load("data/yield_analysis/corr_data.rda")



load("data/group_data.rda")


tgr <- "Sphingoaurantiacus"
snp <- "S2-230389923"
chr <- "2"

abundance_dat <- blup_data %>%
  left_join(dplyr::select(group_data, ASV, tax_group)) %>%
  filter(tax_group == tgr) %>%
  left_join(filter(hmp_data_candidates, rep_snp == snp))

test <- filter(hmp_data_candidates, chrom == chr)




p <- ggplot(mean_blups, aes(x=nitrogen, y=mean_blup, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~snp_id, nrow = 1) +
  stat_compare_means(aes(group = haplotype), label = "p.signif") +
  ylab("log(microbe relative abundance)") +
  #ylab("relative microbe abundance") +
  #ylab("ASV count") +
  scale_fill_manual(values = colors) +
  #geom_point(data=NAM_tax, aes(x=nitrogen, y=log_relab, fill=allele), pch = 21, size = 4, position = position_dodge(width = 1.5)) +
  #geom_text(data=NAM_tax[NAM_tax$allele=='major allele',], aes(label=genotype), size=3, hjust=1, fontface="bold", nudge_x = -0.45) +
  #geom_text(data=NAM_tax[NAM_tax$allele=='minor allele',], aes(label=genotype), size=3, hjust=0, fontface="bold", nudge_x = 0.45) +
  geom_point(data=NAM_tax, aes(x=nitrogen, y=log_relab, fill=allele), pch = 21, size = 4, position = position_dodge(width = 0.75)) +
  geom_text(data=NAM_tax[NAM_tax$allele=='major allele',], aes(label=genotype), size=3, hjust=1, fontface="bold", nudge_x = -0.25) +
  geom_text(data=NAM_tax[NAM_tax$allele=='minor allele',], aes(label=genotype), size=3, hjust=0, fontface="bold", nudge_x = 0.25) +
  
  #scale_color_manual(values = cols) +
  ggtitle(paste(taxgrp, "|", tax$GWAS_nitrogen[1])) +
  theme_bw() +
  theme(axis.title.x = element_blank())








