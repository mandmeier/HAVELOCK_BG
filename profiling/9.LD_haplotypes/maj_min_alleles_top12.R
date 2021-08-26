

library("tidyverse")
library("ggpubr")
library("phyloseq")
library("ggrepel")


unique(top_snps$tax_group)


load("data/top_snps_haplotypes.rda")

load("data/group_data.rda")

top_12_MAPLs <- read_csv("data/top_12_MAPLs.csv")

top12_taxa <- sort(unique(top_12_MAPLs$tax_group))


MAPLs <- unique(top_12_MAPLs$rep_snp)[2: length(unique(top_12_MAPLs$rep_snp))]



hmp_top <- filter(hmp, alleles == "dummy")
hmp <- read_delim('largedata/Hapmaps/hapmap_chr7.txt', delim="\t", na = c("", "NA", "NN"))
colnames(hmp)[1] <- "snp_id"
top <- filter(hmp, snp_id %in% MAPLs)
hmp_top <- rbind(hmp_top, top)


#save(hmp_top, file = "cache/hmp_top_12_MAPLs.rda")

load("cache/hmp_top_12_MAPLs.rda")


BG_panel <- colnames(hmp_top)[12: length(colnames(hmp_top))]

names <- read_csv("data/BG_MM_Gen_names.csv")





BG_panel <- colnames(hmp_top)[12: length(colnames(hmp_top))]


hmp_data <- hmp_top %>%
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



traits <- group_data %>%
  dplyr::select(tax_group, trait) %>%
  filter(tax_group %in% top12_taxa)



load("data/data_summary_150_traits.rda")




### get microbe counts for stdN and lowN
load("cache/counts_per_genotype.rda")

count_data_stdN <- filter(count_data, nitrogen == "+N")
top_snps_stdN <- hmp_data %>%
  #rename(GWAS_nitrogen = nitrogen) %>%
  mutate(nitrogen = "+N") %>%
  left_join(count_data_stdN)

count_data_lowN <- filter(count_data, nitrogen == "-N")
top_snps_lowN <- hmp_data %>%
  #rename(GWAS_nitrogen = nitrogen) %>%
  mutate(nitrogen = "-N") %>%
  left_join(count_data_lowN)

top_snps <- rbind(top_snps_stdN, top_snps_lowN) %>%
            drop_na()



#top_snps <- hmp_data %>%
 # left_join(blups) %>%
  #drop_na()
 






### mark B73 and Nam parent to compare


work_on <- unique(dplyr::select(top_12_MAPLs, tax_group, GWAS_nitrogen))[12, ]
taxgrp <- work_on[1, 1]
nitr <- work_on[1, 2]


tax <- top_snps %>%
  #filter(snp_id == "S8_119562405") %>%
  filter(tax_group %in% taxgrp) #%>%
  #filter(!genotype %in% NAM_set)




NAM_set <- c("B73", "SD40", "K55", "W153R") # RNA seq lines


maj_min <- top_snps %>%
  filter(genotype %in% NAM_set)

save(maj_min, file="cache/maj_min.rda")


NAM_top_snps <- top_snps %>%
  filter(genotype %in% NAM_set)


NAM_tax <- NAM_top_snps %>%
  #filter(snp_id == "S8_119562405") %>%
  filter(tax_group %in% taxgrp) %>%
  mutate(col = ifelse(genotype == "B73", "green", "red"))


colors <- c("major allele"="#ecb602", "minor allele"="#a45ee5")
tax$group <- ifelse( tax$allele == "major allele", "major allele", "minor allele")
tax$group <- ifelse( tax$genotype == "B73", "B73", tax$group)
tax$group <- ifelse( tax$genotype == "SD40", "SD40", tax$group)
tax$group <- ifelse( tax$genotype == "K55", "K55", tax$group)
tax$group <- ifelse( tax$genotype == "W153R", "W153R", tax$group)



p <- ggplot(tax, aes(x=nitrogen, y=log_relab, fill=allele)) +
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

p







test <- tax %>%
  filter(genotype %in% NAM_set)


