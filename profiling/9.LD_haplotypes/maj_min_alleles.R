

library("tidyverse")
library("ggpubr")
library("phyloseq")



load("data/top_snps_haplotypes.rda")



top10_taxa <- c("Acinetobacter nosocomialis", "Candidatus Udaeobacter copiosus",
  "f_A21b", "f_Comamonadaceae Unknown Genus", "Filimonas sp 2", "Ilumatobacter",
  "Massilia niabensis", "Niabella yanshanensis", "Rhizobium daejeonense", "Sphingobium herbicidovorans 1")


colors <- c("major allele"="#ecb602", "minor allele"="#a45ee5")



NAM_set <- c("B73", "B97", "CML103", "CML228", "CML277", "CML322",
  "CML333", "CML52", "CML69", "HP301", "IL14H", "KI11",
  "M162W", "M37W", "MS71", "NC350", "NC358", "OH43",
  "OH7B", "P39", "TZI8")



### get BLUP data +N and -N


load("data/group_data.rda")

traits <- group_data %>%
  dplyr::select(tax_group, trait) %>%
  filter(tax_group %in% top10_taxa)


blup_stdN <- read_csv("data/blup_stdN_150_tax_groups.csv")[, c("genotype", traits$trait)] %>%
  #filter(genotype %in% NAM_set) %>%
  pivot_longer(cols=starts_with("T"), names_to="trait", values_to="blup_stdN")

blup_lowN <- read_csv("data/blup_lowN_150_tax_groups.csv")[, c("genotype", traits$trait)] %>%
  #filter(genotype %in% NAM_set) %>%
  pivot_longer(cols=starts_with("T"), names_to="trait", values_to="blup_lowN")

blups <- blup_stdN %>%
  left_join(blup_lowN) %>% 
  left_join(traits) %>%
  dplyr::select(-trait) %>%
  pivot_longer(cols=starts_with("blup"), names_to="nitrogen", values_to="blup") %>%
  mutate(nitrogen = ifelse(nitrogen == "blup_stdN", "+N","-N"))



top_snps <- top_snps %>%
  left_join(blups)






### mark B73 and Nam parent to compare


taxgrp <- top10_taxa[1]
nitr = "+N"


tax <- top_snps %>%
  #filter(snp_id == "S8_119562405") %>%
  filter(tax_group %in% taxgrp) %>%
  filter(!genotype %in% NAM_set)


NAM_top_snps <- top_snps %>%
  filter(genotype %in% NAM_set)


NAM_tax <- NAM_top_snps %>%
  #filter(snp_id == "S8_119562405") %>%
  filter(tax_group %in% taxgrp) %>%
  mutate(col = ifelse(genotype == "B73", "green", "red"))


B73_tax <- NAM_top_snps %>%
  #filter(snp_id == "S8_119562405") %>%
  filter(tax_group %in% taxgrp) %>%
  filter(genotype == "B73")

parent_tax <- NAM_top_snps %>%
  #filter(snp_id == "S8_119562405") %>%
  filter(tax_group %in% taxgrp) %>%
  group_by(nitrogen) %>%
  filter(allele == "minor allele") %>%
  filter(count != 0) %>%
  filter(!(is.na(blup))) %>%
  filter(blup == min(blup))

marked <- rbind(B73_tax, parent_tax)

marked <- marked %>%
  filter(nitrogen == nitr) %>%
  mutate(col = ifelse(genotype == "B73", "#0f1314", "#d53e2f"))


marked_col <- marked %>%
  dplyr::select(genotype, col) %>%
  unique()

cols <- marked_col$col
names(cols) <- marked_col$genotype

p <- ggplot(tax, aes(x=nitrogen, y=blup, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~snp_id, nrow = 1) +
  stat_compare_means(aes(group = haplotype), label = "p.signif") +
  ylab("BLUP (microbe relative abundance)") +
  #ylab("relative microbe abundance") +
  #ylab("ASV count") +
  scale_fill_manual(values = colors) +
  geom_point(data=NAM_tax, aes(x=nitrogen, y=blup), pch = 16, size = 2, color = 'red', position = position_jitterdodge(jitter.width = 0.05)) +
  geom_point(data=marked, aes(x=nitrogen, y=blup, color=genotype), pch = 16, size = 4, position = position_jitterdodge(jitter.width = 0)) +
  scale_color_manual(values = cols) +
  ggtitle(paste(taxgrp, "|", tax$GWAS_nitrogen[1])) +
  theme_bw() +
  theme(axis.title.x = element_blank())

p


### show only +N



tax <- top_snps %>%
  #filter(snp_id == "S8_119562405") %>%
  filter(tax_group %in% taxgrp)
  #filter(!genotype %in% NAM_set)



p <- ggplot(tax, aes(x=nitrogen, y=blup, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~snp_id, nrow = 1) +
  stat_compare_means(aes(group = haplotype), label = "p.signif") +
  ylab("Log(relative microbe abundance)") +
  scale_y_continuous(limits = c(-4.7, -3.6)) +
  scale_fill_manual(values = colors) +
  #geom_point(data=filter(NAM_tax, nitrogen == "+N"), aes(x=nitrogen, y=blup), pch = 16, size = 2, color = 'red', position = position_jitterdodge(jitter.width = 0.05)) +
  #geom_point(data=filter(marked, nitrogen == "+N"), aes(x=nitrogen, y=blup, color=genotype), pch = 16, size = 4, position = position_jitterdodge(jitter.width = 0)) +
  scale_color_manual(values = cols) +
  #ggtitle(paste(taxgrp, "|", tax$GWAS_nitrogen[1])) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


p

?stat_compare_means


## hypothetical

hypdat <- marked %>%
  ungroup() %>%
  filter(nitrogen == "+N") %>%
  dplyr::select(genotype, blup, allele)

## add dummy data

hypdat <- rbind(hypdat, c("RIL1", -3.97, "major allele"))
hypdat <- rbind(hypdat, c("RIL1", -4.16, "major allele"))
hypdat <- rbind(hypdat, c("RIL2", -4.36, "minor allele"))
hypdat <- rbind(hypdat, c("RIL2", -4.16, "minor allele"))


hypdat$genotype <- factor(hypdat$genotype, levels = c("B73", "RIL1", "RIL2", "CML333"))
hypdat$blup <- as.numeric(hypdat$blup)

str(hypdat)


cols <- c("major allele"="#ecb602", "minor allele"="#a45ee5")

cols <- c("B73"="#ecb602", "RIL1"="#bbbbbb", "RIL2"="#bbbbbb", "CML333"="#a45ee5")
cols <- c("B73"="#0f1314", "RIL1"="#aaaaaa", "RIL2"="#aaaaaa", "CML333"="#d53e2f")



p <- ggplot(hypdat, aes(x=genotype, y=blup, color=genotype)) +
  geom_point(pch = 16, size = 4) +
  geom_point(pch = 1, size = 4, color="black") +
  ylab("BLUP (microbe relative abundance)") +
  scale_y_continuous(limits = c(-4.7, -3.6)) +
  scale_color_manual(values = cols) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  
p






