# yield and fitness comparisons
# compare Massilia abundance to plant performance traits
# for +N and -N for major and minor allele

library("tidyverse")
library("ggpubr")



## get haplotype data for each genotype for top interesting SNPs

load("data/top_snps_haplotypes.rda")



# Get mean fitness traits for each genotype


CC_raw <- read_csv("data/yield_analysis/Table_S1.canopy_coverage_by_dates.csv")

CC <- CC_raw %>%
  #filter(date == "Aug12") %>%
  rename(GX_name = Genotype) %>%
  rename(nitrogen = Treatment) %>%
  mutate(date = paste0("CC_", date)) %>%
  group_by(GX_name, nitrogen,date) %>%
  dplyr::summarize(mean_CC = mean(Canopy_Coverage, na.rm=TRUE)) %>%
  pivot_wider(names_from=date, values_from=mean_CC) %>%
  mutate(nitrogen = ifelse(nitrogen == "Nitrogen", "+N", "-N"))
  

VI_raw <- read_csv("data/yield_analysis/Table_S2.VIs_by_dates.csv")

VI <- VI_raw %>%
  rename(GX_name = Genotype) %>%
  #filter(date == "Aug12") %>%
  rename(nitrogen = Treatment) %>%
  mutate(date = paste0("ExG_", date)) %>%
  group_by(GX_name, nitrogen,date) %>%
  dplyr::summarize(mean_ExG = mean(ExG, na.rm=TRUE)) %>%
  pivot_wider(names_from=date, values_from=mean_ExG) %>%
  mutate(nitrogen = ifelse(nitrogen == "Nitrogen", "+N", "-N"))






# combine data
plant_fitness_data <- top_snps %>%
  left_join(CC) %>%
  left_join(VI)





### look at top 10 taxa and associated SNPs

top10_taxa <- c("Acinetobacter nosocomialis", "Candidatus Udaeobacter copiosus",
                "f_A21b", "f_Comamonadaceae Unknown Genus", "Filimonas sp 2", "Ilumatobacter",
                "Massilia niabensis", "Niabella yanshanensis", "Rhizobium daejeonense", "Sphingobium herbicidovorans 1")

taxgrp <- top10_taxa[7]

tax <- plant_fitness_data %>%
  filter(tax_group %in% taxgrp)



colors <- c("major allele"="#ecb602", "minor allele"="#a45ee5")




p <- ggplot(tax, aes(x=logrel, y=ExG_Sept1, color=allele)) +
  facet_wrap(~nitrogen, nrow = 1) +
  #scale_color_manual(values = colors) +
  geom_point()

p


ExG <- tax %>%
  ungroup() %>%
  dplyr::select(genotype, nitrogen, allele, logrel, starts_with("ExG_")) %>%
  pivot_longer(cols = starts_with("ExG_"), names_to = "date", values_to = "ExG")


ExG$date <- factor(ExG$date, levels = c("ExG_July6", "ExG_Aug12", "ExG_Aug14", "ExG_Aug16", "ExG_Aug20", "ExG_Aug22", "ExG_Aug23", "ExG_Aug26", "ExG_Aug30", "ExG_Sept1", "ExG_Sept3", "ExG_Sept5" ))


p_ExG <- ggplot(ExG, aes(x=date, y=ExG, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~nitrogen, nrow = 2) +
  stat_compare_means(aes(group = allele), label = "p.signif") +
  ylab("ExG") +
  #ylab("relative microbe abundance") +
  #ylab("ASV count") +
  scale_fill_manual(values = colors) +
  theme_bw()

p_ExG




CC <- tax %>%
  ungroup() %>%
  dplyr::select(genotype, nitrogen, allele, logrel, starts_with("CC_")) %>%
  pivot_longer(cols = starts_with("CC_"), names_to = "date", values_to = "CC")


CC$date <- factor(CC$date, levels = c("CC_July6", "CC_Aug12", "CC_Aug14", "CC_Aug16", "CC_Aug20", "CC_Aug22", "CC_Aug23", "CC_Aug26", "CC_Aug30", "CC_Sept1", "CC_Sept3", "CC_Sept5" ))

p_CC <- ggplot(CC, aes(x=date, y=CC, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~nitrogen, nrow = 2) +
  stat_compare_means(aes(group = allele), label = "p.signif") +
  ylab("Canopy Cover") +
  #ylab("relative microbe abundance") +
  #ylab("ASV count") +
  scale_fill_manual(values = colors) +
  theme_bw()

p_CC





### plot microbe abundance and canopy cover




MNab <- ggplot(tax, aes(x=nitrogen, y=logrel, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~snp_id, nrow = 1) +
  stat_compare_means(aes(group = haplotype), label = "p.signif") +
  ylab("log(relative microbe abundance)") +
  #ylab("relative microbe abundance") +
  #ylab("ASV count") +
  scale_fill_manual(values = colors) +
  ggtitle(paste(taxgrp, "|", tax$GWAS_nitrogen[1])) +
  theme_bw() +
  theme(axis.title.x = element_blank())

MNab



MNcc <- ggplot(tax, aes(x=nitrogen, y=CC_Sept5, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~snp_id, nrow = 1) +
  stat_compare_means(aes(group = haplotype), label = "p.signif") +
  ylab("Canopy Cover (Sep 5)") +
  scale_fill_manual(values = colors) +
  ggtitle(paste(taxgrp, "|", tax$GWAS_nitrogen[1])) +
  theme_bw() +
  theme(axis.title.x = element_blank())

MNcc




