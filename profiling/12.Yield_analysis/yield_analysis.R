# yield and fitness comparisons
# compare Massilia abundance to plant performance traits
# for +N and -N for major and minor allele

library("tidyverse")
library("ggpubr")
library("lme4")
library("g3tools")



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




### Yield component traits




names <- read_csv("data/BG_MM_Gen_names.csv")


cob_length <- read_csv("data/yield_analysis/pheno2019_cob_length.csv")
cob_width <- read_csv("data/yield_analysis/pheno2019_cob_width.csv")
cob_weight <- read_csv("data/yield_analysis/pheno2019_cob_weight.csv")


group_by(grp) %>% summarise(across(everything(), list(mean)))

yield_components <- cob_length %>%
  left_join(cob_width) %>%
  left_join(cob_weight) %>%
  dplyr::rename(GX_name = Genotype) %>%
  left_join(dplyr::select(names, MM_name, GX_name)) %>%
  filter(GX_name != "Check") %>%
  rename(genotype = MM_name) %>%
  rename(cob_length = Average.cob.length) %>%
  rename(cob_width = Average.cob.width) %>%
  rename(cob_weight = Average.cob.weight) %>%
  mutate(block = ifelse(row < 2999, "N", "S")) %>%
  dplyr::select(genotype, nitrogen, block, sb, spb, row, cob_width, cob_length, cob_weight)
  


#save(yield_components, file="data/yield_analysis/yield_components.rda")



## calculate BLUPs of yield components



traits <- c("cob_width", "cob_length", "cob_weight")

### BLUP std N

if(!(exists("blup_stdN"))){
  blup_stdN <- data.frame(genotype = sort(unique(yield_components$genotype)))
}

for (trait in traits){

  df <- yield_components %>%
    filter(nitrogen == "+N") %>%
    select_("genotype", "block", "sb", "spb", "row", trait)
  
  colnames(df)[ncol(df)] <- "value"
  
  get_BLUP(data = df, model = value ~ (1 | genotype) + (1 | block) + (1 | sb) + (1 | spb), which.factor = "genotype",
           outfile = "cache/BLUP/tmp_blup.csv")
  
  blup <- read.csv("cache/BLUP/tmp_blup.csv")
  colnames(blup)[1] <- "genotype"
  colnames(blup)[2] <- trait
  
  blup_stdN <- left_join(blup_stdN, blup)
  
}



if(!(exists("blup_lowN"))){
  blup_lowN <- data.frame(genotype = sort(unique(yield_components$genotype)))
}

for (trait in traits){
  
  df <- yield_components %>%
    filter(nitrogen == "-N") %>%
    select_("genotype", "block", "sb", "spb", "row", trait)
  
  colnames(df)[ncol(df)] <- "value"
  
  get_BLUP(data = df, model = value ~ (1 | genotype) + (1 | block) + (1 | sb) + (1 | spb), which.factor = "genotype",
           outfile = "cache/BLUP/tmp_blup.csv")
  
  blup <- read.csv("cache/BLUP/tmp_blup.csv")
  colnames(blup)[1] <- "genotype"
  colnames(blup)[2] <- trait
  
  blup_lowN <- left_join(blup_lowN, blup)
  
}

#combine BLUP data

blup_stdN$nitrogen <- "+N"
blup_lowN$nitrogen <- "-N"


yield_components_BLUPs <- rbind(blup_stdN, blup_lowN) %>%
  select(genotype, nitrogen, cob_width, cob_length, cob_weight)



#save(yield_components_BLUPs, file="data/yield_analysis/yield_components_BLUPs.rda")


load("data/yield_analysis/yield_components_BLUPs.rda")

yield_components_BLUPs


# combine yield data for each genotype

names <- read_csv("data/BG_MM_Gen_names.csv")
names <- names %>%
  rename(genotype=BG_original) %>%
  dplyr::select(genotype, GX_name)



yield_per_genotype <- yield_components_BLUPs %>%
    left_join(names) %>%
    left_join(CC) %>%
    left_join(VI) %>%
    dplyr::select(genotype, GX_name, nitrogen, cob_length, cob_width, cob_weight, CC_Aug12, ExG_Aug12)
  
  
save(yield_per_genotype, file="data/yield_analysis/yield_per_genotype.rda")



top_12_MAPLs <- read_csv("data/top_12_MAPLs.csv")


# combine data
plant_fitness_data <- top_snps %>%
  left_join(yield_components_BLUPs) %>%
  left_join(CC) %>%
  left_join(VI)



#save(plant_fitness_data, file="data/yield_analysis/plant_fitness_data.rda")


load("data/yield_analysis/plant_fitness_data.rda")


unique(plant_fitness_data$tax_group)


### look at top 10 taxa and associated SNPs

top10_taxa <- c("Acinetobacter nosocomialis", "Candidatus Udaeobacter copiosus",
                "f_A21b", "f_Comamonadaceae Unknown Genus", "Filimonas sp 2", "Ilumatobacter",
                "Massilia niabensis", "Niabella yanshanensis", "Rhizobium daejeonense", "Sphingobium herbicidovorans 1")

taxgrp <- top10_taxa[7]

tax <- plant_fitness_data %>%
  filter(tax_group %in% taxgrp)


colors <- c("major allele"="#ecb602", "minor allele"="#a45ee5")

unique(plant_fitness_data$genotype)


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



MNcc <- ggplot(tax, aes(x=nitrogen, y=CC_Aug14, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~snp_id, nrow = 1) +
  stat_compare_means(aes(group = haplotype), label = "p.signif") +
  ylab("Canopy Cover (Aug 14)") +
  scale_fill_manual(values = colors) +
  ggtitle(paste(taxgrp, "|", tax$GWAS_nitrogen[1])) +
  theme_bw() +
  theme(axis.title.x = element_blank())

MNcc





MNcwt <- ggplot(tax, aes(x=nitrogen, y=cob_weight, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~snp_id, nrow = 1) +
  stat_compare_means(aes(group = haplotype), label = "p.signif") +
  ylab("Cob Weight") +
  scale_fill_manual(values = colors) +
  ggtitle(paste(taxgrp, "|", tax$GWAS_nitrogen[1])) +
  theme_bw() +
  theme(axis.title.x = element_blank())

MNcwt




MNcwd <- ggplot(tax, aes(x=nitrogen, y=cob_width, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~snp_id, nrow = 1) +
  stat_compare_means(aes(group = haplotype), label = "p.signif") +
  ylab("Cob Width") +
  scale_fill_manual(values = colors) +
  ggtitle(paste(taxgrp, "|", tax$GWAS_nitrogen[1])) +
  theme_bw() +
  theme(axis.title.x = element_blank())

MNcwd





MNlen <- ggplot(tax, aes(x=nitrogen, y=cob_length, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~snp_id, nrow = 1) +
  stat_compare_means(aes(group = haplotype), label = "p.signif") +
  ylab("Cob Length") +
  scale_fill_manual(values = colors) +
  ggtitle(paste(taxgrp, "|", tax$GWAS_nitrogen[1])) +
  theme_bw() +
  theme(axis.title.x = element_blank())

MNlen

unique(top10_data$snp_id)


top10_data <- plant_fitness_data %>%
  filter(tax_group %in% top10_taxa) %>%
  select(tax_group, snp_id, genotype, nitrogen, GWAS_nitrogen, allele, cob_length, cob_width, cob_weight, tax_group, logrel, CC_Aug12) %>%
  filter(snp_id %in% c("S3_187388952", "S9_103545234", "S3_168717026", "S8_119562405")) %>%
  pivot_longer(cols = cob_length:CC_Aug12, names_to = "trait", values_to = "value")

top10_data$trait <- factor(top10_data$trait, levels = c("logrel", "CC_Aug12", "cob_length", "cob_width", "cob_weight"))

unique(top10_data$trait)

all12 <- ggplot(filter(top10_data, tax_group %in% c("Massilia niabensis", "Niabella yanshanensis")), aes(x=nitrogen, y=value, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  #facet_wrap(snp_id~trait, nrow = 6) +
  facet_wrap(snp_id~trait, nrow = 2, scales = "free") +
  stat_compare_means(aes(group = allele), label = "p.signif") +
  #ylab(trait) +
  scale_fill_manual(values = colors) +
  #ggtitle(paste(taxgrp, "|", tax$GWAS_nitrogen[1])) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

all12


unique(top10_data$tax_group)



all23 <- ggplot(filter(top10_data, tax_group %in% c("Acinetobacter nosocomialis", "Candidatus Udaeobacter copiosus")), aes(x=nitrogen, y=value, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  #facet_wrap(snp_id~trait, nrow = 6) +
  facet_wrap(snp_id~trait, nrow = 2, scales = "free") +
  stat_compare_means(aes(group = allele), label = "p.signif") +
  #ylab(trait) +
  scale_fill_manual(values = colors) +
  #ggtitle(paste(taxgrp, "|", tax$GWAS_nitrogen[1])) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

all23


#plant_fitness_data <- plant_fitness_data %>%
  #rename(asv_count = count) %>%
  #select(-X1, -assembly., -center, -protLSID, -assayLSID, -panelLSID, -QCcode)


#save(plant_fitness_data, file="data/yield_analysis/plant_fitness_data.rda")

