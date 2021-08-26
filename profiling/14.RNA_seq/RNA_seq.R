## gene_expression


library("tidyverse")
library("ggpubr")



fpkm <- read_csv("data/RNA_seq/FPKM_Old_new_four_lines2.csv")
fpkm <- rename(fpkm, gene=X1)

top_12_MAPLs <- read_csv("data/top_12_MAPLs.csv")


mapl_genes <- fpkm %>%
  filter(gene %in% top_12_MAPLs$gene) %>%
  ##dplyr::select(gene, contains("Root")) %>%
  pivot_longer(-gene, names_to = "Sample", values_to = "fpkm") %>%
  separate(Sample, c("nitrogen", "genotype", "tissue", "rep"), sep = "_") %>%
  mutate(nitrogen = ifelse(nitrogen == "HN", "+N", "-N"))
  #mutate(genotype, ifelse(genotype == "W153", "W153R", genotype))
mapl_genes[mapl_genes$genotype == "W153", ]$genotype <- "W153R" # fix name



#sort(unique(top_12_MAPLs$tax_group))

### add major/minor allele info for 4 genotypes

load("cache/maj_min.rda")


maj_min_set <- maj_min %>%
  ungroup() %>%
  dplyr::select(snp_id, genotype, nitrogen, allele, tax_group) %>%
  unique()

plot_dat <- mapl_genes %>%
  left_join(unique(dplyr::select(top_12_MAPLs, gene, snp_id))) %>%
  left_join(maj_min_set) 

sort(unique(plot_dat$tax_group))



work_on <- unique(dplyr::select(plot_dat, tax_group, nitrogen))[11, ]
taxgrp <- as.character(work_on[1, 1])
nitr <- as.character(work_on[1, 2])
taxgrp

taxgrp <- "Rhizobium daejeonense"
nitr <- "-N"

genes <- filter(plot_dat, tax_group == taxgrp)$gene
genes


unique(plot_dat$genotype)

#min_allele <- c("SD") # filimonas
min_allele <- c("B73") # rhizobium daejeonense

MAPL <- plot_dat %>%
  filter(gene %in% genes) %>%
  mutate(allele = ifelse(genotype %in% min_allele, "minor allele", "major allele")) %>%
  #left_join(maj_min_set) 
  #ungroup() %>%
  #unique() %>%
  group_by(gene, nitrogen, tissue, allele) %>%
  summarize(mean_fpkm = mean(fpkm, na.rm=TRUE), n = n(), 
            sd = sd(fpkm, na.rm=TRUE), se = sd/sqrt(n)) #%>%
  #mutate(N_genotype = paste(genotype, nitrogen))



#MAPL$N_genotype <- factor(MAPL$N_genotype, levels= c("B73 -N", "K55 -N", "W153R -N","SD -N", "B73 +N", "K55 +N", "W153R +N","SD +N")) 


colors <- c("major allele"="#ecb602", "minor allele"="#a45ee5")


#colors <- c("#ecb602", "#ecb602", "#ecb602", "#a45ee5","#ecb602", "#ecb602", "#ecb602", "#a45ee5")

ggplot(MAPL, aes(x=nitrogen, y=mean_fpkm, fill=allele)) +
  geom_bar(stat = "identity", position=position_dodge(), color="black") +
  #facet_wrap(~tissue + gene, nrow = 2, scales = "free_y") +
  facet_grid(tissue~gene) +
  #stat_compare_means(aes(group = allele), label = "p.signif") +
  geom_errorbar(aes(ymin=mean_fpkm-se, ymax=mean_fpkm+se), width=.2, position=position_dodge(.9)) +
  ggtitle(paste(taxgrp, "|", nitr)) +
  ylab("gene expression (FPKM)") +
  theme_bw() +
  scale_fill_manual( values = colors) +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust= 0.5))#+
  #theme(axis.title.x = element_blank())





#correlation gene expression with abundance


