#LD plots

library("genetics")
library("LDheatmap")
library("stringr")
library("ggpubr")
library("tidyverse")



### find top bins

#load("cache/top_taxa_bins_genes_10.rda")
#unique(top_bins$tax_group)

#hmp_bins <- top_bins %>%
  #ungroup() %>%
  #dplyr::select(chr, bin, count) %>%
  #unique()

#write_csv(hmp_bins, file = "cache/hmp_bins.csv")



### get trait data 
load("cache/top_taxa_bins_genes.rda")

trait_dat <- top_bins %>%
  ungroup() %>%
  dplyr::select(trait, tax_group) %>%
  unique()


### get haplotype data for top bins

haplotypes <- read.table("largedata/top_snps_hmp_3.txt",head=T,com="",sep="\t",row=1,na.string="NN")
haplotypes <- haplotypes %>% 
  dplyr::select(-index) %>%
  rownames_to_column(var = "snp_id")
  



### get significant SNPs
load("largedata/combined_SNP_data.rda")
gwas_dat <- gwas_dat %>%
  filter(trait %in% trait_dat$trait) %>% ## top 22 microbes
  left_join(trait_dat) %>%
  mutate(snp_id = paste0("S", chr, "_", ps)) %>%
  mutate(log10p = -log10(p_wald))


#RhM <- gwas_dat %>%
  #filter(tax_group == "Rhizobium mesosinicum")

#ano <- gwas_dat %>%
  #filter(tax_group == "Acinetobacter nosocomialis")

#ano$log10p <- -log10(ano$p_wald)


### function to find SNPs to plot
### for given microbe, N treatment and range


hapmap <- function(microbe, N_treatment, chromosome, bins, draw=TRUE) {
  
  
  #microbe <- "Rhizobium mesosinicum"
  #N_treatment <- "+N"
  #chromosome <- 3
  #bins <- c(18741)
  
  
  LD_snps <- gwas_dat %>%
    filter(tax_group == microbe) %>%
    filter(nitrogen == N_treatment) %>%
    filter(chr == chromosome)
  
  #print(LD_snps$snp_id)
  
  
  LD_haplotypes <- haplotypes %>%
    filter(snp_id %in% LD_snps$snp_id)

  selected_haplotypes <- LD_haplotypes %>%
    filter(snp_id == 0)
  
  for (bin in bins) {
    start <- as.numeric(paste0(as.character(bin-1), "0000"))
    end <- as.numeric(paste0(as.character(bin), "0000"))
    add_haplotypes <- LD_haplotypes %>%
      filter(pos > start & pos <= end)
    selected_haplotypes <- rbind(selected_haplotypes, add_haplotypes)
  }
  
  if (draw==TRUE) {
    
    ## make filename
    binrange <- ifelse(length(bins) > 1, paste0(min(bins), "-", max(bins)), max(bins))
    ntreat <- ifelse(N_treatment == "+N", "stdN", "lowN")
    fname <- paste(str_replace(microbe," ","_"), ntreat, binrange, sep = "_")
    fpath <- paste0("figures/LD_plots/", fname, ".png")
    
    gene <- data.frame(t(selected_haplotypes[,12:ncol(selected_haplotypes)]))
    gty <- makeGenotypes(gene,sep="")
    
    rgb.palette <- colorRampPalette(rev(c("blue","orange" ,"red")), space = "rgb")
    
    png(file=fpath,res=300,width=1000,height=1000)
    
    myld <- LDheatmap(gty,genetic.distances=selected_haplotypes$pos,flip=TRUE, text=FALSE, color=rgb.palette(20), title=fname)
    
    dev.off()
    
  }
  
  selected_haplotypes <- dplyr::select(selected_haplotypes, snp_id) %>%
    left_join(LD_snps) %>%
    mutate(log10p = -log10(p_wald)) %>%
    left_join(selected_haplotypes)
  
  ##print(selected_haplotypes)
  
  return(selected_haplotypes)
  
}



hapmap("Acinetobacter nosocomialis", "+N", 3, c(18736, 18739, 18740, 18741), draw=TRUE)

hapmap("Candidatus Udaeobacter copiosus", "-N", 9, c(10355, 10356, 10358), draw=TRUE)

hapmap("Massilia niabensis", "-N", 3, c(16872, 16877), draw=TRUE)

hapmap("Niabella yanshanensis", "+N", 8, c(11943, 11954, 11957, 11960, 11969), draw=TRUE)

hapmap("f_Comamonadaceae Unknown Genus", "+N", 2, c(21743, 21744), draw=TRUE)
hapmap("f_Comamonadaceae Unknown Genus", "+N", 5, c(14749, 14750, 14753, 14754, 14755, 14756), draw=TRUE)

hapmap("f_A21b", "+N", 3, c(17437, 17438, 17441), draw=TRUE)
hapmap("f_A21b", "+N", 9, c(9224, 9225, 9229, 9230, 9253, 9255, 9274, 9275), draw=TRUE)

hapmap("Filimonas sp 2", "-N", 1, c(1277, 1278, 1280, 1282, 1285), draw=TRUE)
hapmap("Filimonas sp 2", "-N", 1, c(21325, 21326, 21327, 21331, 21332, 21333, 21334, 21335, 21338, 21340, 21344), draw=TRUE)
hapmap("Filimonas sp 2", "-N", 5, c(3757, 3769, 3771, 3772, 3773, 3774), draw=TRUE)

hapmap("Sphingobium herbicidovorans 1", "-N", 2, c(6979), draw=TRUE)
hapmap("Sphingobium herbicidovorans 1", "-N", 3, c(310, 311), draw=TRUE)
hapmap("Sphingobium herbicidovorans 1", "-N", 8, c(14527, 14529, 14530), draw=TRUE)
hapmap("Sphingobium herbicidovorans 1", "-N", 9, c(1163), draw=TRUE)

hapmap("Rhizobium daejeonense", "-N", 4, c(21403), draw=TRUE)
hapmap("Rhizobium daejeonense", "-N", 6, c(1818, 1834, 1839, 1843, 1844, 1845, 1846, 1850, 1853, 1854, 1863, 1864, 1867, 1868), draw=TRUE)

hapmap("Ilumatobacter", "-N", 2, c(10621, 10622, 10623, 10624, 10625, 10626, 10633, 10634, 10636), draw=TRUE)
hapmap("Ilumatobacter", "-N", 5, c(8326, 8327, 8328, 8329, 8333, 8343, 8344, 8345, 8346, 8363, 8371, 8372, 8378, 8379, 8380), draw=TRUE)


sort(unique(top_bins$tax_group))




arch <- hapmap("Acinetobacter nosocomialis", "+N", 3, c(18741), draw=FALSE)



top_loci <- top_bins %>%
  ungroup() %>%
  dplyr::select(tax_group, nitrogen, chr, bin, nitrogen) %>%
  unique()



## make dummy table to fill in
top_snps <- hapmap("Acinetobacter nosocomialis", "+N", 3, c(18741), draw=FALSE) %>%
  filter(log10p == max(log10p))

#nrow(top_loci)
for (loc in 1:nrow(top_loci)) {
  microbe <- top_loci[loc,]$tax_group
  N_treatment <- top_loci[loc,]$nitrogen
  chromosome <- top_loci[loc,]$chr
  bins <- c(top_loci[loc,]$bin)
  
  top_snp <- hapmap(microbe, N_treatment, chromosome, bins, draw=FALSE)
  
  if (nrow(top_snp) > 0) {
    top_snp <- filter(top_snp, log10p == max(log10p))
    top_snps <- rbind(top_snps, top_snp)
  }
  
}





## handpick loci from top 10 traits based on LD plot


trait <- gwas_dat %>%
  filter(tax_group == "Niabella yanshanensis") %>%
  filter(ps >= 119590000 & ps <= 119600000) %>%
  arrange(p_wald)
trait$log10p <- -log10(trait$p_wald)


## pivot longer and fix genotype nomenclature

genotype_names_all <- read_csv("data/BG_MM_Gen_names.csv")
gnames <- genotype_names_all %>%
  dplyr::select(MM_name, GX_name) %>%
  rename(genotype = MM_name)



### find top SNPs
handpicked_snps <- read_csv("cache/handpicked_snps.csv")


top_snps <- handpicked_snps %>%
  left_join(gwas_dat) %>%
  left_join(haplotypes) %>%
  #dplyr::select(-X1, -assembly., -protLSID, -assayLSID, -panelLSID, -QCcode, -center, -strand) %>%
  pivot_longer(
    cols = X33.16:Yu796,
    names_to = "GX_name",
    values_to = "haplotype"
  ) %>%
  mutate(GX_name = str_replace_all(GX_name, "[.]", "-")) %>%
  mutate(GX_name = str_replace_all(GX_name, "^X", "")) %>%
  mutate(GX_name = str_replace_all(GX_name, "MO1W", "Mo1W")) %>%
  mutate(GX_name = str_replace_all(GX_name, "Ms71", "MS71")) %>%
  mutate(GX_name = str_replace_all(GX_name, "OH7B", "Oh7B")) %>%
  left_join(gnames) %>%
  drop_na(genotype) %>%
  mutate(haplotype = ifelse(haplotype %in% c("TT", "AA", "GG", "CC"), haplotype, NA)) %>%
  drop_na(haplotype) %>%
  group_by(tax_group, nitrogen, snp_id, haplotype) %>% mutate(obs_per_allele = n()) %>%
  group_by(tax_group, nitrogen, snp_id) %>% mutate(total_obs = n()) %>%
  mutate(allele = ifelse(obs_per_allele >= total_obs/2, "major allele", "minor allele"))
 


### get microbe counts for stdN and lowN
load("cache/counts_per_genotype.rda")

count_data_stdN <- filter(count_data, nitrogen == "+N")
top_snps_stdN <- top_snps %>%
  rename(GWAS_nitrogen = nitrogen) %>%
  mutate(nitrogen = "+N") %>%
  left_join(count_data_stdN)

count_data_lowN <- filter(count_data, nitrogen == "-N")
top_snps_lowN <- top_snps %>%
  rename(GWAS_nitrogen = nitrogen) %>%
  mutate(nitrogen = "-N") %>%
  left_join(count_data_lowN)

top_snps <- rbind(top_snps_stdN, top_snps_lowN)



#genotype_names_hap <- sort(unique(top_snps$genotype))

## missing genotypes
#genotype_names_hap_GX <- sort(unique(top_snps$GX_name))
#genotype_names_hap_MM <- sort(unique(gnames$genotype))

## missing 16 genotypes
#genotype_names_hap_MM[!(genotype_names_hap_MM %in% genotype_names_hap)]


top_snps$logrel <- log(top_snps$relab + 0.001)


save(top_snps, file="data/top_snps_haplotypes.rda")

#load("data/top_snps_haplotypes.rda")


top10_taxa <- c("Acinetobacter nosocomialis", "Candidatus Udaeobacter copiosus",
  "f_A21b", "f_Comamonadaceae Unknown Genus", "Filimonas sp 2", "Ilumatobacter",
  "Massilia niabensis", "Niabella yanshanensis", "Rhizobium daejeonense", "Sphingobium herbicidovorans 1")


#colors <- c("AA"="#d2e2ef", "CC"="#fcd2c2", "GG"="#dfe8d5", "TT"="#ffeabf")


colors <- c("major allele"="#ecb602", "minor allele"="#a45ee5")


taxgrp <- top10_taxa[1]

tax <- top_snps %>%
  filter(tax_group %in% taxgrp)


p <- ggplot(tax, aes(x=nitrogen, y=logrel, fill=allele)) +
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

p

#dev.off()

#compare_means(count ~ haplotype, data = tax, 
  #group.by = c("nitrogen","snp_id"))



taxgrp <- top10_taxa[10]

snp_pvalues <- top_snps %>%
  dplyr::select(tax_group, GWAS_nitrogen, snp_id, log10p) %>%
  unique() %>%
  filter(tax_group == taxgrp)



### select NAM parents


NAM_set <- c("B73", "B97", "CML103", "CML228", "CML277", "CML322",
  "CML333", "CML52", "CML69", "HP301", "IL14H", "KI11",
  "M162W", "M37W", "MS71", "NC350", "NC358", "OH43",
  "OH7B", "P39", "TZI8")




NAM_top_snps <- top_snps %>%
  filter(genotype %in% NAM_set)


taxgrp <- top10_taxa[1]

tax <- NAM_top_snps %>%
  filter(tax_group %in% taxgrp)


p <- ggplot(tax, aes(x=genotype, y=count, fill=haplotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  #geom_boxplot(outlier.shape = NA) +
  #geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~nitrogen + snp_id, nrow = 2) +
  ylab("ASV count") +
  scale_fill_manual(values = colors) +
  ggtitle(paste(taxgrp, "|", tax$GWAS_nitrogen[1])) +
  theme_bw()

p












### mark B73 and Nam parent to compare


taxgrp <- top10_taxa[9]
nitr = "-N"


tax <- top_snps %>%
  #filter(snp_id == "S8_119562405") %>%
  filter(tax_group %in% taxgrp) %>%
  filter(!genotype %in% NAM_set)

NAM_tax <- NAM_top_snps %>%
  filter(snp_id == "S8_119562405") %>%
  filter(tax_group %in% taxgrp) %>%
  mutate(col = ifelse(genotype == "B73", "green", "red"))


B73_tax <- NAM_top_snps %>%
  filter(snp_id == "S8_119562405") %>%
  filter(tax_group %in% taxgrp) %>%
  filter(genotype == "B73")

parent_tax <- NAM_top_snps %>%
  filter(snp_id == "S8_119562405") %>%
  filter(tax_group %in% taxgrp) %>%
  group_by(nitrogen) %>%
  filter(allele == "minor allele") %>%
  filter(count != 0) %>%
  filter(logrel == min(logrel))

marked <- rbind(B73_tax, parent_tax)

marked <- marked %>%
  filter(nitrogen == nitr) %>%
  mutate(col = ifelse(genotype == "B73", "blue", "green"))


marked_col <- marked %>%
  dplyr::select(genotype, col) %>%
  unique()

cols <- marked_col$col
names(cols) <- marked_col$genotype

p <- ggplot(tax, aes(x=nitrogen, y=logrel, fill=allele)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  facet_wrap(~snp_id, nrow = 1) +
  stat_compare_means(aes(group = haplotype), label = "p.signif") +
  ylab("log(relative microbe abundance)") +
  #ylab("relative microbe abundance") +
  #ylab("ASV count") +
  scale_fill_manual(values = colors) +
  geom_point(data=NAM_tax, aes(x=nitrogen, y=logrel), pch = 16, size = 2, color = 'red', position = position_jitterdodge(jitter.width = 0.1)) +
   geom_point(data=marked, aes(x=nitrogen, y=logrel, color=genotype), pch = 16, size = 4, position = position_jitterdodge(jitter.width = 0.1)) +
  scale_color_manual(values = cols)
  ggtitle(paste(taxgrp, "|", tax$GWAS_nitrogen[1])) +
  theme_bw() +
  theme(axis.title.x = element_blank())

p







