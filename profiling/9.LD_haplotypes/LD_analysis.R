#LD plots

library("tidyverse")
library("genetics")
library("LDheatmap")
library("stringr")
library("ggpubr")



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


test <- filter(top_bins, tax_group == "Acinetobacter nosocomialis")



### get haplotype data for top bins

haplotypes <- read.table("cache/top_snps_hmp_2.txt",head=T,com="",sep="\t",row=1,na.string="NN")
haplotypes <- haplotypes %>% 
  dplyr::select(-index) %>%
  rownames_to_column(var = "snp_id")
  





### get significant SNPs
load("largedata/combined_SNP_data.rda")
gwas_dat <- gwas_dat %>%
  filter(trait %in% trait_dat$trait) %>% ## top 22 microbes
  left_join(trait_dat) %>%
  mutate(snp_id = paste0("S", chr, "_", ps))


#RhM <- gwas_dat %>%
  #filter(tax_group == "Rhizobium mesosinicum")

ano <- gwas_dat %>%
  filter(tax_group == "Acinetobacter nosocomialis")

ano$log10p <- -log10(ano$p_wald)


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



hapmap("Acinetobacter nosocomialis", "+N", 3, c(18736, 18739, 18741), draw=TRUE)


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






## pivot longer and fix genotype nomenclature

genotype_names_all <- read_csv("data/BG_MM_Gen_names.csv")
gnames <- genotype_names_all %>%
  dplyr::select(MM_name, GX_name) %>%
  rename(genotype = MM_name)


top_snps <- top_snps %>%
  dplyr::select(-X1, -assembly., -protLSID, -assayLSID, -panelLSID, -QCcode, -center, -strand) %>%
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
  drop_na(genotype) 




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



genotype_names_hap <- sort(unique(top_snps$genotype))

## missing genotypes
genotype_names_hap_GX <- sort(unique(top_snps$GX_name))
genotype_names_hap_MM <- sort(unique(gnames$genotype))

## missing 16 genotypes
genotype_names_hap_MM[!(genotype_names_hap_MM %in% genotype_names_hap)]


## remove heterozygous and NA

string <- "AA"
string <- "TT"
string <- "AT"

str_detect(string, "TT|AA|GG|CC")


top_snps <- top_snps %>%
  mutate(haplotype = ifelse(str_detect(haplotype, "TT|AA|GG|CC"), haplotype, NA)) %>%
  drop_na(haplotype)


#save(top_snps, file="data/top_snps_haplotypes.rda")



top10_taxa <- c("Acinetobacter nosocomialis", "Candidatus Udaeobacter copiosus",
  "f_A21b", "f_Comamonadaceae Unknown Genus", "Filimonas sp 2", "Ilumatobacter",
  "Massilia niabensis", "Niabella yanshanensis", "Rhizobium daejeonense", "Sphingobium herbicidovorans 1")


colors <- c("AA"="#d2e2ef", "CC"="#fcd2c2", "GG"="#dfe8d5", "TT"="#ffeabf")



taxgrp <- top10_taxa[10]

tax <- top_snps %>%
  filter(tax_group %in% taxgrp)

p <- ggplot(tax, aes(x=nitrogen, y=count, fill=haplotype)) +
    #geom_boxplot() +
    geom_boxplot(outlier.shape = NA) +
    geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
    facet_wrap(~snp_id, nrow = 1) +
    stat_compare_means(aes(group = haplotype), label = "p.signif") +
    ylab("ASV count") +
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


