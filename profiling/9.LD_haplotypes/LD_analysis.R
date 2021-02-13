#LD plots

library("tidyverse")
library("genetics")
library("LDheatmap")
library("stringr")



### find top bins

#load("cache/top_taxa_bins_genes.rda")
#unique(top_bins$tax_group)

#hmp_bins <- top_bins %>%
#  ungroup() %>%
#  dplyr::select(chr, bin, sign_snps) %>%
#  unique()

#write_csv(hmp_bins, file = "cache/hmp_bins.csv")



### get trait data 
load("cache/top_taxa_bins_genes.rda")

trait_dat <- top_bins %>%
  ungroup() %>%
  dplyr::select(trait, tax_group) %>%
  unique()


### get haplotype data for top bins

haplotypes <- read.table("cache/top_snps_hmp.txt",head=T,com="",sep="\t",row=1,na.string="NN")
haplotypes <- haplotypes %>% 
  dplyr::select(-index) %>%
  rownames_to_column(var = "snp_id")
  



### get significant SNPs
load("largedata/combined_SNP_data.rda")
gwas_dat <- gwas_dat %>%
  filter(trait %in% trait_dat$trait) %>% ## top 22 microbes
  left_join(trait_dat) %>%
  mutate(snp_id = paste0("S", chr, "_", ps))


RhM <- gwas_dat %>%
  filter(tax_group == "Rhizobium mesosinicum")




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



#hapmap("Acinetobacter nosocomialis", "+N", 3, c(18736, 18741), draw=TRUE)


arch <- hapmap("Acinetobacter nosocomialis", "+N", 3, c(18741), draw=FALSE)

top_loci <- top_bins %>%
  ungroup() %>%
  dplyr::select(tax_group, nitrogen, chr, bin, nitrogen) %>%
  unique()


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


sort(unique())






sort(unique(gwas_dat$tax_group))


tax <- gwas_dat %>%
  filter(tax_group == "Acinetobacter nosocomialis") %>%
  mutate(log10p = -log10(p_wald))




### find most significant SNP for each GWAS peak to determine haplotypes

test <- top_bins

test$top_snp













