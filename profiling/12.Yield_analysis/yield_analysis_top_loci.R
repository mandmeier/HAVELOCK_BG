# yeld_top_traits

library("tidyverse")
library("ggpubr")

## 1) get top loci in genome overview

load("cache/top_loci_traits_all_snps.rda")

top_loci_traits_all_snps$snp_id <- paste0("S", top_loci_traits_all_snps$chr, "_", top_loci_traits_all_snps$ps)

## test SNP "S3_168717026" massilia niabensis
## correlate microbe abundance vs. canopy cover



microbe_counts <- read_csv("data/microbe_counts_per_genotype.csv")
load("data/yield_analysis/yield_per_genotype.rda")
microbe_counts$nitrogen <- ifelse(microbe_counts$nitrogen == "stdN", "+N", "-N")


microbe_counts_and_yield <- microbe_counts %>%
  left_join(dplyr::select(yield_per_genotype, -GX_name))

#save(microbe_counts_and_yield, file="data/yield_analysis/microbe_counts_and_yield.rda")

load("data/yield_analysis/microbe_counts_and_yield.rda")

### determine major vs minor allele genotypes at S3_168717026

# subset huge hapmap data for interesting SNPs

unique(top_loci_traits_all_snps$psabs)
unique(top_loci_traits_all_snps$tax_group)
length(unique(top_loci_traits_all_snps$snp_id))

lociXtax <- length(unique(paste(top_loci_traits_all_snps$psabs, top_loci_traits_all_snps$tax_group)))



top_loci_snps <- unique(top_loci_traits_all_snps$snp_id) # 357 loci, 127 tax groups, 907 combinations

snp_12699 <- data.frame("snp_id"= unique(top_loci_traits_all_snps$rs))


write_csv(snp_12699, file="cache/snp_12699.csv")




hmp_top <- filter(hmp, alleles == "dummy")
#file <- paste0('largedata/Hapmaps/hapmap_chr',9,'.txt')
#hmp <- read_delim(file, delim="\t", na = c("", "NA", "NN"))
for(i in c(1:10)){
  file <- paste0('largedata/Hapmaps/hapmap_chr',i,'.txt')
  hmp <- read_delim(file, delim="\t", na = c("", "NA", "NN"))
  colnames(hmp)[1] <- "snp_id"
  top <- filter(hmp, snp_id %in% top_loci_snps)
  hmp_top <- rbind(hmp_top, top)
}

hmp_357_loci <- hmp_top
#save(hmp_357_loci, file="largedata/Hapmaps/hmp_357_loci.rda")

load("largedata/Hapmaps/hmp_357_loci.rda")


### calculate NAs and MAF,  filter SNPs with MAF < 5

names <- read_csv("data/BG_MM_Gen_names.csv")
BG_panel <- colnames(hmp_357_loci)[12: length(colnames(hmp_357_loci))]


hmp_data <- hmp_357_loci %>%
  pivot_longer(BG_panel, names_to = "genotype", values_to = "haplotype") %>%
  dplyr::select(-starts_with("assembly"), -center, -protLSID, -assayLSID, -panelLSID, -QCcode) %>%
  rename(GX_name = genotype) %>%
  left_join(dplyr::select(names, GX_name, MM_name)) %>%
  mutate(haplotype_NN = ifelse(haplotype %in% c("AA", "CC", "GG", "TT"), haplotype, "NN")) %>% # anything not homozygpus is NN
  group_by(snp_id) %>%
  add_count(name = "total_allele_obs") %>%
  group_by(snp_id, haplotype_NN) %>%
  add_tally(name = "allele_freq")


NAs <- hmp_data %>%
  dplyr::select(snp_id, GX_name, haplotype_NN, allele_freq) %>%
  filter(haplotype_NN == "NN") %>%
  ungroup() %>%
  dplyr::select(-GX_name, -haplotype_NN) %>%
  rename(NA_freq = allele_freq) %>%
  unique()

MAF <- hmp_data %>%
  dplyr::select(snp_id, GX_name, haplotype_NN, allele_freq) %>%
  ungroup() %>%
  dplyr::select(snp_id, allele_freq, haplotype_NN) %>%
  unique() %>%
  filter(haplotype_NN %in% c("AA", "CC", "GG", "TT")) %>%
  group_by(snp_id) %>%
  mutate(MAF = min(allele_freq)) %>%
  dplyr::select(snp_id, MAF) %>%
  unique()
  

hmp_data <- hmp_data %>%
  left_join(NAs) %>%
  left_join(MAF) %>%
  mutate(allele = "major allele") %>%
  mutate(allele = ifelse(allele_freq == MAF, "minor allele", allele)) %>%
  mutate(allele = ifelse(allele_freq == NA_freq, NA, allele)) %>%
  dplyr::select(-GX_name) %>%
  rename(genotype=MM_name) %>%
  filter(!(is.na(genotype))) %>%
  filter(!(is.na(allele))) %>%
  filter(MAF >= 5)
  
#save(hmp_data, file="largedata/Hapmaps/hmp_data.rda")

load("largedata/Hapmaps/hmp_data.rda")

# retain loci for which we have hapmap data with MAF >= 5
# and select highest SNP in each locus

## filter highest snp for which hapmap data is available
top_loci_traits_top_snps <- top_loci_traits_all_snps %>%
  filter(snp_id %in% unique(hmp_data$snp_id)) %>%
  group_by(nitrogen, psabs, tax_group) %>%
  filter(log10p == max(log10p))

unique(top_loci_traits_top_snps$psabs)
unique(top_loci_traits_top_snps$tax_group)
unique(top_loci_traits_top_snps$snp_id)

nrow(top_loci_traits_top_snps)

## hmp data available for 716 representative snps in 355 loci, 123 tax groups, 804 combinations

### save. This is the loci data above the 23 threshold for which hapmap data is available
## and which have minor allele frequency >= 5

#save(top_loci_traits_top_snps, file="cache/Hapmaps/top_loci_traits_top_snps.rda")

load("cache/Hapmaps/top_loci_traits_top_snps.rda")


## 2) do yield analysis for each of 921 snps with corresponding microbe

nrow(top_loci_traits_top_snps)


yield_ttest <- filter(top_loci_traits_top_snps, tax_group == "dummy")


test <- top_loci_traits_top_snps %>%
  filter(snp_id %in% unique(hmp_data$snp_id))


for (phenotype in c("CC_Aug12", "cob_length", "cob_width", "cob_weight")) {
  for (i in c(1:nrow(top_loci_traits_top_snps))){
    #phenotype <- "CC_Aug12"
    #i <- 590 # massilia niabensis
    #i <- 3 # massilia niabensis
    snp <- top_loci_traits_top_snps[i,]
    print(paste(i, snp$snp_id, snp$tax_group))
    
    ## get abundance and yield data
    
    maj_min <- hmp_data %>%
      #ungroup() %>%
      filter(snp_id == snp$snp_id) %>%
      dplyr::select(genotype, allele)
    
    ydat <- microbe_counts_and_yield %>%
      filter(tax_group == snp$tax_group & nitrogen == snp$nitrogen) %>%
      #filter(tax_group == snp$tax_group) %>%
      left_join(maj_min) %>%
      filter( haplotype_NN != "NN")
    
    maj_dat <- c(filter(ydat, allele=="major allele")[, phenotype])[[1]]
    min_dat <- c(filter(ydat, allele=="minor allele")[, phenotype])[[1]]
   
    ttest <- tryCatch(
      {
        t.test(maj_dat, min_dat)
      },
      error = function(e){
        ttest <- "failed"
      }
    )
    
    
    if(ttest != "failed") {
      p <- ttest$p.value
      mean_maj <- ttest$estimate[[1]]
      mean_min <- ttest$estimate[[2]]
    } else {
      p <- NA
      mean_maj <- NA
      mean_min <- NA
    }
    
    snp$phenotype = phenotype
    snp$mean_maj = mean_maj
    snp$mean_min = mean_min
    snp$p = p
    
    yield_ttest <- rbind(yield_ttest, snp)
  }

}


#save(yield_ttest, file="cache/yield_ttest.rda")



### plot yield data p values

load("cache/yield_ttest.rda")


load("largedata/combined_SNP_data.rda")
chr_lengths <- unique(dplyr::select(gwas_dat, chr, chindex))
#chr_len <- c(306970668, 244420837, 235651057, 246967222, 223706090, 173536910, 181718119, 181046068, 159686117, 150929986)
chr_len_10k <- c(306980000, 244430000, 235660000, 246970000, 223710000, 173540000, 181720000, 181050000, 159690000, 150930000)
#vlines <- c(0,cumsum(chr_len)) ## draw vertical lines at chomosome ends
nCHR <- length(unique(gwas_dat$chr))
chindex <- c(0) ## this number needs to be added to position within each chromosome
for (i in c(1:nCHR)){
  #print(i)
  chindex[i+1] <- sum(chr_len_10k[1:i])
}
chindex <- chindex[1:10]
bpadd <- data.frame("chr" = c(1:nCHR), "chindex" = chindex)


# add bins


plot_dat <- yield_ttest %>%
  mutate(log10p_welch = -log10(p)) #%>%
  #left_join(bpadd, by = "chr") %>%
  #mutate(psabs = chindex + ps) %>%
  #mutate(psabs = as.numeric(gsub('.{4}$', '', psabs))) %>% ## pos. in 10kb units
  #mutate(bin = as.numeric(gsub(paste0(".{",4,"}$"), '', ps)) + 1) %>% # absolute bin position
  #group_by(chr, bin, nitrogen) %>%
  
  #left_join(dplyr::select(data_summary_150_traits, trait, tax_group)) %>%
  #mutate(color = ifelse(chr %in% c(1,3,5,7,9), "#276FBF", "#183059")) %>%
 # mutate(chr = as.character(chr)) #%>%
#mutate(ropp_sign_snps = frollmean(sign_snps, n = 10, align = "center")) 





colors = c("1" = "#276FBF", "3" = "#276FBF", "5" = "#276FBF", "7" = "#276FBF", "9" = "#276FBF",
           "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059", "sig" = "#ff0000")


## find position for labels
axis_set <- plot_dat %>% 
  group_by(chr) %>% 
  summarize(center = (max(psabs) + min(psabs)) / 2) # find middle of chromosome (to place labels)



chr_boundaries <- c(0)

chr_len_10k <- c(306980000, 244430000, 235660000, 246970000, 223710000, 173540000, 181720000, 181050000, 159690000, 150930000)

for (i in c(1:10)) {
  sum <- sum(chr_len_10k[1:i])
  chr_boundaries <- c(chr_boundaries, sum/10000)
}




overview_plot <-  ggplot(plot_dat, aes(x = psabs, y = log10p_welch, color = phenotype)) +
  #geom_area() +
  #geom_bar(stat="identity", position="dodge") +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") + 
  #scale_x_continuous(label = axis_set$chr, breaks = axis_set$center, limits = c(1273000000, 1273300000)) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  #geom_vline(xintercept = vlines) +
  #scale_color_manual(values = colors) +
  #scale_size_continuous(breaks = c(1, 5, 10)) +
  #scale_size_continuous(range = c(0.5,3)) +
  #scale_y_continuous(trans='log10') +
  #scale_y_log10() +
  ##coord_trans(y='log10') +
  scale_color_manual(values = c("CC_Aug12"="#008000", "cob_weight"="#ffa701", "cob_length"="#b03b37", "cob_width"="#6e3144")) +
  facet_wrap(~nitrogen, nrow = 2) +
  labs(x = NULL, y = "association with yield component traits \n [-log10p major allele vs. minor allele]") + 
  #ggtitle("total signif. SNPs for 399 10kb bins with nearby annotated genes") +
  geom_vline(xintercept =  chr_boundaries, color = "black", size=0.1) + # to combine plots
  #geom_hline(yintercept = 0, size=0.1) + 
  theme_bw() +
  theme( 
    #legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

overview_plot







### find genes in loci that respond to yield traits, do GO enrichment analysis


### import gene map
gene_pos <- read.table("data/gene_positions.txt", sep="\t")
colnames(gene_pos) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gene_pos$gene_len = gene_pos$end-gene_pos$start



shortlist <- plot_dat %>%
  mutate(ext_bin_start = (bin-2)*10000+1, ext_bin_end = (bin+1)*10000)

### find all genes in there. for each bin

shortlist$genes <- NA
gene_list <- c()
for(i in c(1:nrow(shortlist))){
  #i <- 4
  chr <- shortlist$chr[i]
  print(i)
  range_low <- (shortlist$bin[i]*10000) - 100000
  range_hi <- (shortlist$bin[i]*10000) + 100000
  gene_window <- gene_pos %>%
    filter(seqid == chr) %>%
    filter(start > range_low & end < range_hi)
  
  if (nrow(gene_window) > 0){
    genes <- c()
    for (gene in c(1:nrow(gene_window))){
      
      overlap <- sum(c(gene_window$start[gene]:gene_window$end[gene]) %in% c(shortlist$ext_bin_start[i]:shortlist$ext_bin_end[i]))
      #print(overlap)
      if(overlap > 0){
        #print("gene is near bin")
        gene <- str_match(gene_window$attributes[gene], "gene_id=(.+);")[2]
        #genes <- c(genes, str_match(gene_window$attributes[gene], "gene_id=(.+);")[2])
        row <- shortlist[i,]
        row$genes <- gene
        print(gene)
        gene_list <- rbind(gene_list, row)
      }
    }
    shortlist$genes[i] <- ifelse(is.null(genes), NA, paste(genes, collapse = ';'))
  }
}


### look up gene info ###

gene_info <- read.table("data/Zea_mays.AGPv4_anno.txt", sep="\t", fill=TRUE, header=TRUE)
gene_info <- rename(gene_info, gene = ensembl_gene_id_v4)

gene_list_annotated <- gene_list %>%
  group_by(genes, nitrogen, chr, bin) %>%
  add_tally(name = "assoc_traits") %>%
  rename(gene = genes)

gene_list_annotated <- left_join(gene_list_annotated, gene_info)


#save(gene_list_annotated, file="data/yield_analysis/gene_list_annotated.rda")


load("data/yield_analysis/gene_list_annotated.rda")


yield_genes_stdN <- gene_list_annotated %>%
  ungroup() %>%
  filter(nitrogen == "+N") %>%
  dplyr::select(gene) %>%
  unique()
 
write_csv(yield_genes_stdN, file = "data/yield_analysis/yield_genes_stdN.csv")


yield_genes_lowN <- gene_list_annotated %>%
  ungroup() %>%
  filter(nitrogen == "-N") %>%
  dplyr::select(gene) %>%
  unique()

write_csv(yield_genes_lowN, file = "data/yield_analysis/yield_genes_lowN.csv")


yield_genes_all <- gene_list_annotated %>%
  ungroup() %>%
  dplyr::select(gene) %>%
  unique()

write_csv(yield_genes_all, file = "data/yield_analysis/yield_genes_all.csv")





agrigo <- read_delim("data/yield_analysis/agrigo_all.txt", delim = "\t")



###

##plot log10p_welch vs differential abundance

load("data/group_data.rda")

diffab <- group_data %>%
  dplyr::select(tax_group, log2FoldChange)

diffab_vs_p <- plot_dat %>%
  left_join(diffab)


ggplot(diffab_vs_p, aes(x=log2FoldChange, y=log10p_welch, color=phenotype)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") + 
  scale_color_manual(values = c("CC_Aug12"="#008000", "cob_weight"="#ffa701", "cob_length"="#b03b37", "cob_width"="#6e3144")) +
  facet_wrap(~nitrogen, nrow = 2) +
  theme_bw() +
  theme( 
    #legend.position = "none",
    panel.border = element_blank(),
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    #panel.grid.major.y = element_blank(),
    #panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )




### plot yield vs abundance







### plot selected


short_dat <- plot_dat %>%
  filter(nitrogen == "+N") %>%
  filter(chr == "8")


top_SNPs_lowN<- c("S1_213327012", "S2_27278662", "S3_59034734", "S3_225393869", "S5_36861907", "S7_106821585", "S8_67482384", "S10_9474481", "S10_9499309")

top_SNPs_stdN<- c("S1_244463548", "S2_236874406", "S5_50336991", "S8_119534612")
  

top_yield_signals_lowN <- plot_dat %>%
  filter(nitrogen == "-N") %>%
  filter(snp_id %in% top_SNPs_lowN) %>%
  filter(p <= 0.05)

colors <- c("major allele"="#ecb602", "minor allele"="#a45ee5")

snp <- filter(top_loci_traits_top_snps, snp_id == "S8_119534612")

maj_min <- hmp_data %>%
  #ungroup() %>%
  filter(snp_id == snp$snp_id) %>%
  dplyr::select(genotype, allele)


ydat <- microbe_counts_and_yield %>%
  filter(tax_group == snp$tax_group) %>%
  #filter(tax_group == snp$tax_group & nitrogen == snp$nitrogen) %>%
  #filter(tax_group == snp$tax_group) %>%
  left_join(maj_min) %>%
  filter( haplotype_NN != "NN") %>%
  dplyr::select(-ExG_Aug12) %>%
  pivot_longer(cols = c("cob_length", "cob_width", "cob_weight", "CC_Aug12"), names_to = "yield_trait", values_to = "value")


ggplot(ydat, aes(x=nitrogen, y=value, fill=allele)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(aes(group = allele), label = "p.signif") +
  #ylab("relative microbe abundance") +
  #ylab("ASV count") +
  scale_fill_manual(values = colors) +
  facet_wrap(~yield_trait, scales = "free", nrow = 1) +
  ggtitle(paste0(snp$snp_id," | ", snp$nitrogen)) +
  theme_bw()

## plot microbe abundance



countdat <- microbe_counts_and_yield %>%
  filter(tax_group %in% snp$tax_group) %>%
  #filter(tax_group == snp$tax_group & nitrogen == snp$nitrogen) %>%
  #filter(tax_group == snp$tax_group) %>%
  left_join(maj_min) %>%
  filter( haplotype_NN != "NN") %>%
  dplyr::select(-ExG_Aug12) #%>%
  #pivot_longer(cols = c("cob_length", "cob_width", "cob_weight", "CC_Aug12"), names_to = "yield_trait", values_to = "value")


ggplot(countdat, aes(x=nitrogen, y=log_relab, fill=allele)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(aes(group = allele), label = "p.signif") +
  #ylab("relative microbe abundance") +
  #ylab("ASV count") +
  scale_fill_manual(values = colors) +
  facet_wrap(~tax_group, scales = "free", nrow = 1) +
  ggtitle(paste0(snp$snp_id," | ", snp$nitrogen)) +
  theme_bw()







# annotate major and minor allele genotypes














