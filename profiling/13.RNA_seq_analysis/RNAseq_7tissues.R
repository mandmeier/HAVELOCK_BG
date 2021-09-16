#RNAseq_7tissues.R

#df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_fpm_rounded_origNames_and_Altnames.txt

setwd("~/Desktop/Labwork/Havelock_2019/HAVELOCK_BG")

library("tidyverse")
library("stringr")
library("ggpubr")


## import maize reference genes v4 from maizegdb
gene_list_v4 <- read_delim("largedata/RNAseq/Zea_mays.B73_RefGen_v4.46.chr.txt", delim = "\t")
colnames(gene_list_v4) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")


gene_pos <- gene_list_v4 %>%
  #filter(substr(seqid,1,3) == "chr") %>%
  filter(!(is.na(seqid))) %>%
  mutate(gene_len= end-start) %>%
  filter(type == "gene") #%>%
  #mutate(seqid = substr(seqid, 4, nchar(seqid)))






## find all genes in MAPLs
# get short list of all MAPLs
load("cache/GWAS2/v5_MAPLs.rda")


shortlist <- v5_MAPLs %>%
  ungroup() %>%
  dplyr::select(nitrogen, chr, bin) %>%
  unique() %>%
  mutate(ext_bin_start = (bin-2)*10000+1, ext_bin_end = (bin+1)*10000)

LN <- filter(shortlist, nitrogen == "LN")
HN <- filter(shortlist, nitrogen == "HN")





test <- v5_MAPLs %>%
  filter(tax_group == "Niabella yanshanensis")


### find all genes in there. for each bin

shortlist$genes <- NA
gene_list <- c()
for(i in c(1:nrow(shortlist))){
  #i <- 1
  chr <- shortlist$chr[i]
  print(paste("MAPL",i))
  range_low <- (shortlist$bin[i]*10000) - 100000
  range_hi <- (shortlist$bin[i]*10000) + 100000
  gene_window <- gene_pos %>%
    filter(seqid == chr) %>%
    filter(start > range_low & end < range_hi)
  
  if (nrow(gene_window) > 0){
    genes <- c()
    for (gene in c(1:nrow(gene_window))){
      #gene <- 1
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



MAPL_genes <- gene_list

write_csv(MAPL_genes, file="cache/MAPL_genes.csv")



unique(MAPL_genes$gene) # 395 genes


unique(MAPL_genes$gene_AC) # 307 genes have RNA seq data




## for all genes check expression in root vs other tissues


## translate Zm gene names into AC gene names
gene_names <- read_csv("largedata/RNAseq/V3_V4.csv", col_names = FALSE)

gene_names <- gene_names %>%
  rename(gene_AC=X1, gene=X2) %>%
  dplyr::select(gene, gene_AC)

MAPL_genes <- MAPL_genes %>%
  rename(gene=genes) %>%
  left_join(gene_names)

save(MAPL_genes, file="cache/GWAS2/MAPL_genes.rda")

load("cache/GWAS2/MAPL_genes.rda")


#get gene expression data
rnadat <- read_delim(file="largedata/RNAseq/rna_seq.txt", delim="\t")


## these MAPL genes have RNA seq data available
MAPL_genes_rna <- unique(MAPL_genes$gene_AC)[unique(MAPL_genes$gene_AC) %in% colnames(rnadat)]



## select MAPL gene columns from rna dat


MAPL_rna <- rnadat %>%
  dplyr::select(TissueWODate,RawPhenotypeNames, MAPL_genes_rna) %>%
  rename(tissue = TissueWODate, genotype = RawPhenotypeNames) %>%
  pivot_longer(-c(tissue, genotype), names_to = "gene_AC", values_to = "expression") %>%
  left_join(gene_names) %>%
  group_by(tissue, genotype) %>%
  summarize(mean_expression_MAPL = mean(expression))




colnames(rnadat)[18000]

rev(colnames(rnadat))[1:20]

## calculate mean gene expression by tissue type
ALL_rna <- rnadat %>%
  filter(Tissue != "L3Mid") %>%
  dplyr::select(TissueWODate,RawPhenotypeNames, AC148152.3_FG001:zma.MIR529) %>%
  #dplyr::select(TissueWODate,RawPhenotypeNames, AC148152.3_FG001:GRMZM2G009080) %>%
  #dplyr::select(TissueWODate,RawPhenotypeNames, AC148152.3_FG001:GRMZM2G123843) %>%
  rename(tissue = TissueWODate, genotype = RawPhenotypeNames) %>%
  pivot_longer(-c(tissue, genotype), names_to = "gene_AC", values_to = "expression") %>%
  left_join(gene_names) %>%
  mutate(group=ifelse(gene %in% MAPL_genes$gene,"MAPL", "other")) %>%
  group_by(tissue, genotype, group) %>%
  summarize(mean_expression = mean(expression), sd_expression = sd(expression))



unique(filter(ALL_rna, group == "MAPL")$gene) # 305 MAPL genes
unique(filter(ALL_rna, group == "other")$gene) # 29540 other genes


unique(ALL_rna$genotype)

unique(ALL_rna$tissue)
ALL_rna$tissue <- factor(ALL_rna$tissue, levels = c("GRoot", "GShoot", "L3Base", "L3Tip", "LMAD", "LMAN", "Kern"))

colors <- c(MAPL="#9cd6ff", other="#fff2bd")

ggplot(ALL_rna, aes(x=tissue, y=mean_expression, fill=group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(aes(group = group), label = "p.signif") +
  ylab("mean gene expression [FPKM]") +
  ggtitle("mean gene expression by tissue type") +
  scale_fill_manual(name = "", labels = c("MAPL genes", "other genes"), values = colors) +
  #scale_fill_manual(values = colors) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = c(0.85, 0.7))


# MAPL genes(n=305) vs other genes (n=36822)
# FPKM (Fragments Per Kilobase Million)


ALL_rna



### correlation expression vs microbe abundance
### for each gene, tissue, genotype, microbe, N treatment


load("data/yield_analysis/abundance_vs_phenotype.rda") 


## mean of all samples for each genotype
mean_expression <- rnadat %>%
  filter(Tissue != "L3Mid") %>%
  dplyr::select(TissueWODate,RawPhenotypeNames, AC148152.3_FG001:zma.MIR529) %>%
  #dplyr::select(TissueWODate,RawPhenotypeNames, AC148152.3_FG001:GRMZM2G009080) %>%
  #dplyr::select(TissueWODate,RawPhenotypeNames, AC148152.3_FG001:GRMZM2G123843) %>%
  rename(tissue = TissueWODate, genotype = RawPhenotypeNames) %>%
  pivot_longer(-c(tissue, genotype), names_to = "gene_AC", values_to = "expression") %>%
  left_join(gene_names) %>%
  mutate(group=ifelse(gene %in% MAPL_genes$gene,"MAPL", "other")) %>%
  filter(group == "MAPL") %>%
  group_by(tissue, genotype, gene_AC, gene) %>%
  summarize(mean_expression = mean(expression))

abundance <- abundance_vs_phenotype %>%
  dplyr::select(nitrogen, MM_name, tax_group, blup_logrel) %>%
  rename(genotype = MM_name) %>%
  unique()

expression_vs_abundance <- mean_expression %>%
  left_join(abundance) %>%
  drop_na() 


unique(expression_vs_abundance$tissue)

unique(expression_vs_abundance$gene)

unique(expression_vs_abundance$genotype)


save(expression_vs_abundance, file = "largedata/expression_vs_abundance.rda")



### investigate microbes in top 90% quantile MAPLs

load("cache/RNAseq/top90pc_taxa.rda")

top90pc_taxa


# expression vs abundance

eva <- expression_vs_abundance %>%
  filter(tissue == "GRoot" & tax_group %in% top90pc_taxa) %>%
  filter(mean_expression > 10) %>%
  group_by(tissue, gene, nitrogen, tax_group) %>%
  add_tally(name = "no_genotypes") %>%
  filter(no_genotypes >= 50)



save(eva, file="cache/RNAseq/eva.rda")


corr_data <- eva %>%
  ungroup() %>%
  dplyr::select(nitrogen,tax_group, tissue, gene) %>%
  unique()

corr_data$pearson <- 0
corr_data$p_value <- 0

get_correlation <- function(nitr, gen, tgr, tiss){
  #nitr <- "LN"
  #tiss <- "GRoot"
  #tgr <- "Sphingoaurantiacus"
  #var <- "blup_logrel"
  #gene <- corr_data[i,]$gene
  #gen <- "Zm00001d014467"
  nmp <- data.frame(filter(eva,  nitrogen == nitr & gene == gen & tax_group == tgr & tissue == tiss))
  cor <- cor.test(nmp[, "blup_logrel"], nmp[, "mean_expression"], method=c("pearson"))
  return(cor)
}


c(1:nrow(corr_data))


for (i in c(1:nrow(corr_data))){
  print(i)
  cor <- get_correlation(corr_data[i,]$nitrogen, corr_data[i,]$gene, corr_data[i,]$tax_group, corr_data[i,]$tissue)
  #print(cor)
  pears <- cor$estimate[[1]]
  #print(pears)
  p_value <- cor$p.value
  #print(p_value)
  corr_data[i,]$pearson <- pears
  corr_data[i,]$p_value <- p_value
}


corr_exp_vs_abundance_GRoot <- corr_data

save(corr_exp_vs_abundance_GRoot, file="cache/RNAseq/corr_exp_vs_abundance_GRoot.rda")



ggplot(corr_exp_vs_abundance_GRoot, aes(x=pearson, y=-log10(p_value), color=nitrogen)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 2, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#000080", "#C00001")) +
  facet_wrap(~nitrogen, nrow = 1) +
  xlab("pearson correlation") +
  theme_bw() +
  theme(panel.spacing.x = unit(1.5, "lines"))


filter(corr_exp_vs_abundance_GRoot, nitrogen == "LN")

ggplot(corr_exp_vs_abundance_GRoot, aes(x=pearson, y=-log10(p_value), color=nitrogen)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 2, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#000080", "#C00001")) +
  facet_wrap(~tax_group) +
  xlab("pearson correlation") +
  theme_bw() +
  theme(panel.spacing.x = unit(1.5, "lines"))

sign_LN <- corr_exp_vs_abundance_GRoot %>%
  filter(p_value <= 0.01, nitrogen == "LN")

sign_HN <- corr_exp_vs_abundance_GRoot %>%
  filter(p_value <= 0.01, nitrogen == "HN")


Steroidobacter_LN <- eva %>%
  filter(tax_group == "Steroidobacter" & nitrogen == "LN", gene %in% c("Zm00001d033633", "Zm00001d039079", "Zm00001d041489"))

#### plot correlation
ggscatter(Steroidobacter_LN, x = "blup_logrel", y = "mean_expression", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", alpha=0.25) +
  ylab("mean expression") +
  xlab("abundance") +
  ggtitle("Steroidobacter | LN") +
  facet_wrap(~gene, scales = "free", nrow = 1) +
  theme_bw()


Massilia_putida_HN <- eva %>%
  filter(tax_group == "Massilia putida" & nitrogen == "HN", gene %in% c("Zm00001d013860", "Zm00001d027941", "Zm00001d040697"))

#### plot correlation
ggscatter(Massilia_putida_HN, x = "blup_logrel", y = "mean_expression", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", alpha=0.25) +
  ylab("mean expression") +
  xlab("abundance") +
  ggtitle("Massilia putida | HN") +
  facet_wrap(~gene, scales = "free", nrow = 1) +
  theme_bw()


RB41_sp_1_LN <- eva %>%
  filter(tax_group == "RB41 sp 1" & nitrogen == "LN", gene %in% c("Zm00001d027944", "Zm00001d046322"))

#### plot correlation
ggscatter(RB41_sp_1_LN, x = "blup_logrel", y = "mean_expression", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", alpha=0.25) +
  ylab("mean expression") +
  xlab("abundance") +
  ggtitle("RB41 sp 1 | LN") +
  facet_wrap(~gene, scales = "free", nrow = 1) +
  theme_bw()









cor <- cor.test(Pelomonas_LN$blup_logrel, Pelomonas_LN$mean_expression, method=c("pearson"))


get_correlation <- function(nitr, gen, tgr, tiss){
  nitr <- "LN"
  tiss <- "GRoot"
  tgr <- "Pelomonas saccharophila"
  var <- "blup_logrel"
  gen <- "Zm00001d002476"
  nmp <- data.frame(filter(eva,  nitrogen == nitr & gene == gen & tax_group == tgr & tissue == tiss))
  if(nrow(nmp) >= 5){
    cor <- cor.test(nmp[, "blup_logrel"], nmp[, "mean_expression"], method=c("pearson"))
  } else {
    cor = "missing"
  }
  return(cor)
}



pelo_data <- eva %>%
  ungroup() %>%
  dplyr::select(nitrogen,tax_group, tissue, gene) %>%
  unique() %>%
  filter(tax_group == "Pelomonas saccharophila")

pelo_data$pearson <- 0
pelo_data$p_value <- 0



get_correlation <- function(nitr, gen, tgr, tiss){
  #nitr <- "LN"
  #tiss <- "GRoot"
  #tgr <- "Sphingoaurantiacus"
  #var <- "blup_logrel"
  #gene <- corr_data[i,]$gene
  #gen <- "Zm00001d014467"
  nmp <- data.frame(filter(eva,  nitrogen == nitr & gene == gen & tax_group == tgr & tissue == tiss))
  if(nrow(nmp) >= 5){
    cor <- cor.test(nmp[, "blup_logrel"], nmp[, "mean_expression"], method=c("pearson"))
  } else {
    cor = "missing"
  }
  return(cor)
}


c(1:nrow(pelo_data))


for (i in c(1:nrow(pelo_data))){
  print(i)
  cor <- get_correlation(pelo_data[i,]$nitrogen, pelo_data[i,]$gene, pelo_data[i,]$tax_group, pelo_data[i,]$tissue)
  if(cor == "missing"){
    pears <- NA
    p_value <- NA
  } else {
    pears <- cor$estimate[[1]]
    p_value <- cor$p.value
  }
  pelo_data[i,]$pearson <- pears
  pelo_data[i,]$p_value <- p_value
}






Pelomonas_LN <- eva %>%
  filter(tax_group == "Pelomonas saccharophila" & nitrogen == "LN", gene %in% c("Zm00001d002476"))




#### plot correlation
ggscatter(Pelomonas_LN, x = "blup_logrel", y = "mean_expression", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", alpha=0.25) +
  ylab("mean expression") +
  xlab("abundance") +
  ggtitle("Pelomonas saccharophila | LN") +
  #facet_wrap(~gene, scales = "free", nrow = 1) +
  theme_bw()






#save(corr_data, file="data/yield_analysis/corr_data.rda")









test <- ALL_rna %>%
  dplyr::select(-sd_expression) %>%
  pivot_wider(names_from = group, values_from = mean_expression) %>%
  mutate(diff_expression = MAPL/other)




ggplot(test, aes(x=tissue, y=diff_expression)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  #stat_compare_means(aes(group = group), label = "p.signif") +
  #scale_fill_manual(values = colors) +
  ggtitle("mean gene expression by tissue type, all genes (n=37127) vs. MAPL genes(n=305)") +
  theme_bw()



ALL_rna$group <- factor(ALL_rna$group)
str(ALL_rna)

mean_expr <- MAPL_rna %>%
  left_join(ALL_rna) %>%
  rename(MAPL_genes=mean_expression_MAPL, all_genes=mean_expression) %>%
  pivot_longer(-c(tissue, genotype), names_to = "gene_set", values_to = "mean_expression")
  





ggplot(mean_expr, aes(x=tissue, y=mean_expression, fill=gene_set)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(aes(group = gene_set), label = "p.signif") +
  #scale_fill_manual(values = colors) +
  ggtitle("mean gene expression by tissue type, all genes (n=37127) vs. MAPL genes(n=305)") +
  theme_bw()
  






#save(MAPL_rna, file="cache/RNAseq/MAPL_rna.rda")

load("cache/RNAseq/MAPL_rna.rda")


## Are MAPL genes preferentially expressed in roots?
## plot expression of genes for each tissue type (boxplot)
## MAPL root expression vs mean root expression






gt0 <- filter(MAPL_rna, expression > 0)


ggplot(gt0, aes(x=tissue, y=expression)) +
  geom_bar(stat = "identity")








head(rnadat)


colnames(rnadat)

unique(rnadat$RNA_TaxaName)

sort(unique(rnadat$RawPhenotypeNames))


abundance <- read_csv("data/abundance_logrel.csv")

unique(abundance$genotype)

names <- read_csv("data/BG_MM_Gen_names.csv")


rnaseq_names <- names %>% 
  dplyr::select(MM_name, GX_name, RawPhenotypeNames) %>%
  dplyr::filter(GX_name %in% unique(abundance$genotype)) %>%
  drop_na ()



rnadat_210 <- dplyr::filter(rnadat, RawPhenotypeNames %in% rnaseq_names$RawPhenotypeNames)


colnames(rnadat_210)

#save(rnadat_210, file = "largedata/RNAseq/rnadat_210.rda")
  



#### Validate with Gen's RNA seq data



gendat <- read_csv("data/RNA_seq/FPKM_Old_new_four_lines2.csv")



MAPL_genes_zm <- unique(MAPL_genes$gene)



rna_gen <- gendat %>%
  rename(gene = X1) %>%
  pivot_longer(-gene, names_to = "id", values_to = "expression") %>%
  separate(id, c("nitrogen", "genotype", "tissue", "rep"), "_") %>%
  mutate(group = ifelse(gene %in% MAPL_genes_zm, "MAPL", "other")) %>%
  group_by(nitrogen, genotype, tissue, group) %>%
  summarize(mean_expression = mean(expression), sd_expression = sd(expression))



unique(filter(rna_gen, group == "MAPL")$gene) # 395 MAPL genes
unique(filter(rna_gen, group == "other")$gene) # 45877 other genes


length(unique(filter(rna_gen, group == "other")$gene))
       
       
       
rna_gen$tissue <- factor(rna_gen$tissue, levels = c("Root", "leaf"))





colors <- c(MAPL="#9cd6ff", other="#fff2bd")

ggplot(rna_gen, aes(x=tissue, y=mean_expression, fill=group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(aes(group = group), label = "p.signif") +
  #facet_wrap(~nitrogen, nrow = 1) +
  ylab("mean gene expression [FPKM]") +
  ggtitle("mean gene expression by tissue type") +
  scale_fill_manual(name = "", labels = c("MAPL genes", "other genes"), values = colors) +
  #scale_fill_manual(values = colors) +
  theme_bw() +
  theme(axis.title.x = element_blank())






