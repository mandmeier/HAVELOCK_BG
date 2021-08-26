#RNAseq_7tissues.R

#df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_fpm_rounded_origNames_and_Altnames.txt

library("tidyverse")
library("stringr")
library("ggpubr")


## import maize reference genes v5 from maizegdb
gene_list_v5 <- read_csv("largedata/RNAseq/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.csv")
colnames(gene_list_v5) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")


## import maize reference genes v5 from maizegdb
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


rev(colnames(rnadat))[1:20]

## calculate mean gene expression by tissue type
ALL_rna <- rnadat %>%
  dplyr::select(TissueWODate,RawPhenotypeNames, AC148152.3_FG001:zma.MIR529) %>%
  rename(tissue = TissueWODate, genotype = RawPhenotypeNames) %>%
  pivot_longer(-c(tissue, genotype), names_to = "gene_AC", values_to = "expression") %>%
  left_join(gene_names) %>%
  group_by(tissue, genotype) %>%
  summarize(mean_expression = mean(expression))





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

## plot expression of genes for each tissue type (boxplot)

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
  













