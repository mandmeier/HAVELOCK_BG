
library("tidyverse")
library("stringr")

load("data/group_data.rda")
gwas_dat_stdN <- read_csv("data/merged_short_5_out_traits_150_stdN_201109-103909.csv")
gwas_dat_lowN <- read_csv("data/merged_short_5_out_traits_150_lowN_201109-100251.csv")

### combine both N treatments
gwas_dat_stdN$nitrogen <- "+N"
gwas_dat_lowN$nitrogen <- "-N"

gwas_dat <- rbind(gwas_dat_stdN, gwas_dat_lowN)




### find shortlist of all bins with > 10 signif. snps

shortlist <- gwas_dat %>%
  mutate(bin = as.numeric(gsub('.{4}$', '', ps)) + 1) %>%
  select(trait, nitrogen, chr, bin) %>%
  group_by(trait, nitrogen, chr, bin) %>%
  tally(name="sign_snps") %>%
  filter(sign_snps >= 10) %>%
  mutate(ext_bin_start = (bin-2)*10000+1, ext_bin_end = (bin+1)*10000) %>%
  group_by(nitrogen, chr, bin) %>%
  add_tally(name="assoc_traits")
 
### import gene map
gene_pos <- read.table("data/gene_positions.txt", sep="\t")
colnames(gene_pos) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gene_pos$gene_len = gene_pos$end-gene_pos$start


### find all genes in there. for each bin

shortlist$genes <- NA
gene_list <- c()
for(i in c(1:nrow(shortlist))){
  #i <- 8
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

gene_list_annotated <- left_join(gene_list_annotated, group_data, by="trait")


gene_list_annotated <- left_join(gene_list_annotated, gene_info)

save(gene_list_annotated, file ="cache/gene_list_annotated_10kb_10above5.rda")





