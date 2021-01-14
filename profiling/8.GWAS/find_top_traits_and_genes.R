### plotted total significant snps across all 150 traits and both N treatments
### to find strongest signals
### identified 53 (10kb) bins with > 25 significant snps and that were
### associated with gene models
### 53 top regions correspond to 28 traits and 70 genes


load("cache/gene_list_annotated_10kb_10above5.rda")
load("largedata/combined_SNP_data.rda")

head(gene_list_annotated)

head(gwas_dat)

# find bins with at least 10 sig snps

bin_regex <- paste0(".{",4,"}$")

unique(gwas_dat$trait)

# add bins
plot_dat <- gwas_dat %>%
  filter(-log10(p_wald) >= 1) %>%
  mutate(bin = as.numeric(gsub(bin_regex, '', ps)) + 1) %>%
  group_by(chr, bin, nitrogen, trait) %>%
  add_count(count = sum(-log10(p_wald) >= 5)) %>%
  filter(count >= 10) %>%
  filter(bin %in% unique(gene_list_annotated$bin))



## these are the SNP counts for all 509 bins that had at least 10 SNPs above 5 for at least one trait,
## and that overlapped with at least one gene model within +/- 10 kb



snp_counts <- left_join(plot_dat, gene_list_annotated)


snp_counts$psabs <- as.numeric(gsub('.{4}$', '', snp_counts$psabs))





ggplot(snp_counts, aes(x=psabs, y=sign_snps)) +
  geom_point() +
  facet_wrap(~ nitrogen, ncol = 1)


snp_counts$color <- ifelse(snp_counts$chr %in% c(1,3,5,7,9), "#276FBF", "#183059")
snp_counts$chr <- as.character(snp_counts$chr)

## find position for labels
axis_set <- gwas_dat %>% 
  group_by(chr) %>% 
  summarize(center = (max(psabs) + min(psabs)) / 2) # find middle of chromosome (to place labels)


unique(snp_counts$psabs)




colors = c("1" = "#276FBF", "3" = "#276FBF", "5" = "#276FBF", "7" = "#276FBF", "9" = "#276FBF",
  "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059", "sig" = "#ff0000")


manhplot <- ggplot(snp_counts, aes(x = psabs, y = sign_snps, color = chr, shape=nitrogen)) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = 27, color = "red", linetype = "dashed") + 
  #scale_x_continuous(label = axis_set$chr, breaks = axis_set$center, limits = c(1273000000, 1273300000)) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  #geom_vline(xintercept = vlines) +
  scale_color_manual(values = colors) +
  #scale_size_continuous(range = c(0.5,3)) +
  #facet_wrap(~nitrogen, nrow = 2) +
  labs(x = NULL, y = "Total Sign. SNPs") + 
  ggtitle("total signif. SNPs for 399 10kb bins with nearby annotated genes") +
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )

manhplot




ref <- read_csv("data/microbe_references.csv")



### find threshld for bins
### plot number of traits retained for 10-100 sign. SNPs




threshold_plot_data <- data.frame(matrix(ncol = 4, nrow = 0))

for(threshold in seq(10,100)) {
  top <- snp_counts %>%
    filter(sign_snps >= threshold) %>%
    select(nitrogen, chr, bin, sign_snps, assoc_traits, tax_group, gene) %>%
    distinct()
  
  threshold_plot_data <- rbind(threshold_plot_data, c(threshold, length(unique(top$tax_group)), length(unique(top$gene)), length(unique(top$bin))))
  
}

colnames(threshold_plot_data) <- c("sign_snps", "taxa_retained", "genes_retained", "bins_retained")

ggplot(threshold_plot_data, aes(x=sign_snps, y=taxa_retained)) +
  geom_point() +
  geom_vline(xintercept = 27, color = "red", linetype = "dashed") +
  theme_bw()


ggplot(threshold_plot_data, aes(x=sign_snps, y=genes_retained)) +
  geom_point() +
  geom_vline(xintercept = 27, color = "red", linetype = "dashed") +
  theme_bw()


ggplot(threshold_plot_data, aes(x=sign_snps, y=bins_retained)) +
  geom_point() +
  geom_vline(xintercept = 27, color = "red", linetype = "dashed") +
  theme_bw()


## threshold 27 is good



### find bins with > 25 total sign. 


top_bins <- snp_counts %>%
  filter(sign_snps >= 27) %>%
  select(nitrogen, chr, bin, sign_snps, assoc_traits, tax_group, unique_ASVs, count_stdN, count_lowN, H2_stdN, H2_lowN, gene, description, kegg_enzyme, name_1006, definition_1006) %>%
  distinct()


top_bins <- left_join(top_bins, ref)

save(top_bins, file = "cache/top_taxa_bins_genes.rda")

sort(unique(filter(top_bins, nitrogen == "+N")$tax_group))

sort(unique(filter(top_bins, nitrogen == "-N")$tax_group))


# These have a lot of noise, signals questionable
exclude <- c("Chryseobacterium indologenes", "Chryseobacterium joostei", "Delftia tsuruhatensis", )

# These show particularly nice signals
best <- c("Acinetobacter nosocomialis", "Archangium gephyra","Candidatus Udaeobacter copiosus","f_A21b","f_Vicinamibacteraceae", "Filimonas sp 2",
  "Ilumatobacter",  "Massilia niabensis", "Niabella yanshanensis","Oryzihumus terrae","Pajaroellobacter",
  "Rhizobium daejeonense", "Rhizobium mesosinicum", "Stenotrophomonas maltophilia", )


View(filter(top_bins, tax_group == "Pseudomonas chlororaphis"))

	

