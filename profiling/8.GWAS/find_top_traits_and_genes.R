### plotted total significant snps across all 150 traits and both N treatments
### to find strongest signals
### identified 53 (10kb) bins with > 25 significant snps and that were
### associated with gene models
### 53 top regions correspond to 28 traits and 70 genes
library("tidyverse")
library("grid")


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




View(filter(top_bins, tax_group == "Pseudomonas chlororaphis"))

	

### plot essential data for top 22 tax groups:


#1) total abundance (ASV counts)
#2) differential abundance +N vs -N
#3) top GWAS signals +N / -N
#4) +/- selection under +N/-N

load("data/group_data.rda")





top_stdN <- sort(unique(filter(top_bins, nitrogen == "+N")$tax_group))

top_lowN <- sort(unique(filter(top_bins, nitrogen == "-N")$tax_group))




data_summary <- group_data


## mark 22 traits with strongest GWAS results in +N / -N
data_summary$top_stdN <- NA
data_summary$top_lowN <- NA
data_summary[data_summary$tax_group %in% top_stdN, ]$top_stdN <- "x"
data_summary[data_summary$tax_group %in% top_lowN, ]$top_lowN <- "x"


## add differential abundance

data_summary$relab_stdN <- data_summary$count_stdN/sum(data_summary$count_stdN)
data_summary$relab_lowN <- data_summary$count_lowN/sum(data_summary$count_lowN)
data_summary$log2FC <- log2(data_summary$count_stdN/data_summary$count_lowN)
data_summary$relab_log2FC <- log2(data_summary$relab_stdN/data_summary$relab_lowN)





## add selection data

gtcb_HN <- read_csv("cache/gtcb_HN_n150.csv")
gtcb_HN <- subset(gtcb_HN, par %in% "S")
gtcb_HN <- gtcb_HN[, c("Mean", "trait")]
colnames(gtcb_HN) <- c("stdN_Mean_S", "ASV")


gtcb_LN <- read_csv("cache/gtcb_LN_n150.csv")
gtcb_LN <- subset(gtcb_LN, par %in% "S")
gtcb_LN <- gtcb_LN[, c("Mean", "trait")]
colnames(gtcb_LN) <- c("lowN_Mean_S", "ASV")


data_summary <- left_join(data_summary, gtcb_HN)
data_summary <- left_join(data_summary, gtcb_LN)


data_summary_150_traits <- data_summary

save(data_summary_150_traits, file = "data/data_summary_150_traits.rda")

#load("data/data_summary_150_traits.rda")

# subset to taxa


common_names <- c("ASV", "Kingdom", "Phylum", "Class", "Order",
  "Family", "tax_group", "unique_ASVs", "count",
  "H2", "trait", "mean_blup", "top",
  "relab", "log2FC", "relab_log2FC", "Mean_S", "nitrogen" )


data_summary_stdN <- data_summary_150_traits %>%
  select(ASV:unique_ASVs, count_stdN, H2_stdN, trait, mean_blup_stdN, top_stdN, relab_stdN, log2FC, relab_log2FC, stdN_Mean_S)
data_summary_stdN$nitrogen <- "+N"
colnames(data_summary_stdN) <- common_names


data_summary_lowN <- data_summary_150_traits %>%
  select(ASV:unique_ASVs, count_lowN, H2_lowN, trait, mean_blup_lowN, top_lowN, relab_lowN, log2FC, relab_log2FC, lowN_Mean_S)
data_summary_lowN$nitrogen <- "-N"
colnames(data_summary_lowN) <- common_names


data_summary_long <- rbind(data_summary_stdN, data_summary_lowN)


data_summary_long$top <- ifelse(is.na(data_summary_long$top), 0, 1)


### order by differential abundance
factors <- unique(arrange(data_summary_long, desc(relab_log2FC))[, "tax_group"])

### order by mean
factors <- unique(arrange(data_summary_long, desc(Mean_S))[, "tax_group"])

data_summary_long$tax_group <- factor(data_summary_long$tax_group, levels = factors)


### plot 150 tax groups

plot1 <- ggplot(data_summary_long, aes(x = count/sum(data_summary_long$count)*100, y = tax_group)) +
  geom_bar(stat="identity") +
  xlab("% of core microbiome") +
  theme_minimal() +
  theme(axis.text.y = element_text(size=6), axis.title.y=element_blank())

plot2 <- ggplot(data_summary_long, aes(x = H2, y = tax_group, color=nitrogen)) +
  geom_point() +
  scale_color_manual(values=c("red", "blue")) +
  xlab("Heritability") +
  theme_minimal()+
  #facet_wrap(~nitrogen, ncol=2) +
  theme(legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank())

plot3 <- ggplot(data_summary_long, aes(x = relab_log2FC, y = tax_group, fill = relab_log2FC < 0)) +
  geom_bar(stat="identity") +
  scale_fill_manual(guide = FALSE, breaks = c(TRUE, FALSE), values=c("red", "blue")) +
  xlab("-N <- differential abundance -> +N") +
  theme_minimal() +
  theme(legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank())

plot4 <- ggplot(data_summary_long, aes(x = Mean_S, y = tax_group, color=nitrogen)) +
  geom_point() +
  scale_color_manual(values=c("red", "blue")) +
  xlab("selection (mean S)") +
  theme_minimal() +
  theme(legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank())

plot5 <- ggplot(data_summary_long, aes(nitrogen, tax_group, fill=as.factor(top))) +
  geom_tile(colour="white",size=0.25) +
  coord_equal() +
  theme_minimal() +
  scale_fill_manual(values=c("#f9f9f9", "#008000")) +
  xlab("strong GWAS signal") +
  theme(plot.background=element_blank(),
    panel.border=element_blank(),
    legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank())



grid.newpage()
grid.draw(cbind(
  ggplotGrob(plot1),
  ggplotGrob(plot2),
  ggplotGrob(plot3),
  ggplotGrob(plot4),
  ggplotGrob(plot5),
  size = "last"))


### plot only top 22 taxa


top22_taxa <- unique(as.character(filter(data_summary_long, top == 1)[, "tax_group"]))

data_summary_top <- filter(data_summary_long, tax_group %in% top22_taxa)




plot1 <- ggplot(data_summary_top, aes(x = count/sum(data_summary_long$count)*100, y = tax_group)) +
  geom_bar(stat="identity") +
  xlab("% of core microbiome") +
  theme_minimal() +
  theme(axis.text.y = element_text(size=12), axis.title.y=element_blank())

plot2 <- ggplot(data_summary_top, aes(x = H2, y = tax_group, color=nitrogen)) +
  geom_point() +
  scale_color_manual(values=c("red", "blue")) +
  xlab("Heritability") +
  theme_minimal()+
  #facet_wrap(~nitrogen, ncol=2) +
  theme(legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank())

plot3 <- ggplot(data_summary_top, aes(x = relab_log2FC, y = tax_group, fill = relab_log2FC < 0)) +
  geom_bar(stat="identity") +
  scale_fill_manual(guide = FALSE, breaks = c(TRUE, FALSE), values=c("red", "blue")) +
  xlab("-N <- differential abundance -> +N") +
  theme_minimal() +
  theme(legend.position = "none",
  axis.title.y=element_blank(),
  axis.text.y=element_blank())


plot4 <- ggplot(data_summary_top, aes(x = Mean_S, y = tax_group, color=nitrogen)) +
  geom_point() +
  scale_color_manual(values=c("red", "blue")) +
  xlab("selection (mean S)") +
  theme_minimal() +
  theme(legend.position = "none",
  axis.title.y=element_blank(),
  axis.text.y=element_blank())

plot5 <- ggplot(data_summary_top, aes(nitrogen, tax_group, fill=as.factor(top))) +
  geom_tile(colour="white",size=0.25) +
  coord_equal() +
  theme_minimal() +
  scale_fill_manual(values=c("#f9f9f9", "#008000")) +
  xlab("strong GWAS signal") +
  theme(plot.background=element_blank(),
    panel.border=element_blank(),
    legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank())



grid.newpage()
grid.draw(cbind(
  ggplotGrob(plot1),
  ggplotGrob(plot2),
  ggplotGrob(plot3),
  ggplotGrob(plot4),
  ggplotGrob(plot5),
  size = "last"))



