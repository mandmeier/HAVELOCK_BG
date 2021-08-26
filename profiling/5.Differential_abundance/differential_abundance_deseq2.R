### for each subgroup do differential abundance analysis for Species and Nitrogen in RHZ and SOL
### identify groups that respond the most.


library("phyloseq")
library("tidyverse")
library("ggplot2")
library("DESeq2")


load("data/ps_grp.rda")


## avoid special characters
sample_data(ps_grp)$nitrogen <- ifelse(sample_data(ps_grp)$nitrogen == "+N", "stdN", "lowN")

## add pseudocount
ps_grp <- transform_sample_counts(ps_grp, function(x) x + 1)

## make asv table for deseq
asv_table <- as.data.frame(t(as.matrix(otu_table(ps_grp))))

### make deseq2 model
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = asv_table,
  colData = as.matrix(sample_data(ps_grp)),
  design = as.formula(paste("~", "nitrogen")))



## fit deseq2 model
diagdds <-  DESeq(ddsFullCountTable, test="Wald", fitType="parametric")
res <- results(diagdds, cooksCutoff = FALSE)
res <- as.data.frame(res)




## add diffab to group data
load("data/group_data.rda")
res <- res %>%
  select(log2FoldChange, lfcSE, padj) %>%
  rownames_to_column(var = "ASV")

#group_data <- left_join(group_data, res)

#save(group_data, file="data/group_data.rda")


plot_data <- group_data %>%
  mutate(color = ifelse(log2FoldChange < 0, "negative", "positive")) %>%
  mutate(color = ifelse(padj >= 0.05, "n.s.", color))


label_pos <- "more abundant under Std N"
label_neg = "more abundant under low N"
label_ns = "no significant difference"

sum(plot_data$diffab_group == "positive")
sum(plot_data$diffab_group == "negative")
sum(plot_data$diffab_group == "n.s.")



p <- ggplot(plot_data, aes(x = tax_group, y = log2FoldChange)) +
  geom_bar(stat = "identity", aes(fill = color), width=0.8) +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width=0.4) +
  #geom_text(aes(label = signif,  vjust=0.5, hjust = ifelse(plot_data$log2FoldChange >= 0, -1, 2)), colour = "black", fontface = "bold", size = 3) +
  scale_fill_manual(values=c(positive="#000080", negative="#C00001", n.s.="#999999"), labels = c("no significant difference", "more abundant under -N", "more abundant under +N")) +
  theme_minimal() +
  xlab("") +
  #ggtitle(title) +
  labs(fill = "") +
  #scale_y_continuous(limits = c(-10,10)) +
  coord_flip()

p





