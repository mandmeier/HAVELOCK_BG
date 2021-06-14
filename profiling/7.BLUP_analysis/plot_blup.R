#plot_blup

library("tidyverse")
library("phyloseq")
library("ggrepel")

load("data/group_data.rda")

blups <- read_csv("data/blup_stdN_150_tax_groups.csv")
blups <- gather(blups, key = "ASV", value = "counts", -genotype)
blups$nitrogen <- "+N"



blupl <- read_csv("data/blup_lowN_150_tax_groups.csv")
blupl <- gather(blupl, key = "ASV", value = "counts", -genotype)
blupl$nitrogen <- "-N"



taxa <- group_data[, c("trait", "tax_group")]



blup <- blups %>%
  rbind(blupl) %>%
  rename(trait = ASV) %>%
  left_join(taxa) %>%
  select(nitrogen, genotype, tax_group, counts)


blup$label <- paste0(ifelse(blup$nitrogen=="+N", "stdN", "lowN"), "_", blup$genotype)




blup_df <- blup[, c("label", "tax_group","counts")] %>%
  pivot_wider(names_from = tax_group, values_from = counts ) %>%
  column_to_rownames(var = "label")

blup_matrix <- na.omit(as.matrix(blup_df))

PCs <- data.frame(prcomp(blup_matrix, scale. = TRUE, center = TRUE)$x)

PC_123 <- PCs %>%
  rownames_to_column(var = "label") %>%
  select(label, PC1, PC2, PC3)

blup_pc <- blup %>%
  left_join(PC_123) %>%
  select(nitrogen, genotype, PC1, PC2, PC3) %>%
  unique()


## get subpopulation data
load("data/ps_grp.rda")
sdat <- data.frame(sample_data(ps_grp))
subpop <- unique(select(sdat, genotype, subpopulation))
rownames(subpop) <- NULL
subpop <- unique(subpop)

blup_pc <- blup_pc %>%
  left_join(subpop)


### plot PC
pc_plot <- ggplot(blup_pc, aes(x=PC1, y=PC2, color=nitrogen)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values = c("#C00001", "#000080")) +
  theme_classic()
pc_plot






### plot +N vs -N mean BLUP

# see also  group_level_H2_230genotypes.R


mean_blup <- blup %>%
  group_by(nitrogen, tax_group) %>%
  summarize(mean_blup = mean(counts, na.rm=TRUE)) %>%
  pivot_wider(names_from = nitrogen, values_from = mean_blup) %>%
  rename(lowN = "-N") %>%
  rename(stdN = "+N")



blup_plot <- ggplot(mean_blup, aes(x=lowN, y=stdN)) +
  geom_point() +
  geom_abline(coef = c(0,1), color = "red") +
  geom_label_repel(aes(label = mean_blup$tax_group),
                   box.padding   = 0.05, 
                   point.padding = 0.5,
                   size = 2,
                   min.segment.length = 0.5,
                   segment.color = 'grey50') +
  theme_classic() +
  ylab("Log(relative abundance) +N") +
  xlab("Log(relative abundance) -N")


blup_plot






