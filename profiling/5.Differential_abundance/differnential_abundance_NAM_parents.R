### Compare abundance of top 22 taxa in B73 vs NAM founder lines
### where differences exist, find recombinant inbred lines RIL for follow-up experiment

library("tidyverse")

### load countd by genotype data

load("cache/counts_per_genotype.rda")


### subset B73 qnd NAM founder genotypes

## from https://science.sciencemag.org/content/325/5941/737/tab-figures-data

NAM_missing <- c("CML247", "KI3", "Ky21 = KI2021?", "Mo18W", "Tx303" )

NAM_set <- c("B73", "B97", "CML103", "CML228", "CML277", "CML322",
  "CML333", "CML52", "CML69", "HP301", "IL14H", "KI11",
  "M162W", "M37W", "MS71", "NC350", "NC358", "OH43",
  "OH7B", "P39", "TZI8")


NAM <- count_data %>%
  filter(genotype %in% NAM_set)

unique(NAM$genotype)

### subset top22 tax groups

load("data/data_summary_150_traits.rda")

top_taxa <- data_summary_150_traits %>%
  filter(top_stdN == "x" | top_lowN == "x")

top_taxa <- unique(top_taxa$tax_group)


NAM_top <- filter(NAM, tax_group %in% top_taxa)




### calculate abundance under stdN and lowN as percentage of total counts

relab_data_top <- NAM_top %>%
  mutate(nitrogen = ifelse(nitrogen == "+N", "stdN", "lowN")) %>%
  pivot_wider(names_from = nitrogen, values_from = count) %>%
  mutate(total_count = stdN + lowN) %>%
  mutate(stdN_norm = stdN/total_count) %>%
  mutate(lowN_norm = -lowN/total_count) %>%
  pivot_longer(cols= ends_with("norm"), names_to="nitrogen", values_to="abundance") %>%
  mutate(pedigree = paste(genotype, paste0("(",subpopulation,")")))


sort(unique(relab_data_top$pedigree))

order <- c("B73", ## ss
  "B97", "M162W", "MS71", "OH7B", "OH43", ## nss
  "CML103", "CML228", "CML277", "CML322", "CML333", "CML52", ## ts
  "CML69", "KI11", "NC350",  "NC358", "TZI8", ## ts
  "IL14H", "P39", ## sweet corn
  "HP301", ## popcorn
  "M37W") ## mixed


relab_data_top$genotype <- factor(relab_data_top$genotype, levels = order)


## plot abundance


taxgrp <- "Acinetobacter nosocomialis"

taxa <- relab_data_top %>%
  filter(tax_group %in% c(taxgrp))
  


p <- ggplot(relab_data_top) +
  geom_bar(aes(x=genotype, y=abundance, fill=nitrogen), stat='identity', position = position_dodge(width = 0)) +
  #facet_wrap(~tax_group, ncol = 1, scales = "free") +
  #scale_x_reordered() +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  #ggtitle(taxgrp) +
  facet_wrap(~ tax_group, ncol=3) +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

p




