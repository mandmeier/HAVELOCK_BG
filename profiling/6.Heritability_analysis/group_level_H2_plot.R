library("tidyverse")


load("data/group_data.rda")

#### plot H2 vs. BLUP abundance ####



colors <- c("Burkholderia" = "#0000bd",
  "Sphingobium" = "#f20019", "Ralstonia" = "#70fe00",
  "f__Enterobacteriaceae" = "#8f00c7", "Sphingobacterium" = "#0086fe",
  "Massilia" = "#008000", "Fulvimonas" = "#00fefe",
  "Chitinophaga" = "#fe68fe", "Labrys" = "#fe8420",
  "Chryseobacterium" = "#827800", "f__Rhodanobacteraceae" = "#fefe00",
  "other" = "#eeeeee")




### h2 vs BLUP

## mark anything in 90% quantile for H2 or BLUP

### stdN
q90_H2_stdN <- quantile(group_data$H2_stdN[group_data$H2_stdN!=min(group_data$H2_stdN)], 0.90)
q90_blup_stdN <- quantile(group_data$mean_blup_stdN[group_data$mean_blup_stdN!=min(group_data$mean_blup_stdN)], 0.90)
q90_stdN <- group_data %>%
  filter(H2_stdN >= q90_H2_stdN | mean_blup_stdN >= q90_blup_stdN)
top_stdN <- sort(unique(q90_stdN$tax_group))

### lowN
q90_H2_lowN <- quantile(group_data$H2_lowN[group_data$H2_lowN!=min(group_data$H2_lowN)], 0.90)
q90_blup_lowN <- quantile(group_data$mean_blup_lowN[group_data$mean_blup_lowN!=min(group_data$mean_blup_lowN)], 0.90)
q90_lowN <- group_data %>%
  filter(H2_lowN >= q90_H2_lowN | mean_blup_lowN >= q90_blup_lowN)
top_lowN <- sort(unique(q90_lowN$tax_group))

top_taxa <- unique(c(top_stdN, top_lowN))


### add label if in top_taxa

stdN <- group_data %>%
  select(tax_group, H2_stdN, mean_blup_stdN) %>%
  mutate(nitrogen = "+N") %>%
  rename(H2 = H2_stdN, mean_blup = mean_blup_stdN)

lowN <- group_data %>%
  select(tax_group, H2_lowN, mean_blup_lowN) %>%
  mutate(nitrogen = "-N") %>%
  rename(H2 = H2_lowN, mean_blup = mean_blup_lowN)


combined <- rbind(stdN,lowN)

po <- combined %>%
  filter(H2 > 0)


## assign colors ## currently colors don't work together with ggrepel...
po$tax_group_other <- ifelse(po$tax_group %in% top_taxa, po$tax_group, "other")
colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
top_taxa_other <- c(top_taxa, "other")
testcolors <- names(top_taxa_other) <- colors37[1:length(top_taxa_other)]



library("ggrepel")

h2_vs_blup <- ggplot(po, aes(x = mean_blup, y = H2, label=tax_group_other)) +
  geom_point() +
  geom_label_repel(aes(label=ifelse(tax_group_other %in% top_taxa, tax_group_other, '')), size=2.5, min.segment.length=0.1) +
  scale_color_manual(values = testcolors) +
  facet_wrap(~nitrogen, ncol = 2) +
  theme_bw()
  
h2_vs_blup







