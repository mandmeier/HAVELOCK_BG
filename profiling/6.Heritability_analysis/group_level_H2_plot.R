library("tidyverse")


load("data/group_data.rda")

#### plot H2 vs. BLUP abundance ####


### h2 vs BLUP

## mark anything in 95% quantile for H2 or BLUP

### stdN
q95_H2_stdN <- quantile(group_data$H2_stdN[group_data$H2_stdN!=min(group_data$H2_stdN)], 0.95)
q95_blup_stdN <- quantile(group_data$mean_blup_stdN[group_data$mean_blup_stdN!=min(group_data$mean_blup_stdN)], 0.95)
q95_stdN <- group_data %>%
  filter(H2_stdN >= q95_H2_stdN | mean_blup_stdN >= q95_blup_stdN)
top_stdN <- sort(unique(q95_stdN$tax_group))

### lowN
q95_H2_lowN <- quantile(group_data$H2_lowN[group_data$H2_lowN!=min(group_data$H2_lowN)], 0.95)
q95_blup_lowN <- quantile(group_data$mean_blup_lowN[group_data$mean_blup_lowN!=min(group_data$mean_blup_lowN)], 0.95)
q95_lowN <- group_data %>%
  filter(H2_lowN >= q95_H2_lowN | mean_blup_lowN >= q95_blup_lowN)
top_lowN <- sort(unique(q95_lowN$tax_group))

top_taxa <- unique(c(top_stdN, top_lowN))

#save(top_taxa, file = "cache/top_heritable_taxa_150.rda")




load("cache/top_heritable_taxa_150.rda")
load("cache/taxa_colors.rda")

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
po$tax_group_other <- ifelse(po$tax_group %in% names(colors), po$tax_group, "other")
#colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")


#library("ggrepel")



h2_vs_blup <- ggplot(po, aes(x = mean_blup, y = H2, color=tax_group_other)) +
  geom_point() +
  #geom_label_repel(aes(label=ifelse(tax_group_other %in% top_taxa, tax_group_other, '')), size=2.5, min.segment.length=0.1) +
  scale_color_manual("tax group", values = colors) +
  xlab("mean BLUP") +
  facet_wrap(~nitrogen, ncol = 2) +
  theme_bw()
  
h2_vs_blup






po_h2stdN <- po %>%
  filter(nitrogen == "+N") %>%
  select(-mean_blup, -nitrogen) %>%
  rename(H2_stdN = H2)


po_h2lowN <- po %>%
  filter(nitrogen == "-N") %>%
  select(-mean_blup, -nitrogen) %>%
  rename(H2_lowN = H2)


po_h2h2 <- po_h2stdN %>%
  left_join(po_h2lowN) #%>%
  #filter(H2_stdN > 0) %>%
  #filter(H2_lowN > 0)


h2_vs_h2 <- ggplot(po_h2h2, aes(x = H2_stdN, y = H2_lowN, color=tax_group_other)) +
  geom_point() +
  #geom_label_repel(aes(label=ifelse(tax_group_other %in% top_taxa, tax_group_other, '')), size=2.5, min.segment.length=0.1) +
  scale_color_manual("tax group", values = colors) +
  xlab("H2 (+N)") +
  ylab("H2 (-N)") +
  #facet_wrap(~nitrogen, ncol = 2) +
  theme_bw()

h2_vs_h2













