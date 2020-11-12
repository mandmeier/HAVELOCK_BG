library("tidyverse")
library("phyloseq")

load("cache/asv_H2_lowN.rda")
load("cache/asv_H2_stdN.rda")


blup_stdN <- read_csv("data/blup_stdN_3618_asvs.csv")
mean_blup_stdN <- colMeans(blup_stdN[-1], na.rm = TRUE)
mean_blup_stdN <- data.frame("ASV"=names(mean_blup_stdN), "mean_blup" = mean_blup_stdN)
mean_blup_stdN$nitrogen <- "+N"
mean_blup_stdN <- mean_blup_stdN %>%
  left_join(H2_stdN) %>%
  rename(H2 = H2_stdN)


blup_lowN <- read_csv("data/blup_lowN_3618_asvs.csv")
mean_blup_lowN <- colMeans(blup_lowN[-1], na.rm = TRUE)
mean_blup_lowN <- data.frame("ASV"=names(mean_blup_lowN), "mean_blup" = mean_blup_lowN)
mean_blup_lowN$nitrogen <- "-N"
mean_blup_lowN <- mean_blup_lowN %>%
  left_join(H2_lowN) %>%
  rename(H2 = H2_lowN)


mean_blup <- rbind(mean_blup_stdN, mean_blup_lowN)

### add taxonomy

load("data/ps_asv.rda")

taxtab <- data.frame(tax_table(ps_asv)) %>%
  rownames_to_column(var = "ASV")


asv_data <- left_join(taxtab, mean_blup)

save(asv_data, file = "data/asv_data.rda")


load("data/asv_data.rda")

## #get top heritable taxa determined at group level
load("cache/top_heritable_taxa_150.rda")

## assign colors 
asv_data$tax_group_other <- ifelse(asv_data$tax_group %in% top_taxa, asv_data$tax_group, "other")

colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")



top_taxa_other <- sort(unique(asv_data$tax_group_other))
group_colors <- data.frame("tax_group_other" = top_taxa_other, "group_color" = colors37[1:length(top_taxa_other)])
## assign grey to "other"
group_colors$group_color[group_colors$tax_group_other == "other"] <- "#eeeeee"
group_colors$group_color[group_colors$tax_group_other == "Sphingobium herbicidovorans 1"] <- "#ff0000"
group_colors$group_color[group_colors$tax_group_other == "Sphingobium herbicidovorans 2"] <- "#800000"
group_colors$group_color[group_colors$tax_group_other == "Niabella yanshanensis "] <- "#ffa500"

asv_data <- left_join(asv_data, group_colors, by = "tax_group_other")

colors <- asv_data$group_color
names(colors) <-asv_data$tax_group_other

colors <- colors[unique(names(colors))]

save(colors, file = "cache/taxa_colors.rda")


po <- asv_data %>%
  filter(H2 > 0) %>%
  filter(mean_blup > median(asv_data$mean_blup))

h2_vs_blup <- ggplot(po, aes(x = mean_blup, y = H2, color=tax_group_other)) +
  geom_point() +
  #geom_label_repel(aes(label=ifelse(tax_group_other %in% top_taxa, tax_group_other, '')), size=2.5, min.segment.length=0.1) +
  scale_color_manual("tax group", values = colors) +
  xlab("mean BLUP") +
  facet_wrap(~nitrogen, ncol = 2) +
  theme_bw()

h2_vs_blup



po_h2stdN <- asv_data %>%
  filter(nitrogen == "+N") %>%
  select(-mean_blup, -nitrogen) %>%
  rename(H2_stdN = H2)


po_h2lowN <- asv_data %>%
  filter(nitrogen == "-N") %>%
  select(-mean_blup, -nitrogen) %>%
  rename(H2_lowN = H2)
  




po_h2h2 <- po_h2stdN %>%
  left_join(po_h2lowN) %>%
  filter(H2_stdN > 0) %>%
  filter(H2_lowN > 0)


h2_vs_h2 <- ggplot(po_h2h2, aes(x = H2_stdN, y = H2_lowN, color=tax_group_other)) +
  geom_point() +
  #geom_label_repel(aes(label=ifelse(tax_group_other %in% top_taxa, tax_group_other, '')), size=2.5, min.segment.length=0.1) +
  scale_color_manual("tax group", values = colors) +
  xlab("H2 (+N)") +
  ylab("H2 (-N)") +
  #facet_wrap(~nitrogen, ncol = 2) +
  theme_bw()

h2_vs_h2





