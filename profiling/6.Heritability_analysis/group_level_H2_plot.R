


load("data/group_data.rda")

#### plot H2 vs. abundance ####






colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")


colors <- c("Burkholderia" = "#0000bd",
  "Sphingobium" = "#f20019", "Ralstonia" = "#70fe00",
  "f__Enterobacteriaceae" = "#8f00c7", "Sphingobacterium" = "#0086fe",
  "Massilia" = "#008000", "Fulvimonas" = "#00fefe",
  "Chitinophaga" = "#fe68fe", "Labrys" = "#fe8420",
  "Chryseobacterium" = "#827800", "f__Rhodanobacteraceae" = "#fefe00",
  "other" = "#eeeeee")




### h2 vs BLUP



po <- group_data %>%
  filter(H2_stdN > 0)

h2_vs_count_stdN <- ggplot(po, aes(x = count_stdN, y = H2_stdN)) +
  geom_point() #+
#scale_color_manual(values = colors) +
theme_bw()#+
#facet_wrap(~genus2, nrow = 3)

h2_vs_count_stdN



po <- group_data %>%
  filter(H2_lowN > 0)

h2_vs_count_lowN <- ggplot(po, aes(x = count_lowN, y = H2_lowN)) +
  geom_point() #+
#scale_color_manual(values = colors) +
theme_bw()#+
#facet_wrap(~genus2, nrow = 3)

h2_vs_count_lowN



po <- group_data %>%
  filter(H2_lowN > 0 & H2_stdN > 0)


h2_vs_h2 <- ggplot(po, aes(x = H2_stdN, y = H2_lowN)) +
  geom_point() #+
#scale_color_manual(values = colors) +
theme_bw()#+
#facet_wrap(~genus2, nrow = 3)

h2_vs_h2





