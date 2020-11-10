library("tidyverse")


load("data/group_data.rda")

#### plot H2 vs. BLUP abundance ####


colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")


colors <- c("Burkholderia" = "#0000bd",
  "Sphingobium" = "#f20019", "Ralstonia" = "#70fe00",
  "f__Enterobacteriaceae" = "#8f00c7", "Sphingobacterium" = "#0086fe",
  "Massilia" = "#008000", "Fulvimonas" = "#00fefe",
  "Chitinophaga" = "#fe68fe", "Labrys" = "#fe8420",
  "Chryseobacterium" = "#827800", "f__Rhodanobacteraceae" = "#fefe00",
  "other" = "#eeeeee")




### h2 vs BLUP


## mark anything in 75% quantile for H2 and BLUP

### median H2
q75_H2_stdN <- quantile(group_data$H2_stdN[group_data$H2_stdN!=min(group_data$H2_stdN)], 0.75)
### median blup
q75_blup_stdN <- quantile(group_data$mean_blup_stdN[group_data$mean_blup_stdN!=min(group_data$mean_blup_stdN)], 0.75)


q75_stdN <- group_data %>%
  filter(H2_stdN >= q75_H2_stdN | mean_blup_stdN >= q75_blup_stdN)
 

unique(q75_stdN$tax_group)


## mark anything in 75% quantile for H2 or BLUP

### median H2
q90_H2_stdN <- quantile(group_data$H2_stdN[group_data$H2_stdN!=min(group_data$H2_stdN)], 0.90)
### median blup
q90_blup_stdN <- quantile(group_data$mean_blup_stdN[group_data$mean_blup_stdN!=min(group_data$mean_blup_stdN)], 0.90)

q90_stdN <- group_data %>%
  filter(H2_stdN >= q90_H2_stdN | mean_blup_stdN >= q90_blup_stdN)


unique(q90_stdN$tax_group)





min(group_data$mean_blup_stdN)


median_H2_lowN <- median(group_data$H2_lowN[group_data$H2_lowN>0])




po <- group_data %>%
  filter(H2_stdN > 0)

h2_vs_blup_stdN <- ggplot(po, aes(x = mean_blup_stdN, y = H2_stdN)) +
  geom_point() #+
#scale_color_manual(values = colors) +
theme_bw()#+
h2_vs_blup_stdN



po <- group_data %>%
  filter(H2_lowN > 0)

h2_vs_blup_lowN <- ggplot(po, aes(x = mean_blup_lowN, y = H2_lowN)) +
  geom_point() #+
#scale_color_manual(values = colors) +
theme_bw()#+
h2_vs_blup_lowN







po <- group_data %>%
  filter(H2_lowN > 0 & H2_stdN > 0)


h2_vs_h2 <- ggplot(po, aes(x = H2_stdN, y = H2_lowN)) +
  geom_point() #+
#scale_color_manual(values = colors) +
theme_bw()#+
#facet_wrap(~genus2, nrow = 3)

h2_vs_h2





