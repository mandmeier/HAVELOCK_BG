#plot_blup

library("tidyverse")
library("phyloseq")

load("data/group_data.rda")

blups <- read_csv("data/blup_stdN_150_tax_groups.csv")
blups <- gather(blups, key = "ASV", value = "counts", -genotype)
blups$nitrogen <- "+N"



blupl <- read_csv("data/blup_lowN_150_tax_groups.csv")
blupl <- gather(blupl, key = "ASV", value = "counts", -genotype)
blupl$nitrogen <- "-N"

blup <- rbind(blups, blupl)
po <- left_join(blup, group_data)



### get tree tip order


## helper function to find tree tip order
tip_order <- function(tree){
  d <- fortify(tree)
  d <- subset(d, isTip)
  ord <- with(d, label[order(y, decreasing=T)])
  return(ord)
}


load("data/ps_grp.rda")

ps_grp

phy_tree(ps_grp)

grp_order <- data.frame("ASV"=tip_order(phy_tree(ps_grp)))

plot_order <- left_join(grp_order, group_data)

plot_order$tax_group


po$tax_group <- factor(po$tax_group, levels = plot_order$tax_group)


blup_plot <- ggplot(po, aes(x=tax_group, y= counts, fill = nitrogen)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#ff0000","#0000ff")) +
  theme_bw() +
  ylab("BLUP log(relative abundance)") +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

blup_plot




?gather
