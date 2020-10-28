library("phyloseq")
library("tidyverse")
library("ggtree")
library("ggstance")


load("cache/ps_core5t.rda")


stdN <- subset_samples(ps_core5t, nitrogen == "+N")
lowN <- subset_samples(ps_core5t, nitrogen == "-N")

# get total ASV counts for +N and -N
total_counts <- rownames_to_column(data.frame("stdN" = colSums(otu_table(stdN)), "lowN" = colSums(otu_table(lowN))), var = "ASV")
total_counts$total_counts <- total_counts$stdN + total_counts$lowN
total_counts$ratioN <- total_counts$stdN/total_counts$lowN

taxa <- data.frame(tax_table(ps_core5t))


plot_data <- left_join(total_counts, taxa)




tree_plot <- function(ps, po, family){
  
  #ps <- ps_core5t
  #po <- plot_data
  #family <- "Nitrosomonadaceae"
  
  ### Subset family
  ps_fam <- subset_taxa(ps, Family == family)
  po_fam <- subset(plot_data, Family == family)
  
  colnames(po_fam)[1] <- "label"
  
  p1 <- ggtree(phy_tree(ps_fam), branch.length = "none") +
    geom_tiplab(align=TRUE, size = 2) #+
  #scale_x_continuous(expand=expand_scale(0.8)) #+
  #theme(plot.margin = unit(c(0,2,0,0), "cm"))
  
  p2 <- facet_plot(p1, panel="log10(total ASV counts)", data=po_fam, geom=geom_point, aes(x=log10(total_counts)), color="firebrick") + theme_tree2()
  
  p3 <- facet_plot(p2, panel="count ratio log2(+N/-N)", data=po_fam, geom=geom_point, aes(x=log2(ratioN)), color="blue") + theme_tree2()
  
  p4 <- p3 %<+% po_fam +
    geom_tiplab(aes(label=Genus), align=T, linetype=NA, size=2, offset=8, hjust=0.5) +
    geom_tiplab(aes(label=Species), align=T, linetype=NA, size=2, offset=14, hjust=0.5) +
    ggtitle(family)
  
  p4
  
  return(p4)
  
}

fam



#### draw plot for all families, find subgroups manually ####

## helper function to find tree tip order
tip_order <- function(tree){
  d <- fortify(tree)
  d <- subset(d, isTip)
  ord <- with(d, label[order(y, decreasing=T)])
  return(ord)
}

## list 58 families in frequency table
fam <- arrange(plyr::count(plot_data$Family), desc(freq))

asv_order <- c()
for (f in fam$x){
  ### Subset family
  ps_fam <- subset_taxa(ps_core5t, Family == f)
  asv_order <- c(asv_order, tip_order(phy_tree(ps_fam)))
}

asv_order <- data_frame("ASV"=asv_order)

plot_data_ord <- left_join(asv_order, plot_data)

### use this csv file to manually assign best matching taxonomic group to each ASV
write_csv(plot_data_ord, file = "cache/tax_groups_man.csv")

tree_plot(ps_core5t, plot_data, "Chitinophagaceae")

### load file with manually curated taxonomic groups
tax_groups_man <- read_csv("cache/tax_groups_man.csv")

unique(tax_groups_man$man_taxa)


### add manual taxonomic groups to ps object

tax_groups <- select(tax_groups_man, ASV, tax_group)

taxtab <- as.data.frame(tax_table(ps_core5t))

taxtab <- taxtab %>%
  left_join(tax_groups) %>%
  select(-taxa) ## remove taxa column, no longer needed
rownames(taxtab) <- taxtab$ASV

tax_table(ps_core5t) <- as.matrix(taxtab)

save(ps_core5t, file = "cache/ps_core5t.rda")



library("phyloseq")
library("tidyverse")

load("cache/ps_core5t.rda")

tax_table(ps_core5t)


### agglomerate ASV counts for each taxonomic group

ps_glom <- tax_glom(ps_core5t, taxrank="tax_group")

save(ps_glom, file = "cache/ps_glom.rda")



#### draw phylogenetic tree ####






### agglomerate taxa at subgroup level





ps_glom <- tax_glom(ps_final_S, taxrank=rank_names(ps_final_S)[7])








