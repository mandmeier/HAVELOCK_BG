library("phyloseq")
library("tidyverse")
library("ggtree")
library("ggstance")
library("ape")


load("cache/ps_core5t.rda")


stdN <- subset_samples(ps_core5t, nitrogen == "+N")
lowN <- subset_samples(ps_core5t, nitrogen == "-N")

# get total ASV counts for +N and -N
total_counts <- rownames_to_column(data.frame("stdN" = colSums(otu_table(stdN)), "lowN" = colSums(otu_table(lowN))), var = "ASV")
total_counts$total_counts <- total_counts$stdN + total_counts$lowN
total_counts$ratioN <- total_counts$stdN/total_counts$lowN

taxa <- rownames_to_column(data.frame(tax_table(ps_core5t)), var = "ASV")


plot_data <- left_join(total_counts, taxa)




tree_plot <- function(ps, po, family){
  
  #ps <- ps_core5t
  #po <- plot_data
  #family <- "Burkholderiaceae"
  
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


burk_short <- tree_plot(ps_core5t, plot_data, "Burkholderiaceae")


### save img to ppt as vector graphic

library("officer")
library("rvg")
doc <- read_pptx()
doc <- add_slide(doc, "Title and Content", "Office Theme")
doc <- ph_with(doc, dml(ggobj = burk_short), location = ph_location_fullsize())
print(doc, target = 'temp/burk_short.pptx')









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

# remove species column to get tax_glom to work

taxtab <- taxtab %>%
  relocate(Genus, .after = tax_group)

tax_table(ps_core5t) <- as.matrix(taxtab)

save(ps_core5t, file = "cache/ps_core5t.rda")




### agglomerate ASV counts for each taxonomic group
# did this on HCC
#load("cache/ps_core5t.rda")


ps_glom <- tax_glom(ps_core5t, taxrank="tax_group")
save(ps_glom, file = "cache/ps_glom.rda")


# find number of unique ASVs for each taxonomic group

unique_asvs <- plot_data %>%
  group_by(tax_group) %>%
  tally(name = "unique_ASVs") %>%
  arrange(desc(unique_ASVs)) %>%
  filter(unique_ASVs >= 5)


## remove taxonomic groups with fewer than 5 ASVs
## this will be the final table at the group level
ps_grp <- subset_taxa(ps_glom, tax_group %in% unique_asvs$tax_group)
save(ps_grp, file = "data/ps_grp.rda")

### subset the same groups in the ASV level ps object
## this will be the final table at the ASV level
ps_asv <- subset_taxa(ps_core5t, tax_group %in% unique_asvs$tax_group)
save(ps_asv, file = "data/ps_asv.rda")






#### draw phylogenetic tree ####



### find node number between archaea and nearest bacteria
find_asvs <- data.frame(tax_table(ps_grp))
### Archaea: Candidatus Nitrocosmicus oleophilus asv_001204	
### Proteobacteria: Burkholderia pseudomallei asv_000031
getMRCA(phy_tree(ps_grp), c("asv_001204", "asv_004558", "asv_000552")) ## node number 201


### root tree
phy_tree(ps_grp) <- ape::root.phylo(phy_tree(ps_grp), node = 201)


#### plot circular tree with annotation ####


unique_asvs <- unique_asvs %>%
  arrange(tax_group)


plot_data <- data.frame(tax_table(ps_grp)) %>%
  arrange(tax_group) %>%
  select(-Genus,-Species) %>%
  rownames_to_column(var = "ASV")

plot_data$unique_ASVs <- unique_asvs$unique_ASVs


### save data for later in group data summary
group_data <- plot_data
save(group_data, file = "data/group_data.rda")


mycols <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075')

class <- tax_table(ps_grp)[,3]
#class[class == "Blastocatellia_(Subgroup_4)"] <- "Blastocatellia" ### fix that long name

p <- ggtree(phy_tree(ps_grp), layout = "circular",  branch.length = "none") +
  scale_x_continuous(expand=expand_scale(0.5))

p2 <-gheatmap(p, class, offset=3, width=0.1, font.size=3, colnames=F) +
  scale_fill_manual(values = mycols)

p3 <- p2 %<+% plot_data +
  geom_tiplab2(aes(angle = angle, label=unique_ASVs), size = 2, offset = 1) +
  geom_tiplab2(aes(angle = angle, label=tax_group), size = 2.1, offset = 7) +
  theme(plot.margin = unit(c(2,0,2,0), "cm"), legend.margin=margin(0,0,0,0, unit = "cm"), legend.box.margin=margin(0,0,0,1, unit = "cm")) +
  labs(fill="Class") 

p3


### save img to ppt as vector graphic

library("officer")
library("rvg")
doc <- read_pptx()
doc <- add_slide(doc, "Title and Content", "Office Theme")
doc <- ph_with(doc, dml(ggobj = p3), location = ph_location_fullsize())
print(doc, target = 'figures/find_core_taxa_tree.pptx')









