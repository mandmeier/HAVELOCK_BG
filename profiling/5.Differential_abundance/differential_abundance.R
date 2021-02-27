### Differential Abundance by genotype
library("tidyverse")
library("phyloseq")
library("tidytext")

load("data/data_summary_150_traits.rda")

load("data/ps_grp.rda")

load("data/group_data.rda")

## get absolute counts and relative abundance for each genotype
### get abundance by genotype

counts <- data.frame(otu_table(ps_grp))
counts <- rownames_to_column(counts, var="Sample_ID")

sdat <- data.frame(sample_data(ps_grp))


count_data <- sdat %>%
  filter(genotype != "CHECK") %>%
  left_join(counts) %>%
  pivot_longer(cols = starts_with("asv_"), names_to="ASV", values_to="count") %>%
  left_join(group_data[, c("ASV", "tax_group")]) %>%
  group_by(genotype, tax_group, nitrogen, subpopulation) %>%
  summarise(count = sum(count)) %>%
  group_by(genotype, nitrogen) %>%
  mutate(total_count = sum(count)) %>%
  mutate( relab = count/total_count)



#save(count_data, file = "cache/counts_per_genotype.rda")

### filter top22 tax groups

load("data/data_summary_150_traits.rda")

top_taxa <- data_summary_150_traits %>%
  filter(top_stdN == "x" | top_lowN == "x")

top_taxa <- unique(top_taxa$tax_group)


count_data_top <- filter(count_data, tax_group %in% top_taxa)




### violin plot to find bimodal distributions of counts fro each microbe


count_data_top <- filter(count_data, tax_group %in% top_taxa)


viol <- ggplot(count_data_top, aes(x=nitrogen, y=count, color=nitrogen)) +
  geom_violin() +
  scale_color_manual(values = c("red", "blue")) +
  facet_wrap(~tax_group, ncol = 5, scales = "free_y") +
  theme_bw()

viol



count_data_lesser <- filter(count_data, !(tax_group %in% top_taxa))


lesser_1 <- unique(count_data_lesser$tax_group)[c(1:25)]

lesser_2 <- unique(count_data_lesser$tax_group)[c(26:50)]



count_data_lesser1 <- filter(count_data, tax_group %in% lesser_1)

viol <- ggplot(count_data_lesser1, aes(x=nitrogen, y=count, color=nitrogen)) +
  geom_violin() +
  scale_color_manual(values = c("red", "blue")) +
  facet_wrap(~tax_group, ncol = 5, scales = "free_y") +
  theme_bw()

viol




count_data_lesser2 <- filter(count_data, tax_group %in% lesser_2)

viol <- ggplot(count_data_lesser2, aes(x=nitrogen, y=count, color=nitrogen)) +
  geom_violin() +
  scale_color_manual(values = c("red", "blue")) +
  facet_wrap(~tax_group, ncol = 5, scales = "free_y") +
  theme_bw()

viol





relab_data_top <- count_data_top %>%
  mutate(nitrogen = ifelse(nitrogen == "+N", "stdN", "lowN")) %>%
  pivot_wider(names_from = nitrogen, values_from = count) %>%
  mutate(total_count = stdN + lowN) %>%
  mutate(stdN_norm = stdN/total_count) %>%
  mutate(lowN_norm = -lowN/total_count) %>%
  pivot_longer(cols= ends_with("norm"), names_to="nitrogen", values_to="abundance")


## move NA to "mixed", i.e. "others"
relab_data_top$subpopulation <- as.character(relab_data_top$subpopulation)
relab_data_top$subpopulation <- ifelse(is.na(relab_data_top$subpopulation), "mixed", relab_data_top$subpopulation)




?pivot_longer


addColor <- function(x) {
  if (x =="nss"){
    return("#00fefe")
  } else if (x == "popcorn") {
    return("#70fe00")
  } else if (x == "mixed") {
    return("#fe68fe")
  } else if (x == "ss") {
    return("#fed38b")
  } else if (x == "ts") {
    return("#fefe00")
  } else {
    return("#0086fe")
  }
}

taxa <- relab_data_top %>%
  #filter(tax_group %in% c("Acinetobacter nosocomialis", "Massilia niabensis", "Rhizobium daejeonense")) %>%
  #filter(tax_group %in% c("Acinetobacter nosocomialis")) %>%
  #filter(tax_group %in% c("Massilia niabensis")) %>%
  #filter(tax_group %in% c("Rhizobium daejeonense")) %>%
  filter(tax_group %in% c("Candidatus Udaeobacter copiosus")) %>%
  mutate(subpopColor = addColor(subpopulation))
  #mutate(tax_group = as.factor(tax_group), genotype = reorder_within(genotype, abundance, tax_group))
  #mutate(genotype =  str_split(genotype, "___")[[1]][1])
  

sort(unique(relab_data_top$tax_group))

test <- Rd %>%
  filter(nitrogen == "stdN_norm") %>%
  select(-nitrogen, -subpopulation) #%>%
  #pivot_longer(cols = ends_with("N"), names_to = "nitrogen", values_to = "count")



p <-  ggplot(test, aes(x=reorder(genotype, -lowN), y = lowN, fill = "red")) +
      geom_bar(stat = "identity", position="stack") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
    
p



p <-  ggplot(test, aes(x=reorder(genotype, -total_count), y = count, fill = nitrogen)) +
        geom_bar(stat = "identity", position="stack") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))


p





unique(taxa$subpopulation)


subpops <- c("#00fefe", "#70fe00", "#fe68fe", "#fed38b", "#fefe00", "#0086fe" )

names(subpops) <- unique(taxa$subpopulation)




subpops <- c("#00fefe", "#70fe00", "#fe68fe", "#fed38b", "#fefe00", "#0086fe" )

names(subpops) <- unique(taxa$subpopulation)

p <- ggplot(taxa) +
  geom_bar(aes(x=reorder(genotype, -abundance), y=abs(abundance), fill=nitrogen), stat='identity', position ="stack") +
  #facet_wrap(~tax_group, ncol = 1, scales = "free") +
  #scale_x_reordered() +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, fill=taxa$subpopColor))

p


p1 <- ggplot(taxa) +
  geom_bar(aes(x=reorder(genotype, -abundance), y=total_count, fill=nitrogen), stat='identity', position ="stack") +
  #facet_wrap(~tax_group, ncol = 1, scales = "free") +
  #scale_x_reordered() +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, color=taxa$subpopColor))

p1








taxa <- relab_data_top %>%
  filter(tax_group %in% c("Acinetobacter nosocomialis")) %>%
  mutate(subpopColor = addColor(subpopulation))


p <- ggplot(taxa) +
  geom_bar(aes(x=reorder(genotype, -abundance), y=abundance, fill=nitrogen), stat='identity', position = position_dodge(width = 0)) +
  #facet_wrap(~tax_group, ncol = 1, scales = "free") +
  #scale_x_reordered() +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

p



taxa <- relab_data_top %>%
  filter(tax_group %in% c("Niabella yanshanensis")) %>%
  mutate(subpopColor = addColor(subpopulation))




p2 <- ggplot(taxa) +
  geom_bar(aes(x=reorder(genotype, -abundance), y=abundance, fill=nitrogen), stat='identity', position = position_dodge(width = 0)) +
  #facet_wrap(~tax_group, ncol = 1, scales = "free") +
  #scale_x_reordered() +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

p2




taxa <- relab_data_top %>%
  filter(tax_group %in% c("Candidatus Udaeobacter copiosus")) %>%
  mutate(subpopColor = addColor(subpopulation))



p3 <- ggplot(taxa) +
  geom_bar(aes(x=reorder(genotype, -abundance), y=abundance, fill=nitrogen), stat='identity', position = position_dodge(width = 0)) +
  #facet_wrap(~tax_group, ncol = 1, scales = "free") +
  #scale_x_reordered() +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

p3



taxa <- relab_data_top %>%
  filter(tax_group %in% c("Rhizobium daejeonense")) %>%
  mutate(subpopColor = addColor(subpopulation))




p4 <- ggplot(taxa) +
  geom_bar(aes(x=reorder(genotype, -abundance), y=abundance, fill=nitrogen), stat='identity', position = position_dodge(width = 0)) +
  #facet_wrap(~tax_group, ncol = 1, scales = "free") +
  #scale_x_reordered() +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

p4






library(scales)
numColors <- length(levels(melted_state$Region)) # How many colors you need
getColors <- scales::brewer_pal('qual') # Create a function that takes a number and returns a qualitative palette of that length (from the scales package)
myPalette <- getColors(numColors)
names(myPalette) <- levels(state_data$Region) # Give every color an appropriate name
p <- p + theme(axis.text.x = element_text(colour=myPalette[state_data$Region])))






ordered <- taxa %>%
  ungroup() %>%
  arrange(tax_group, abundance) %>%
  mutate(order = row_number())






p <- ggplot(ordered) +
  geom_bar(aes(x=order, y=abundance, fill=nitrogen), stat='identity', position = position_dodge(width = 0)) +
  coord_flip() +
  facet_wrap(~subpopulation, nrow = 2, scales = "free")
  
p






  
  geom_bar(aes(x=factor(genotype), y=Amount, fill=Type), stat='identity', position=position_dodge(width = 0)) +
  coord_flip() +
  ylim(c(-100,100))










