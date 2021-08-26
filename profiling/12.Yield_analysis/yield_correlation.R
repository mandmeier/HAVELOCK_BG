# correlation of microbe abundance with yield traits


library("tidyverse")
library("ggpubr")



## summarize abundance BLUPS and yield traits for each genotype

load("data/data_summary_150_traits.rda")

test <- read_csv("data/abundance_logrel.csv")


data_blup_LN <- read_delim(file="data/Blup_and_Heritability/data_blup_LN_result.txt", delim = "\t")

data_blup_LN <- data_blup_LN %>%
  drop_na() %>%
  pivot_longer(cols=c(-id),names_to="trait", values_to = "abundance") %>%
  rename(GX_name = id, ASV = trait, blup_logrel = abundance) %>%
  mutate(nitrogen="LN")




data_blup_HN <- read_delim(file="data/Blup_and_Heritability/data_blup_HN_result.txt", delim = "\t")

data_blup_HN <- data_blup_HN %>%
  drop_na() %>%
  pivot_longer(cols=c(-id),names_to="trait", values_to = "abundance") %>%
  rename(GX_name = id, ASV = trait, blup_logrel = abundance) %>%
  mutate(nitrogen="HN")

names <- read_csv("data/BG_MM_Gen_names.csv")

blup_data <- rbind(data_blup_LN, data_blup_HN)
blup_data <- blup_data %>%
  left_join(dplyr::select(names, MM_name, GX_name, CS_name))

#save(blup_data, file="cache/BLUP/blup_data.rda")


## Semra yield component traits, Eric imaging traits
load("data/yield_analysis/microbe_counts_and_yield.rda")
load("data/group_data.rda")

yield_traits <- microbe_counts_and_yield %>%
  rename(MM_name=genotype) %>%
  mutate(nitrogen = ifelse(nitrogen == "+N", "HN", "LN")) %>%
  left_join(dplyr::select(group_data, ASV, tax_group))



## yufeng nutrient traits
nutrient_traits <- read_delim(file="data/yield_analysis/phenotype.txt", delim="\t")
nutrient_traits_HN <- nutrient_traits %>%
  dplyr::select(ID, starts_with("HN")) %>%
  mutate(nitrogen="HN")
colnames(nutrient_traits_HN) <- c("GX_name", "B", "Ca", "CHL", "Cu", "DW", "Fe", "FW", "K", "LA", "Mg", "Mn", "N", "P", "S", "Zn", "nitrogen")

nutrient_traits_LN <- nutrient_traits %>%
  dplyr::select(ID, starts_with("LN")) %>%
  mutate(nitrogen="LN")

colnames(nutrient_traits_LN) <- c("GX_name", "B", "Ca", "CHL", "Cu", "DW", "Fe", "FW", "K", "LA", "Mg", "Mn", "N", "P", "S", "Zn", "nitrogen")

nutrient_traits <- rbind(nutrient_traits_LN, nutrient_traits_HN)

## combine data

abundance_vs_phenotype <- blup_data %>%
  left_join(yield_traits) %>%
  left_join(nutrient_traits) %>%
  pivot_longer(cob_length:Zn, names_to = "phenotype", values_to = "value") %>%
  select(nitrogen, GX_name, MM_name, CS_name, tax_group, ASV, subpopulation, count, total_count, relab, log_relab, blup_logrel, phenotype, value )


#save(abundance_vs_phenotype, file="data/yield_analysis/abundance_vs_phenotype.rda")




#### for each microbe and phenotype calculate correlation between abundance and phenotype


corr_data <- abundance_vs_phenotype %>%
  dplyr::select(nitrogen, tax_group, phenotype) %>%
  unique()

corr_data$pearson <- 0
corr_data$p_value <- 0


get_correlation <- function(n, m, p, var = "blup_logrel"){
  #n <- "LN"
  #m <- "Burkholderia sp 2"
  #p <- "cob_length"
  #var <- "blup_logrel"
  #p <- corr_data[i,]$phenotype
  nmp <- data.frame(filter(abundance_vs_phenotype, nitrogen == n & tax_group == m & phenotype == p))
  cor <- cor.test(nmp[, var], nmp[, "value"], method=c("pearson"))
  return(cor)
}


for (i in c(1:nrow(corr_data))){
  print(i)
  cor <- get_correlation(corr_data[i,]$nitrogen, corr_data[i,]$tax_group, corr_data[i,]$phenotype)
  pears <- cor$estimate[[1]]
  p_value <- cor$p.value
  corr_data[i,]$pearson <- pears
  corr_data[i,]$p_value <- p_value
}


#save(corr_data, file="data/yield_analysis/corr_data.rda")

load("data/yield_analysis/corr_data.rda")


ggplot(corr_data, aes(x=pearson, y=-log10(p_value), color=nitrogen)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 2, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#000080", "#C00001")) +
  facet_wrap(~phenotype+nitrogen, nrow = 4) +
  xlab("pearson correlation") +
  theme_bw()



CC <- corr_data %>%
  filter(phenotype == "CC_Aug12")


ggplot(CC, aes(x=pearson, y=-log10(p_value))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 2, color = "black", linetype = "dashed") +
  #scale_color_manual(values = c("#000080", "#C00001")) +
  facet_wrap(~nitrogen, nrow = 1) +
  xlab("pearson correlation") +
  theme_bw()



shortlist <- corr_data %>%
  filter(p_value <= 0.01)



sort(unique(shortlist$tax_group))



load("data/yield_analysis/abundance_vs_phenotype.rda")






df <- abundance_vs_phenotype %>%
  filter(nitrogen == "LN" & tax_group == "Steroidobacter") %>%
  filter(phenotype %in% c("CC_Aug12", "ExG_Aug12", "LA", "FW", "S"))
  #filter(phenotype %in% c("CC_Aug12", "S"))



df <- abundance_vs_phenotype %>%
  filter(nitrogen == "LN" & tax_group == "Conexibacter") %>%
  filter(phenotype %in% c("CC_Aug12", "ExG_Aug12", "LA", "FW", "S"))
  #filter(phenotype %in% c("CC_Aug12", "LA", "FW", "S"))




df <- abundance_vs_phenotype %>%
  filter(nitrogen == "HN" & tax_group == "Acidothermus") %>%
  filter(phenotype %in% c("CC_Aug12", "ExG_Aug12", "LA", "FW", "S"))
  #filter(phenotype %in% c("CC_Aug12", "ExG_Aug12", "LA", "S"))


df <- abundance_vs_phenotype %>%
  filter(nitrogen == "HN" & tax_group == "Sphingoaurantiacus") %>%
  filter(phenotype %in% c("CC_Aug12", "ExG_Aug12", "LA", "FW", "S"))
  #filter(phenotype %in% c("CC_Aug12", "ExG_Aug12"))






df <- abundance_vs_phenotype %>%
  filter(nitrogen == "HN" & tax_group == "Sphingoaurantiacus") %>%
  filter(phenotype %in% c("CC_Aug12"))
#filter(phenotype %in% c("CC_Aug12", "ExG_Aug12"))


df <- abundance_vs_phenotype %>%
  filter(nitrogen == "LN" & tax_group == "Acinetobacter nosocomialis") %>%
  filter(phenotype %in% c("CC_Aug12"))
#filter(phenotype %in% c("CC_Aug12", "ExG_Aug12"))




df <- abundance_vs_phenotype %>%
  filter(nitrogen == "LN" & tax_group == "Steroidobacter") %>%
  filter(phenotype %in% c("CC_Aug12"))
#filter(phenotype %in% c("CC_Aug12", "S"))



df <- abundance_vs_phenotype %>%
  filter(nitrogen == "LN" & tax_group == "Conexibacter") %>%
  filter(phenotype %in% c("CC_Aug12"))
#filter(phenotype %in% c("CC_Aug12", "LA", "FW", "S"))




var <- "blup_logrel"


#### plot correlation
ggscatter(df, x = var, y = "value", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "abundance", ylab = var, alpha=0.25) +
  ylab("Canopy Coverage") +
  #facet_wrap(~phenotype, scales = "free_y", nrow = 1) +
  theme_bw()




#### subset microbes that positively correleate with canopy cover


CC_shortlist <- corr_data %>%
  filter(phenotype == "CC_Aug12")






### heritability vs p_value

H2_vs_CC <- CC_shortlist %>%
  left_join(group_data) %>%
  #dplyr::select(nitrogen, tax_group, pearson, p_value, H2_stdN_19, H2_lowN_19) %>%
  mutate(log10p=-log10(p_value)) %>%
  mutate(heritability = ifelse(nitrogen == "HN", H2_stdN_19, H2_lowN_19)) %>%
  mutate(abundance = ifelse(nitrogen == "HN", mean_blup_stdN, mean_blup_lowN)) 



#### plot correlation
ggscatter(H2_vs_CC, x = "heritability", y = "log10p",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "heritability", ylab = "correlation with CC [-log10(p)]", alpha=0.5) +
  scale_color_manual(values = c("#000080", "#C00001")) +
  geom_hline(yintercept = 2, color = "black", linetype = "dashed") +
  facet_wrap(~nitrogen, nrow = 1) +
  theme_bw()

#### plot correlation
ggscatter(H2_vs_CC, x = "abundance", y = "log10p",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "abundance", ylab = "correlation with CC [-log10(p)]", alpha=0.5) +
  scale_color_manual(values = c("#000080", "#C00001")) +
  geom_hline(yintercept = 2, color = "black", linetype = "dashed") +
  facet_wrap(~nitrogen, nrow = 1) +
  theme_bw()






#### plot correlation
ggscatter(H2_vs_CC, x = "H2_stdN_19", y = "log10p",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "heritability", ylab = "correlation with CC [-log10(p)]", alpha=0.5) +
  scale_color_manual(values = c("#000080", "#C00001")) +
  facet_wrap(~nitrogen, nrow = 1) +
  theme_bw()




H2_vs_CC_shortlist <- H2_vs_CC %>%
  filter(p_value <= 0.01)

  



total <- sort(unique(H2_vs_CC_shortlist$tax_group))
taxa_HN <- unique(filter(H2_vs_CC_shortlist, nitrogen == "HN")$tax_group)
taxa_LN <- unique(filter(H2_vs_CC_shortlist, nitrogen == "LN")$tax_group)
common <- sort(Reduce(intersect, list(taxa_HN,taxa_LN)))
only_HN <- sort(taxa_HN[!(taxa_HN %in% common)])
only_LN <- sort(taxa_LN[!(taxa_LN %in% common)])


### find CC associated taxa
CC_associated_taxa <- filter(H2_vs_CC_shortlist, tax_group %in% common & nitrogen == "LN")


### check if PLAMs
load("cache/GWAS2/v5_MAPLs.rda")


PLAMs <- sort(unique(v5_MAPLs$tax_group))


CC_PLAMs <- total[total %in% PLAMs]

total[!(total %in% PLAMs)]

### check heritability

H2_vs_CC_shortlist$sign_lowN <- ifelse(H2_vs_CC_shortlist$perm_p_lowN < 0.05, "sign", "n.s.")
H2_vs_CC_shortlist$sign_stdN <- ifelse(H2_vs_CC_shortlist$perm_p_stdN < 0.05, "sign", "n.s.")

H2_vs_CC_shortlist$sign_group <- "none"
H2_vs_CC_shortlist$sign_group <- ifelse(H2_vs_CC_shortlist$sign_lowN == "sign", "lowN", H2_vs_CC_shortlist$sign_group)
H2_vs_CC_shortlist$sign_group <- ifelse(H2_vs_CC_shortlist$sign_stdN == "sign", "stdN", H2_vs_CC_shortlist$sign_group)
H2_vs_CC_shortlist$sign_group <- ifelse(H2_vs_CC_shortlist$sign_stdN == "sign" & H2_vs_CC_shortlist$sign_lowN == "sign", "both", H2_vs_CC_shortlist$sign_group)
H2_vs_CC_shortlist$sign_group <- factor(H2_vs_CC_shortlist$sign_group, levels = c("stdN", "lowN", "both", "none"))



CC_associated_taxa <- filter(H2_vs_CC_shortlist, tax_group %in% only_LN) %>%
  dplyr::select(nitrogen, tax_group, sign_group)



### check for strong GWAS signals

gwas_sig <- read_csv("data/microbe_references.csv")


CC_associated_taxa <- filter(H2_vs_CC_shortlist, tax_group %in% common & nitrogen == "HN") %>%
  left_join(gwas_sig) %>%
  dplyr::select(tax_group, gwas_wow)




