library("phyloseq")
library("tidyverse")
library("lme4")

#### prepare data to calculate H2 ####

load("data/ps_asv.rda")

# find genotypes common to both years, without check
#sdat <- data.frame(sample_data(ps_asv))
#y2018 <- sdat %>%
#  filter(year == "Y2018")
#y2019 <- sdat %>%
#  filter(year == "Y2019")
#common_genotypes <- Reduce(intersect, list(unique(y2018$genotype),unique(y2019$genotype)))
#common_genotypes <- common_genotypes[-1]


### filter ps object with common genotypes
#ps_common <- subset_samples(ps_asv, genotype %in% common_genotypes)
#asv_table <- rownames_to_column(data.frame(otu_table(ps_common)), var = "Sample_ID")
#counts <- left_join(data.frame(sample_data(ps_common)), asv_table)


### calculate means of ASV counts of all subsamples by group
#mean_counts <- counts %>%
#  group_by(year, genotype, nitrogen, block, sp, spb) %>%
#  summarize_each(funs(mean), -Sample_ID, -raw_seq_count, -filename, -subsample, -collected_by_person, -pedigree, -subpopulation)


### calculate log relative abundance
#asvtab <- mean_counts[, 8:ncol(mean_counts)]
#logrel <- t(apply(asvtab, 1, function(x) log(x/sum(x) + 0.001)))
#mean_counts_logrel <- cbind(data.frame(mean_counts[, 1:7]), logrel)

#print("saving mean_counts_logrel")
#save(mean_counts_logrel, file = "cache/mean_counts_logrel_asv.rda")

load("cache/mean_counts_logrel_asv.rda")



### split stdN, lowN samples, retain only genotypes with complete set of 4 reps

#asv_h2dat_stdN <- mean_counts_logrel %>%
#  filter(nitrogen == "+N") %>%
#  group_by(genotype) %>%
#  dplyr::add_tally(name="count") %>%
#  filter(count == 4) %>%
#  dplyr::select(-count)
# 20 genotypes


#save(asv_h2dat_stdN, file="cache/asv_h2dat_stdN.rda")


#asv_h2dat_lowN <- mean_counts_logrel %>%
#  filter(nitrogen == "-N") %>%
#  group_by(genotype) %>%
#  dplyr::add_tally(name="count") %>%
#  filter(count == 4) %>%
#  dplyr::select(-count)
# 19 genotypes

save(asv_h2dat_lowN, file="cache/asv_h2dat_lowN.rda")



#### calculate H2 ####

load("cache/asv_h2dat_stdN.rda")

## calculate H2 +N
df <- asv_h2dat_stdN
traits <- colnames(df)[8:ncol(df)]

H2_stdN <- data.frame()
for(trait in traits){
  #trait <- "asv_000013"
  f <- formula(paste0(trait, ' ~ (1|genotype) + (1|year) + (1|block) + (1|genotype:year)'))
  fit <- lmer(f, data=df)
  v <- as.data.frame(VarCorr(fit))
  Vg <- v$vcov[v$grp == "genotype"]
  Vgxy <- v$vcov[v$grp == "genotype:year"]
  Ve <- v$vcov[v$grp == "Residual"]
  H2 <- round(Vg/(Vg + Vgxy/2 + Ve/4), 6)
  H2_stdN <- rbind(H2_stdN, data.frame("ASV"=trait, "H2_stdN"=H2))
}

H2_stdN

print("saving H2_stdN")
save(H2_stdN, file = "cache/asv_H2_stdN.rda")




## calculate H2 -N

load("cache/asv_h2dat_lowN.rda")

df <- asv_h2dat_lowN
traits <- colnames(df)[8:ncol(df)]

H2_lowN <- data.frame()
for(trait in traits){
  #trait <- "asv_000013"
  f <- formula(paste0(trait, ' ~ (1|genotype) + (1|year) + (1|block) + (1|genotype:year)'))
  fit <- lmer(f, data=df)
  v <- as.data.frame(VarCorr(fit))
  Vg <- v$vcov[v$grp == "genotype"]
  Vgxy <- v$vcov[v$grp == "genotype:year"]
  Ve <- v$vcov[v$grp == "Residual"]
  H2 <- round(Vg/(Vg + Vgxy/2 + Ve/4), 6)
  H2_lowN <- rbind(H2_lowN, data.frame("ASV"=trait, "H2_lowN"=H2))
}

H2_lowN

print("saving H2_lowN")
save(H2_lowN, file = "cache/asv_H2_lowN.rda")


















