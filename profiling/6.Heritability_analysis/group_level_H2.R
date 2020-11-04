library("phyloseq")
library("tidyverse")
library("lme4")

#### prepare data to calculate H2 ####

load("data/ps_grp.rda")

# find genotypes common to both years, without check
sdat <- data.frame(sample_data(ps_grp))
y2018 <- sdat %>%
  filter(year == "Y2018")
y2019 <- sdat %>%
  filter(year == "Y2019")
common_genotypes <- Reduce(intersect, list(unique(y2018$genotype),unique(y2019$genotype)))
common_genotypes <- common_genotypes[-1]


### filter ps object with common genotypes
ps_common <- subset_samples(ps_grp, genotype %in% common_genotypes)
asv_table <- rownames_to_column(data.frame(otu_table(ps_common)), var = "Sample_ID")
counts <- left_join(data.frame(sample_data(ps_common)), asv_table)


### calculate means of ASV counts of all subsamples by group
mean_counts <- counts %>%
  group_by(year, genotype, nitrogen, block, sp, spb) %>%
  summarize_each(funs(mean), -Sample_ID, -raw_seq_count, -filename, -subsample, -collected_by_person, -pedigree, -subpopulation)
  

### calculate log relative abundance
asvtab <- mean_counts[, 8:ncol(mean_counts)]
logrel <- t(apply(asvtab, 1, function(x) log(x/sum(x) + 0.001)))
mean_counts_logrel <- cbind(data.frame(mean_counts[, 1:7]), logrel)


### split stdN, lowN samples, retain only genotxpes with complete set of 4 reps

h2dat_stdN <- mean_counts_logrel %>%
  filter(nitrogen == "+N") %>%
  group_by(genotype) %>%
  add_tally(name="count") %>%
  filter(count == 4) %>%
  dplyr::select(-count)
# 20 genotypes

h2dat_lowN <- mean_counts_logrel %>%
  filter(nitrogen == "-N") %>%
  group_by(genotype) %>%
  add_tally(name="count") %>%
  filter(count == 4) %>%
  dplyr::select(-count)
# 19 genotypes

## add ASV counts to group data
#load("data/group_data.rda")
#calculate ASV counts
#stdN_counts <- colSums(otu_table(subset_samples(ps_grp, nitrogen == "+N")))
#stdN_counts <- rownames_to_column(data.frame("count_stdN"=stdN_counts), var = "ASV")
#lowN_counts <- colSums(otu_table(subset_samples(ps_grp, nitrogen == "-N")))
#lowN_counts <- rownames_to_column(data.frame("count_lowN"=lowN_counts), var = "ASV")
#group_data <- left_join(group_data, stdN_counts)
#group_data <- left_join(group_data, lowN_counts)




#### calculate H2 ####


## calculate H2 +N
df <- h2dat_stdN
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


## calculate H2 -N
df <- h2dat_lowN
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

### add H2 data to group data
#group_data <- left_join(group_data, H2_stdN)
#group_data <- left_join(group_data, H2_lowN)
#save(group_data, file = "data/group_data.rda")



















