## calculate BLUP for 150 tax groups

library("phyloseq")
library("tidyverse")
library("lme4")

#devtools::install_github("jyanglab/g3tools")
library("g3tools")

#### prepare data to calculate BLUP ####

load("data/ps_grp.rda")


### use 2019 data only, exclude checks
ps19 <- subset_samples(ps_grp, year == "Y2019" & genotype != "CHECK")
asv_table <- rownames_to_column(data.frame(otu_table(ps19)), var = "Sample_ID")
counts <- left_join(data.frame(sample_data(ps19)), asv_table)

### calculate means of ASV counts of all subsamples by group
mean_counts <- counts %>%
  group_by(genotype, nitrogen, block, sp, spb) %>%
  summarize_each(funs(mean), -year, -Sample_ID, -raw_seq_count, -filename, -subsample, -collected_by_person, -pedigree, -subpopulation)


### calculate log relative abundance
asvtab <- mean_counts[, 7:ncol(mean_counts)]
logrel <- t(apply(asvtab, 1, function(x) log(x/sum(x) + 0.001)))
mean_counts_logrel <- cbind(data.frame(mean_counts[, 1:6]), logrel)



### BLUP std N

if(!(exists("blup_stdN"))){
  blup_stdN <- data.frame(genotype = sort(unique(df$genotype)))
}

for (asvid in colnames(asvtab)){
  
df <- mean_counts_logrel %>%
  filter(nitrogen == "+N") %>%
  select_("genotype", "block", "sp", "spb", "row", asvid)

colnames(df)[ncol(df)] <- "asv_abundance"

get_BLUP(data = df, model = asv_abundance ~ (1 | genotype) + (1 | block) + (1 | sp) + (1 | spb), which.factor = "genotype",
  outfile = "cache/BLUP/tmp_blup.csv")

blup <- read.csv("cache/BLUP/tmp_blup.csv")
colnames(blup)[1] <- "genotype"
colnames(blup)[2] <- asvid

blup_stdN <- left_join(blup_stdN, blup)

}



### BLUP low N

if(!(exists("blup_lowN"))){
  blup_lowN <- data.frame(genotype = sort(unique(df$genotype)))
}

for (asvid in colnames(asvtab)){
  
  df <- mean_counts_logrel %>%
    filter(nitrogen == "-N") %>%
    select_("genotype", "block", "sp", "spb", "row", asvid)
  
  colnames(df)[ncol(df)] <- "asv_abundance"
  
  get_BLUP(data = df, model = asv_abundance ~ (1 | genotype) + (1 | block) + (1 | sp) + (1 | spb), which.factor = "genotype",
    outfile = "cache/BLUP/tmp_blup.csv")
  
  blup <- read.csv("cache/BLUP/tmp_blup.csv")
  colnames(blup)[1] <- "genotype"
  colnames(blup)[2] <- asvid
  
  blup_lowN <- left_join(blup_lowN, blup)
  
}






### save data

write_csv(blup_stdN, file="data/blup_stdN_150_tax_groups.csv")
write_csv(blup_lowN, file="data/blup_lowN_150_tax_groups.csv")





