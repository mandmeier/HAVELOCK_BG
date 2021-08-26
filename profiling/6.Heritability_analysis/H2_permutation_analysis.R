# permutation_test
# This permutation test is to determine a threshold to distinguish
# heritable vs non heritable traits
# the threshold is the heritability score you would expect by chance


library("tidyverse")
library("lme4")


# 1) for each trait and N treatment: calculate "observed" H2
# 2) shuffle each genotype's ASV counts of this trait 1000 times
# 3) calculate 1000 permutation H2 scores
# 4) how many of the 1000 are larger than the observed H2?
# 5) p-value = #larger/1000
# 6) this yields 300 p values for each trait and treatment
# 7) mark quadrants: check if p < 0.05 for +N and -N

# get H2 data: log(relative abundance of ASV counts for each trait), traits identified by representative ASV
# data is from group_level_H2_230genotypes.R

load("cache/h2dat_stdN.rda")
load("cache/h2dat_lowN.rda")


# Function to calculate h2

getH2 <- function(h2dat) {
  #h2dat <- h2dat_stdN
  traits <- colnames(h2dat)[7:ncol(h2dat)]
  #print(traits)
  H2_df <- data.frame()
  
  for(trait in traits){
    #trait <- "asv_000013"
    f <- formula(paste0(trait, ' ~ (1|genotype) + (1|block)'))
    fit <- suppressMessages(lmer(f, data=h2dat))
    v <- as.data.frame(VarCorr(fit))
    Vg <- v$vcov[v$grp == "genotype"]
    Ve <- v$vcov[v$grp == "Residual"]
    H2 <- round(Vg/(Vg + Ve/6), 6)
    H2_df <- rbind(H2_df, data.frame("ASV"=trait, "H2"=H2))
  }
  
  return(H2_df)
  
}

# function to shuffle abundances
shuffle <- function(h2dat) {
  h2dat$genotype <- sample(h2dat$genotype)
  return(h2dat)
}
  


#### calculate permutation p_values stdN

set.seed(2021)
H2_permutations_stdN <-  getH2(h2dat_stdN)
colnames(H2_permutations_stdN)[2] <- "obs"

for( i in c(1:1000)){
  print(paste("permutation", i))
  perm <-  getH2(shuffle(h2dat_stdN))
  colnames(perm)[2] <- paste0("p",i)
  H2_permutations_stdN <- suppressMessages(left_join(H2_permutations_stdN, perm))
}

p_values_stdN <- H2_permutations_stdN %>%
  rowwise() %>%
  mutate(perm_p_stdN = (sum(c_across(3:length(H2_permutations_stdN))> obs)+1)/(length(H2_permutations_stdN)-1)) %>%
  dplyr::select(ASV, perm_p_stdN)



#### calculate permutation p_values lowN

set.seed(2021)
H2_permutations_lowN <-  getH2(h2dat_lowN)
colnames(H2_permutations_lowN)[2] <- "obs"

for( i in c(1:1000)){
  print(paste("permutation", i))
  perm <-  getH2(shuffle(h2dat_lowN))
  colnames(perm)[2] <- paste0("p",i)
  H2_permutations_lowN <- suppressMessages(left_join(H2_permutations_lowN, perm))
}

p_values_lowN <- H2_permutations_lowN %>%
  rowwise() %>%
  mutate(perm_p_lowN = (sum(c_across(3:length(H2_permutations_lowN))> obs)+1)/(length(H2_permutations_lowN)-1)) %>%
  dplyr::select(ASV, perm_p_lowN)


### save permutation data

save(H2_permutations_stdN, file="cache/H2_permutations_stdN.rda")
save(H2_permutations_lowN, file="cache/H2_permutations_lowN.rda")




### add pvalues to group data


load("data/group_data.rda")

group_data <- group_data %>%
  left_join(p_values_stdN) %>%
  left_join(p_values_lowN)



#save(group_data, file = "data/group_data.rda")




