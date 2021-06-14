
#group_level_H2_230genotypes.R

# calculate heritability from 2019 data (same data was used for GWAS)

# used 206 genotypes for which balanced data was available

# model: formula(paste0(trait, ' ~ (1|genotype) + (1|block)'))

# H2 formula: H2 <- round(Vg/(Vg + Ve/6), 6) 2 blocks x 3 subsamples = 6

library("phyloseq")
library("tidyverse")
library("ggrepel")



load("data/ps_grp.rda")



### filter asv table for 2019 data without checks
ps_19 <- subset_samples(ps_grp, year == "Y2019" & genotype != "CHECK")
asv_table <- rownames_to_column(data.frame(otu_table(ps_19)), var = "Sample_ID")
counts <- left_join(data.frame(sample_data(ps_19)), asv_table)


### calculate means of ASV counts of all subsamples by group
mean_counts <- counts %>%
  group_by(genotype, nitrogen, block, sp, spb) %>%
  summarize_each(funs(mean), -Sample_ID, -raw_seq_count, -filename, -year, -subsample, -collected_by_person, -pedigree, -subpopulation)


### calculate log relative abundance
asvtab <- mean_counts[, 7:ncol(mean_counts)]
logrel <- t(apply(asvtab, 1, function(x) log(x/sum(x) + 0.001)))
mean_counts_logrel <- cbind(data.frame(mean_counts[, 1:6]), logrel)



### split stdN, lowN samples, retain only genotxpes with complete set of 4 reps

h2dat_stdN <- mean_counts_logrel %>%
  filter(nitrogen == "+N") %>%
  group_by(genotype) %>%
  add_tally(name="count") %>%
  filter(count >= 2) %>%
  dplyr::select(-count)

# unique(as.character(h2dat_stdN$genotype))
# 206 genotypes

h2dat_lowN <- mean_counts_logrel %>%
  filter(nitrogen == "-N") %>%
  group_by(genotype) %>%
  add_tally(name="count") %>%
  filter(count >= 2) %>%
  dplyr::select(-count)

# unique(as.character(h2dat_lowN$genotype))
# 206 genotypes





#### calculate H2 ####




## calculate H2 +N
df <- h2dat_stdN
traits <- colnames(df)[7:ncol(df)]

H2_stdN <- data.frame()

for(trait in traits){
  #trait <- "asv_000013"
  f <- formula(paste0(trait, ' ~ (1|genotype) + (1|block)'))
  fit <- lmer(f, data=df)
  v <- as.data.frame(VarCorr(fit))
  Vg <- v$vcov[v$grp == "genotype"]
  Ve <- v$vcov[v$grp == "Residual"]
  H2 <- round(Vg/(Vg + Ve/6), 6)
  H2_stdN <- rbind(H2_stdN, data.frame("ASV"=trait, "H2_stdN_19"=H2))
}

H2_stdN




## calculate H2 -N
df <- h2dat_lowN
traits <- colnames(df)[7:ncol(df)]

H2_lowN <- data.frame()

for(trait in traits){
  #trait <- "asv_000013"
  f <- formula(paste0(trait, ' ~ (1|genotype) + (1|block)'))
  fit <- lmer(f, data=df)
  v <- as.data.frame(VarCorr(fit))
  Vg <- v$vcov[v$grp == "genotype"]
  Ve <- v$vcov[v$grp == "Residual"]
  H2 <- round(Vg/(Vg + Ve/6), 6)
  H2_lowN <- rbind(H2_lowN, data.frame("ASV"=trait, "H2_lowN_19"=H2))
}

H2_lowN




### H2 for BOTH N treatments


h2dat <- mean_counts_logrel %>%
  group_by(genotype) %>%
  add_tally(name="count") %>%
  filter(count >= 4) %>%
  dplyr::select(-count)




## calculate H2 -N
df <- h2dat
traits <- colnames(df)[7:ncol(df)]

H2_allN <- data.frame()

for(trait in traits){
  #trait <- "asv_000013"
  f <- formula(paste0(trait, ' ~ (1|genotype) + (1|nitrogen) + (1|block)'))
  fit <- lmer(f, data=df)
  v <- as.data.frame(VarCorr(fit))
  Vg <- v$vcov[v$grp == "genotype"]
  Ve <- v$vcov[v$grp == "Residual"]
  H2 <- round(Vg/(Vg + Ve/6), 6)
  H2_allN <- rbind(H2_allN, data.frame("ASV"=trait, "H2_allN_19"=H2))
}

H2_allN



### save data

load("data/group_data.rda")

#H2_19 <- H2_stdN %>%
#  left_join(H2_lowN) %>%
#  left_join(H2_allN)



#group_data <- group_data %>%
#  left_join(H2_19)

#save(group_data, file = "data/group_data.rda")

load("data/group_data.rda")

#### plot H2 +N vs abundance +N ####


stdN_h2ab_plot <- ggplot(group_data, aes(x=mean_blup_stdN, y=H2_stdN_19)) +
  geom_point() +
  geom_abline(coef = c(0,1), color = "red") +
  #geom_label_repel(aes(label = group_data$tax_group),
  #                 box.padding   = 0.05, 
  #                 point.padding = 0.5,
  #                 size = 2,
  #                 min.segment.length = 0.5,
  #                 segment.color = 'grey50') +
  theme_classic() +
  ylab("H2 +N") +
  xlab("Log(relative abundance) +N")


stdN_h2ab_plot




#### plot H2 -N vs abundance -N ####


lowN_h2ab_plot <- ggplot(group_data, aes(x=mean_blup_lowN, y=H2_lowN_19)) +
  geom_point() +
  geom_abline(coef = c(0,1), color = "red") +
  #geom_label_repel(aes(label = group_data$tax_group),
  #                 box.padding   = 0.05, 
  #                 point.padding = 0.5,
  #                 size = 2,
  #                 min.segment.length = 0.5,
  #                 segment.color = 'grey50') +
  theme_classic() +
  ylab("H2 -N") +
  xlab("Log(relative abundance) -N")


lowN_h2ab_plot





#### plot H2 +N vs -N ####


h2_plot <- ggplot(group_data, aes(x=H2_lowN_19, y=H2_stdN_19)) +
  geom_point() +
  geom_abline(coef = c(0,1), color = "red") +
  #geom_label_repel(aes(label = group_data$tax_group),
  #                box.padding   = 0.05, 
  #                point.padding = 0.5,
  #                size = 2,
  #                min.segment.length = 0.5,
  #                segment.color = 'grey50') +
  theme_classic() +
  ylab("H2 +N") +
  xlab("H2 -N")


h2_plot



# permutation test to find threshold

mean(group_data$H2_lowN_19)
mean(group_data$H2_stdN_19)


# lets calculate the absolute diff in means
test.stat1 <- abs(mean(group_data$H2_lowN_19) - 
                    mean(group_data$H2_stdN_19)) 

perm_data <- group_data %>%
  dplyr::select(tax_group, H2_lowN_19, H2_stdN_19) %>%
  pivot_longer(cols = starts_with("H2_"), names_to = "nitrogen", values_to = "H2") %>%
  mutate(nitrogen = ifelse(nitrogen == 'H2_stdN_19', '+N', '-N'))

?pivot_longer

set.seed(2021)

n <- length(perm_data$nitrogen)
P <- 100000
variable <- perm_data$H2

PermSamples <- matrix(0, nrow = n, ncol = P)

for(i in 1:P){
  PermSamples[,i] <- sample(variable, size = n, replace = FALSE)
}

PermSamples[,1:5]


# initialize vectors to store all of the Test-stats:
means_lowN <- means_stdN <- permutationMeans <- rep(0, P)

# loop thru, and calculate the test-stats
for (i in 1:P){
  # calculate the perm-test-stat1 and save it
  means_lowN[i] <- mean(PermSamples[perm_data$nitrogen=="-N",i])
  means_stdN[i] <- mean(PermSamples[perm_data$nitrogen=="+N",i])
  permutationMeans[i] <- mean(PermSamples[,i])
  
}


mean(perm_data$H2) # overall mean 0.335919

threshold_lowN <- mean(means_lowN) # permutation mean 0.3359948 use this as threshold heritable/not heritable
mean(group_data$H2_lowN_19) # data mean 0.352033

threshold_stdN <- mean(means_stdN) # permutation mean 0.3358431 use this as threshold heritable/not heritable
mean(group_data$H2_stdN_19) # data mean 0.3198049





h2_plotperm <- ggplot(group_data, aes(x=H2_lowN_19, y=H2_stdN_19)) +
  geom_density_2d() +
  geom_point(alpha=0.7) +
  geom_abline(coef = c(0,1), linetype="dashed") +
  geom_vline(xintercept = threshold_lowN, color = "red") +
  geom_hline(yintercept = threshold_stdN, color = "red") +
  geom_smooth(method='lm', se=TRUE, color="#008000") +
  #geom_label_repel(aes(label = group_data$tax_group),
  #                box.padding   = 0.05, 
  #                point.padding = 0.5,
  #                size = 2,
  #                min.segment.length = 0.5,
  #                segment.color = 'grey50') +
  theme_classic() +
  ylab("H2 +N") +
  xlab("H2 -N")


h2_plotperm




?geom_smooth







# calculate p value for all 100k permutations

mean(Perm.test.stat1 >= test.stat1)


#### plot H2 +N vs -N, overlay labels for groups with strong GWAS signal ####

microbe_traits <- read_csv("data/microbial_traits.csv")


top_traits <- microbe_traits %>%
  filter(GWAS_signal == "strong") %>%
  select(tax_group, N_treatment) %>%
  rename(GWAS_signal_N = N_treatment) %>%
  left_join(group_data)



group_data <- left_join(group_data, top_traits[, c("tax_group", "GWAS_signal_N")])
group_data$GWAS_signal_N <- ifelse(is.na(group_data$GWAS_signal_N), "none", group_data$GWAS_signal_N)
group_data$GWAS_signal_N <- factor(group_data$GWAS_signal_N, levels=c("stdN", "lowN", "none" ))


h2_plot <- ggplot(group_data, aes(x=H2_lowN_19, y=H2_stdN_19)) +
#h2_plot <- ggplot(group_data, aes(x=H2_lowN_19, y=H2_stdN_19, color=GWAS_signal_N)) +
#h2_plot <- ggplot(group_data, aes(x=H2_lowN, y=H2_stdN, color=GWAS_signal_N)) +
  geom_point(size = 2) +
  geom_abline(coef = c(0,1), color = "red") +
  #geom_text_repel(data = . %>%
  #                   mutate(label = ifelse(tax_group %in% top_traits$tax_group,
  #                                         tax_group, "")),
  #                 aes(label = label),
  #                 #box.padding = 2,
  #                 #size = 2,
  #                 min.segment.length = 0,
  #                 max.overlaps = 150,
  #                 segment.color = 'grey50') +
  #scale_color_manual(values = c("lowN"="#C00001", "stdN"="#008000", "none"="#cccccc")) +
  theme_classic() +
  theme(text = element_text(size=20)) +
  labs(y = bquote(~H^2~'under +N treatment') , x = bquote(~H^2~'under -N treatment'))
  


h2_plot



#### plot abundance (mean BLUP), overlay labels for groups with strong GWAS signal ####



blup_plot <- ggplot(group_data, aes(x=mean_blup_lowN, y=mean_blup_stdN, color=GWAS_signal_N)) +
  geom_point(size = 2) +
  geom_abline(coef = c(0,1), color = "red") +
  geom_text_repel(data = . %>%
                    mutate(label = ifelse(tax_group %in% top_traits$tax_group,
                                          tax_group, "")),
                  aes(label = label),
                  #box.padding = 2,
                  #size = 2,
                  min.segment.length = 0,
                  max.overlaps = 150,
                  segment.color = 'grey50') +
  scale_color_manual(values = c("lowN"="#C00001", "stdN"="#008000", "none"="#cccccc")) +
  theme_classic() +
  ylab("Log(relative abundance under +N treatment)") +
  xlab("Log(relative abundance under -N treatment)")


blup_plot



#### plot heritability vs abundance



h2ab_stdN <- ggplot(group_data, aes(x=mean_blup_stdN, y=H2_stdN_19, color=GWAS_signal_N)) +
  geom_point(size = 2) +
  geom_text_repel(data = . %>%
                    mutate(label = ifelse(tax_group %in% top_traits$tax_group,
                                          tax_group, "")),
                  aes(label = label),
                  #box.padding = 2,
                  #size = 2,
                  min.segment.length = 0,
                  max.overlaps = 150,
                  segment.color = 'grey50') +
  scale_color_manual(values = c("lowN"="#C00001", "stdN"="#008000", "none"="#cccccc")) +
  theme_classic() +
  labs(y = bquote(~H^2~'under +N treatment') , x = "Log(relative abundance under +N treatment)")


h2ab_stdN




h2ab_lowN <- ggplot(group_data, aes(x=mean_blup_lowN, y=H2_lowN_19, color=GWAS_signal_N)) +
  geom_point(size = 2) +
  geom_text_repel(data = . %>%
                    mutate(label = ifelse(tax_group %in% top_traits$tax_group,
                                          tax_group, "")),
                  aes(label = label),
                  #box.padding = 2,
                  #size = 2,
                  min.segment.length = 0,
                  max.overlaps = 150,
                  segment.color = 'grey50') +
  scale_color_manual(values = c("lowN"="#C00001", "stdN"="#008000", "none"="#cccccc")) +
  theme_classic() +
  labs(y = bquote(~H^2~'under -N treatment') , x = "Log(relative abundance under -N treatment)")


h2ab_lowN



#### plot heritability vs S score


### merge S score data
load("data/data_summary_150_traits.rda")

data_summary_150_traits <- left_join(data_summary_150_traits, group_data[, c("tax_group","H2_stdN_19", "H2_lowN_19", "H2_allN_19", "GWAS_signal_N")])

#save(data_summary_150_traits, file = "data/data_summary_150_traits.rda")




h2S_stdN <- ggplot(data_summary_150_traits, aes(x=stdN_Mean_S, y=H2_stdN_19, color=GWAS_signal_N)) +
  geom_point(size = 2) +
  geom_text_repel(data = . %>%
                    mutate(label = ifelse(tax_group %in% top_traits$tax_group,
                                          tax_group, "")),
                  aes(label = label),
                  #box.padding = 2,
                  #size = 2,
                  min.segment.length = 0,
                  max.overlaps = 150,
                  segment.color = 'grey50') +
  scale_color_manual(values = c("lowN"="#C00001", "stdN"="#008000", "none"="#cccccc")) +
  theme_classic() +
  labs(y = bquote(~H^2~'under +N treatment') , x = "selection effect S")


h2S_stdN





h2S_lowN <- ggplot(data_summary_150_traits, aes(x=lowN_Mean_S, y=H2_lowN_19, color=GWAS_signal_N)) +
  geom_point(size = 2) +
  geom_text_repel(data = . %>%
                    mutate(label = ifelse(tax_group %in% top_traits$tax_group,
                                          tax_group, "")),
                  aes(label = label),
                  #box.padding = 2,
                  #size = 2,
                  min.segment.length = 0,
                  max.overlaps = 150,
                  segment.color = 'grey50') +
  scale_color_manual(values = c("lowN"="#C00001", "stdN"="#008000", "none"="#cccccc")) +
  theme_classic() +
  labs(y = bquote(~H^2~'under -N treatment') , x = "selection effect S")


h2S_lowN


