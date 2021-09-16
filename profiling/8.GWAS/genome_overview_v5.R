#genome_overview_v5

## redone gwas with SNP map version 5


library("tidyverse")
library("RcppRoll")
library("LDheatmap")
library("genetics")


# load v5 data
gwas_HN <- read_csv("largedata/GWAS2/merged_HNshort_5.csv")
gwas_LN <- read_csv("largedata/GWAS2/merged_LNshort_5.csv")


gwas_HN$log10p <- -log10(gwas_HN$p_wald)



## v4 data for comparison
load("largedata/combined_SNP_data.rda")
nrow(filter(gwas_dat, nitrogen == "+N"))
# v4: 61339 snps over threshold 5
# v5: 46908 snps over threshold 5
nrow(filter(gwas_dat, nitrogen == "-N"))
# v4: 54774 snps over threshold 5
# v5: 135743 snps over threshold 5

gwas_HN$nitrogen <- "HN"
gwas_LN$nitrogen <- "LN"

gwas_dat_v5 <- rbind(gwas_HN, gwas_LN)
gwas_dat <- gwas_dat_v5


gwas_dat$log10p <- log10(gwas_dat$p_wald)

load("largedata/gwas_dat_v5.rda")

## remove outlier Chryseobacterium sp 2, asv_000483


gwas_dat <- filter(gwas_dat, trait != "asv_000483")

save(gwas_dat, file="largedata/gwas_dat_v5.rda")

load("data/data_summary_150_traits.rda")


chr_len_v5 <- c(308452471, 243675191, 238017767, 250330460, 226353449, 181357234, 185808916, 182411202, 163004744, 152435371)

chr_len_10k <- c(308460000, 243680000, 238020000, 250340000, 226360000, 181360000, 185810000, 182420000, 163010000, 152440000)




nCHR <- length(unique(gwas_dat$chr))
chindex <- c(0) ## this number needs to be added to position within each chromosome
for (i in c(1:nCHR)){
  #print(i)
  chindex[i+1] <- sum(chr_len_10k[1:i])
}
chindex <- chindex[1:10]
bpadd <- data.frame("chr" = c(1:nCHR), "chindex" = chindex)



options(scipen = 999)

# add bins
plot_dat <- gwas_dat %>%
  dplyr::select(-X1) %>%
  mutate(log10p = -log10(p_wald)) %>%
  #filter(log10p >= 5) %>%
  left_join(bpadd, by = "chr") %>%
  mutate(psabs = chindex + ps) %>%
  mutate(psabs = as.numeric(gsub('.{4}$', '', psabs))) %>% ## pos. in 10kb units
  mutate(ps = as.character(ps)) %>%
  mutate(bin = as.numeric(gsub('.{4}$', '', ps)) + 1) %>% # absolute bin position in 10kb units
  group_by(chr, bin, nitrogen) %>%
  add_count(count = sum(-log10(p_wald) >= 5), name="total_sign_snps") %>%
  ##mutate(total_mean_p = mean(log10p)) %>%
  group_by(chr, bin, nitrogen, trait) %>%
  add_count(count = sum(-log10(p_wald) >= 5), name="sign_snps") %>%
  ##mutate(mean_p = mean(log10p)) %>%
  rename(ASV=trait) %>%
  left_join(dplyr::select(data_summary_150_traits, ASV, tax_group)) %>%
  #mutate(color = ifelse(chr %in% c(1,3,5,7,9), "#276FBF", "#183059")) %>%
  mutate(chr = as.character(chr)) 





# get taxa per bin
taxa_per_bin <- plot_dat %>%
  dplyr::select(nitrogen, chr, bin, tax_group) %>%
  unique() %>%
  group_by(nitrogen, chr, bin) %>%
  add_tally(name="taxa_per_bin") 

plot_dat <- plot_dat %>%
  left_join(taxa_per_bin)

plot_dat$nitrogen <- factor(plot_dat$nitrogen, levels = c("LN", "HN"))





###### find outlier ######
total_sign_snps <- plot_dat %>%
  group_by(nitrogen, tax_group) %>%
  dplyr::select(nitrogen, tax_group) %>%
  tally(name="total_sign_snps")

ggplot(total_sign_snps, aes(y=total_sign_snps, x=tax_group)) + 
  geom_bar(stat = "identity") +
  #facet_wrap(~nitrogen, nrow=2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## v4 data for comparison
nrow(filter(plot_dat, nitrogen == "HN"))
# v4: 61339 snps over threshold 5
# v5: 46908 snps over threshold 5
# v5: 46830 snps over threshold 5 (outlier removed)
nrow(filter(plot_dat, nitrogen == "LN"))
# v4: 54774 snps over threshold 5
# v5: 135743 snps over threshold 5
# v5: 43252 snps over threshold 5 (outlier removed)

###### find outlier ######

#save(plot_dat, file="largedata/GWAS2/plot_dat_v5_gt5.rda")



colors = c("1" = "#276FBF", "3" = "#276FBF", "5" = "#276FBF", "7" = "#276FBF", "9" = "#276FBF",
           "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059", "sig" = "#ff0000")



total_snp_dat <- plot_dat %>%
  filter(sign_snps >= 1) %>%
  ungroup() %>%
  dplyr::select(nitrogen, chr, bin, psabs, total_sign_snps, taxa_per_bin) %>%
  unique()


short_snp_dat <- plot_dat %>%
  filter(total_sign_snps >= 23) %>%
  ungroup() %>%
  dplyr::select(nitrogen, chr, bin, psabs, total_sign_snps, taxa_per_bin) %>%
  unique()


## find position for labels
axis_set <- plot_dat %>% 
  group_by(chr) %>% 
  summarize(center = (max(psabs) + min(psabs)) / 2) # find middle of chromosome (to place labels)




overview_plot <-  ggplot(total_snp_dat, aes(x = psabs, y = total_sign_snps, color = chr)) +
  #geom_area() +
  geom_bar(stat="identity", position="dodge") +
  geom_point(data=short_snp_dat, aes(size = taxa_per_bin), alpha = 0.3) +
  #geom_hline(yintercept = 6, color = "red", linetype = "dashed") + # 50% quantile
  geom_hline(yintercept = 26, color = "red", linetype = "dashed") + # 75% quantile
  geom_hline(yintercept = 91, color = "red", linetype = "dashed") + # 90% quantile
  #scale_x_continuous(label = axis_set$chr, breaks = axis_set$center, limits = c(1273000000, 1273300000)) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  #geom_vline(xintercept = vlines) +
  scale_color_manual(values = colors) +
  scale_size_continuous(breaks = c(1, 5, 10)) +
  #scale_size_continuous(range = c(0.5,3)) +
  #scale_y_continuous(trans='log10') +
  #scale_y_log10() +
  ##coord_trans(y='log10') +
  facet_wrap(~nitrogen, nrow = 2) +
  labs(x = NULL, y = "Total sign. SNPs per 10kb bin") + 
  #ggtitle("total signif. SNPs for 399 10kb bins with nearby annotated genes") +
  geom_vline(xintercept =  c(0, 212850), color = "red", size=0.1) + # to combine plots
  geom_hline(yintercept = 0, size=0.1) + 
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

overview_plot




test <- plot_dat %>%
  filter(total_sign_snps >= 26)


top75pc_taxa <- unique(test$tax_group)

save(top75pc_taxa, file="cache/RNAseq/top75pc_taxa.rda")

test <- plot_dat %>%
  filter(total_sign_snps >= 91)


top90pc_taxa <- unique(test$tax_group)

save(top90pc_taxa, file="cache/RNAseq/top90pc_taxa.rda")



### find threshold: plot sign snps vs taxa retained

### plot number of traits retained for 10-100 sign. SNPs


threshold_plot_data <- data.frame(matrix(ncol = 4, nrow = 0))

for(threshold in seq(1,100)) {
  top <- plot_dat %>%
    filter(total_sign_snps >= threshold) %>%
    dplyr::select(nitrogen, chr, bin, total_sign_snps, taxa_per_bin, tax_group) %>%
    distinct()
  
  threshold_plot_data <- rbind(threshold_plot_data, c(threshold, length(unique(top$tax_group)), length(unique(top$gene)), length(unique(top$bin))))
  
}

colnames(threshold_plot_data) <- c("total_sign_snps", "taxa_retained", "genes_retained", "bins_retained")


quantile(filter(plot_dat, nitrogen == "LN")$total_sign_snps, c(.25, .5, .75, .9))
mean(filter(plot_dat, nitrogen == "LN")$total_sign_snps)


quantile(filter(plot_dat, nitrogen == "HN")$total_sign_snps, c(.25, .5, .75, .9))
mean(filter(plot_dat, nitrogen == "HN")$total_sign_snps)





ggplot(threshold_plot_data, aes(x=total_sign_snps, y=taxa_retained)) +
  geom_point() +
  geom_vline(xintercept = 6, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 26, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 91, color = "red", linetype = "dashed") +
  theme_bw()


ggplot(threshold_plot_data, aes(x=total_sign_snps, y=bins_retained)) +
  geom_point() +
  geom_vline(xintercept = 6, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 26, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 91, color = "red", linetype = "dashed") +
  theme_bw()





v5_MAPLs <- plot_dat %>%
  filter(total_sign_snps >= 26)



save(v5_MAPLs, file="cache/GWAS2/v5_MAPLs.rda")





