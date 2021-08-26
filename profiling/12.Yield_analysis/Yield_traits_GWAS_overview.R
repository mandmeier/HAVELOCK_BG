## yield analysis using 13 yield traits
library("tidyverse")



### combine GWAS data for all traits and treatments

files <- list.files(path="largedata/Yield_Traits_GWAS", pattern="*.txt", full.names=TRUE, recursive=FALSE)
dummy <- read_delim(file = "largedata/Yield_Traits_GWAS/canopy-HN-12-Aug.assoc.txt", delim = "\t")

Yield_GWAS <- filter(dummy, ps == "x")
Yield_GWAS$trait <- ""
Yield_GWAS$nitrogen <- ""

for(file in files){
  #file <- "largedata/Yield_Traits_GWAS/canopy-HN-12-Aug.assoc.txt"
  if (grepl("LN", file, fixed=TRUE)){
    N <- "-N"
  } else {
    N <- "+N"
  }
  ext <- str_split(file, "-")[[1]][2]
  trait <- str_split(ext, "\\.")[[1]][1]
  
  gwasdat <- read_delim(file = file, delim = "\t")
  gwasdat$trait <- trait
  gwasdat$nitrogen <- N
  Yield_GWAS <- rbind(Yield_GWAS, gwasdat)
  
}



#save(Yield_GWAS, file="cache/Yield_traits_GWAS/Yield_GWAS.rda")


unique(Yield_GWAS$trait)


##load("cache/Yield_traits_GWAS/Yield_GWAS.rda")


length(unique(Yield_GWAS$rs))


### add chromosome index and absolute position


gwas_dat <- Yield_GWAS

chr_lengths <- unique(dplyr::select(gwas_dat, chr, chindex))
#chr_len <- c(306970668, 244420837, 235651057, 246967222, 223706090, 173536910, 181718119, 181046068, 159686117, 150929986)
chr_len_10k <- c(306980000, 244430000, 235660000, 246970000, 223710000, 173540000, 181720000, 181050000, 159690000, 150930000)
#vlines <- c(0,cumsum(chr_len)) ## draw vertical lines at chomosome ends
nCHR <- length(unique(gwas_dat$chr))
chindex <- c(0) ## this number needs to be added to position within each chromosome
for (i in c(1:nCHR)){
  #print(i)
  chindex[i+1] <- sum(chr_len_10k[1:i])
}
chindex <- chindex[1:10]
bpadd <- data.frame("chr" = c(1:nCHR), "chindex" = chindex)

# add bins
plot_dat <- gwas_dat %>%
  mutate(log10p = -log10(p_wald)) %>%
  filter(log10p >= 3) %>%
  left_join(bpadd, by = "chr") %>%
  mutate(psabs = chindex + ps) %>%
  mutate(psabs = as.numeric(gsub('.{4}$', '', psabs))) %>% ## pos. in 10kb units
  mutate(bin = as.numeric(gsub(paste0(".{",4,"}$"), '', ps)) + 1) %>% # absolute bin position
  group_by(chr, bin, nitrogen) %>%
  add_count(count = sum(log10p >= 3), name="total_sign_snps") %>%
  mutate(total_mean_p = mean(log10p)) %>%
  group_by(chr, bin, nitrogen, trait) %>%
  add_count(count = sum(log10p >= 3), name="sign_snps") %>%
  #mutate(color = ifelse(chr %in% c(1,3,5,7,9), "#276FBF", "#183059")) %>%
  mutate(chr = as.character(chr)) #%>%
#mutate(ropp_sign_snps = frollmean(sign_snps, n = 10, align = "center")) 

# get traits per bin
traits_per_bin <- plot_dat %>%
  dplyr::select(nitrogen, chr, bin, trait) %>%
  unique() %>%
  group_by(nitrogen, chr, bin) %>%
  add_tally(name="traits_per_bin") 

plot_dat <- plot_dat %>%
  left_join(traits_per_bin)




colors = c("1" = "#276FBF", "3" = "#276FBF", "5" = "#276FBF", "7" = "#276FBF", "9" = "#276FBF",
           "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059", "sig" = "#ff0000")


colors = c("1" = "#fcae1e", "3" = "#fcae1e", "5" = "#fcae1e", "7" = "#fcae1e", "9" = "#fcae1e",
           "2" = "#dd571c", "4" = "#dd571c", "6" = "#dd571c", "8" = "#dd571c", "10" = "#dd571c", "sig" = "#ff0000")



## find position for labels
axis_set <- plot_dat %>% 
  group_by(chr) %>% 
  summarize(center = (max(psabs) + min(psabs)) / 2) # find middle of chromosome (to place labels)





chr_boundaries <- c(0)

chr_len_10k <- c(306980000, 244430000, 235660000, 246970000, 223710000, 173540000, 181720000, 181050000, 159690000, 150930000)

for (i in c(1:10)) {
  sum <- sum(chr_len_10k[1:i])
  chr_boundaries <- c(chr_boundaries, sum/10000)
}


yield_gwas_overview <- plot_dat


save(yield_gwas_overview, file="cache/Yield_traits_GWAS/yield_gwas_overview.rda")


overview_plot <-  ggplot(plot_dat, aes(x = psabs, y = total_sign_snps, color = chr)) +
  #geom_area() +
  geom_bar(stat="identity", position="dodge") +
  geom_point(data=plot_dat, aes(size = traits_per_bin), alpha = 0.3) +
  #geom_hline(yintercept = 23, color = "red", linetype = "dashed") + 
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
  #geom_vline(xintercept =  c(0, 210460), color = "red", size=1) + # to combine plots
  geom_vline(xintercept =  chr_boundaries, color = "black", size=0.1) + # to combine plots
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



### compare to microbe traits gwas overview

#load("cache/Yield_traits_GWAS/microbe_gwas_overview.rda")


# +N chr 8
test <- microbe_gwas_overview %>%
  filter(psabs == 165149)

# -N chr 8
test <- microbe_gwas_overview %>%
  filter(psabs == 164803)


# -N chr 2
test <- microbe_gwas_overview %>%
  filter(psabs == 33736)




