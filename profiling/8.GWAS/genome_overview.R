library("tidyverse")
library("RcppRoll")
library("LDheatmap")
library("genetics")



load("largedata/combined_SNP_data.rda")


gwas_dat$log10p <- log10(gwas_dat$p_wald)

load("data/data_summary_150_traits.rda")

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
  dplyr::select(-chindex, -psabs) %>%
  mutate(log10p = -log10(p_wald)) %>%
  filter(log10p >= 5) %>%
  left_join(bpadd, by = "chr") %>%
  mutate(psabs = chindex + ps) %>%
  mutate(psabs = as.numeric(gsub('.{4}$', '', psabs))) %>% ## pos. in 10kb units
  mutate(bin = as.numeric(gsub(paste0(".{",4,"}$"), '', ps)) + 1) %>% # absolute bin position
  group_by(chr, bin, nitrogen) %>%
  add_count(count = sum(-log10(p_wald) >= 5), name="total_sign_snps") %>%
  mutate(total_mean_p = mean(log10p)) %>%
  group_by(chr, bin, nitrogen, trait) %>%
  add_count(count = sum(-log10(p_wald) >= 5), name="sign_snps") %>%
  mutate(mean_p = mean(log10p)) %>%
  left_join(dplyr::select(data_summary_150_traits, trait, tax_group)) %>%
  #mutate(color = ifelse(chr %in% c(1,3,5,7,9), "#276FBF", "#183059")) %>%
  mutate(chr = as.character(chr)) #%>%
  #mutate(ropp_sign_snps = frollmean(sign_snps, n = 10, align = "center")) 

# get taxa per bin
taxa_per_bin <- plot_dat %>%
  dplyr::select(nitrogen, chr, bin, tax_group) %>%
  unique() %>%
  group_by(nitrogen, chr, bin) %>%
  add_tally(name="taxa_per_bin") 

plot_dat <- plot_dat %>%
  left_join(taxa_per_bin)



microbe_gwas_overview  <- plot_dat

#save(microbe_gwas_overview, file="cache/Agro_traits_GWAS/microbe_gwas_overview.rda")





colors = c("1" = "#276FBF", "3" = "#276FBF", "5" = "#276FBF", "7" = "#276FBF", "9" = "#276FBF",
           "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059", "sig" = "#ff0000")


total_snp_dat <- plot_dat %>%
  filter(sign_snps >= 1) %>%
  ungroup() %>%
  dplyr::select(nitrogen, chr, bin, psabs, total_sign_snps, taxa_per_bin, total_mean_p) %>%
  unique() %>%
  mutate(sxp = total_sign_snps*total_mean_p) %>%
  mutate(sxt = total_sign_snps*taxa_per_bin)
  

short_snp_dat <- plot_dat %>%
  filter(total_sign_snps >= 23) %>%
  ungroup() %>%
  dplyr::select(nitrogen, chr, bin, psabs, total_sign_snps, taxa_per_bin, total_mean_p) %>%
  unique() %>%
  mutate(sxp = total_sign_snps*total_mean_p) %>%
  mutate(sxt = total_sign_snps*taxa_per_bin)


## find position for labels
axis_set <- plot_dat %>% 
  group_by(chr) %>% 
  summarize(center = (max(psabs) + min(psabs)) / 2) # find middle of chromosome (to place labels)




overview_plot <-  ggplot(total_snp_dat, aes(x = psabs, y = total_sign_snps, color = chr)) +
  #geom_area() +
  geom_bar(stat="identity", position="dodge") +
  geom_point(data=short_snp_dat, aes(size = taxa_per_bin), alpha = 0.3) +
  geom_hline(yintercept = 23, color = "red", linetype = "dashed") + 
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
  geom_vline(xintercept =  c(0, 165149, 210460), color = "red", size=0.1) + # to combine plots
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


### how many loci above threshold? ## 360 loci over threshold 23 (red line)

summary <- total_snp_dat %>%
  filter(total_sign_snps >= 23)

## How many traits? 921 locus x trait combinations, 127 traits

### get associated traits and pick largest snp in each locus

top_loci_traits_all_snps <- plot_dat %>%
  filter(psabs %in% summary$psabs) %>%
  group_by(nitrogen, psabs, tax_group)

save(top_loci_traits_all_snps, file="cache/top_loci_traits_all_snps.rda")


top_loci_traits <- plot_dat %>%
  filter(psabs %in% summary$psabs) %>%
  group_by(nitrogen, psabs, tax_group) %>%
  filter(log10p == max(log10p))

save(top_loci_traits, file="cache/top_loci_traits.rda")

unique(top_loci_traits$trait)





### do yield analysis for traits

### check for significant differences







### get bins with signals in both N treatments

stdN_bins <- unique(filter(short_snp_dat, nitrogen == "+N")$psabs)

lowN_bins <- unique(filter(short_snp_dat, nitrogen == "-N")$psabs)

common_bins <- Reduce(intersect, list(stdN_bins,lowN_bins))

# common bins are 81168 107174 107204
# different taxa in -N/+N in the common bins!

common_bin_data <- short_snp_dat %>%
  filter(psabs %in% common_bins)

common_bin_taxa <- plot_dat %>%
  filter(psabs %in% common_bins) %>%
  dplyr::select(nitrogen, psabs, tax_group, total_sign_snps, sign_snps) %>%
  unique()



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

ggplot(threshold_plot_data, aes(x=total_sign_snps, y=taxa_retained)) +
  geom_point() +
  geom_vline(xintercept = 23, color = "red", linetype = "dashed") +
  theme_bw()


ggplot(threshold_plot_data, aes(x=total_sign_snps, y=bins_retained)) +
  geom_point() +
  geom_vline(xintercept = 23, color = "red", linetype = "dashed") +
  theme_bw()



#### zoom into interesting signals

## microbes dominating signals:

top_bins_lowN <- c(153820, 33428, 21332, 72012, 69707, 127620)
top_bins_stdN <- c(15907, 161129, 73876, 188680, 67065, 41433)

### get top taxa in interesting bins:

taxa_dat_lowN <- plot_dat %>%
  filter(nitrogen == "-N" & psabs %in% top_bins_lowN) %>%
  dplyr::select(nitrogen, psabs, tax_group, total_sign_snps, sign_snps) %>%
  unique() %>%
  group_by(psabs) %>%
  filter(sign_snps == max(sign_snps))


taxa_dat_stdN <- plot_dat %>%
  filter(nitrogen == "+N" & psabs %in% c(15907, 161129, 73876, 188680, 67065, 41433)) %>%
  dplyr::select(nitrogen, psabs, tax_group, total_sign_snps, sign_snps) %>%
  unique() %>%
  group_by(psabs) %>%
  filter(sign_snps == max(sign_snps))


## use SNP data of top taxa, zoom into bin region

load("largedata/combined_SNP_data.rda")
load("data/data_summary_150_traits.rda")


plot_dat_all_snps <- gwas_dat %>%
  #dplyr::select(-chindex, -psabs) %>%
  mutate(log10p = -log10(p_wald)) %>%
  left_join(bpadd, by = "chr") %>%
  mutate(psabs = chindex + ps) %>%
  mutate(psabs = as.numeric(gsub('.{4}$', '', psabs))) %>% ## pos. in 10kb units
  mutate(bin = as.numeric(gsub(paste0(".{",4,"}$"), '', ps)) + 1) %>% # absolute bin position
  group_by(chr, bin, nitrogen) %>%
  add_count(count = sum(-log10(p_wald) >= 5), name="total_sign_snps") %>%
  left_join(dplyr::select(data_summary_150_traits, trait, tax_group)) %>%
  group_by(chr, bin, nitrogen, trait) %>%
  add_count(count = sum(-log10(p_wald) >= 5), name="sign_snps")






### LD analysis



### plot zoom

plotZoom <- function(){
  
  
  colors = c("1" = "#276FBF", "3" = "#276FBF", "5" = "#276FBF", "7" = "#276FBF", "9" = "#276FBF",
             "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059", "sig" = "#ff0000")
  
  zoom_plot <-  ggplot(zoom_dat, aes(x = ps, y = log10p, color = chr)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = 5, color="red", linetype="dashed") +
    labs(x = NULL, y = "-log(p_wald)") +
    scale_x_continuous(breaks = axis_range, labels= c(min(zoom_dat$bin)-1, zoom_range)) +
    theme_bw() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)
    )
  return(zoom_plot)
  
}
  
  
### find representative snp

repSNP <- function(){
  rep_snp <- zoom_dat %>%
    filter(bin == zoom_bin) %>%
    arrange(-log10p) %>%
    dplyr::select(snp_id, log10p)
  return(rep_snp)
}



### plot hapmap




nitr = "-N"



zoom_taxon <- "Filimonas sp 2"
zoom_bin <- 21333
zoom_chr <- 1
gwas <- read_delim('largedata/top_12_traits/T75.assoc.txt', delim="\t")
hmp <- read_delim('largedata/Hapmaps/hapmap_chr1.txt', delim="\t", na = c("", "NA", "NN"))



zoom_taxon <- "Oryzihumus terrae"
zoom_bin <- 2731
zoom_chr <- 2
gwas <- read_delim('largedata/top_12_traits/T110.assoc.txt', delim="\t")
hmp <- read_delim('largedata/Hapmaps/hapmap_chr2.txt', delim="\t", na = c("", "NA", "NN"))



zoom_taxon = "Mesorhizobium huakuii"
zoom_bin = 14567
zoom_chr <- 3
gwas <- read_delim('largedata/top_12_traits/T30.assoc.txt', delim="\t")
hmp <- read_delim('largedata/Hapmaps/hapmap_chr3.txt', delim="\t", na = c("", "NA", "NN"))



zoom_taxon = "Massilia niabensis"
zoom_bin = 16872
zoom_chr <- 3
gwas <- read_delim('largedata/top_12_traits/T130.assoc.txt', delim="\t")
hmp <- read_delim('largedata/Hapmaps/hapmap_chr3.txt', delim="\t", na = c("", "NA", "NN"))



zoom_taxon = "Rhizobium daejeonense"
zoom_bin = 1846
zoom_chr <- 6
gwas <- read_delim('largedata/top_12_traits/T34.assoc.txt', delim="\t")
hmp <- read_delim('largedata/Hapmaps/hapmap_chr6.txt', delim="\t", na = c("", "NA", "NN"))



zoom_taxon = "Niabella yanshanensis"
zoom_bin = 10692
zoom_chr <- 7
gwas <- read_delim('largedata/top_12_traits/T83.assoc.txt', delim="\t")
hmp <- read_delim('largedata/Hapmaps/hapmap_chr7.txt', delim="\t", na = c("", "NA", "NN"))






nitr = "+N"

zoom_taxon = "Archangium gephyra"
zoom_bin = 15908
zoom_chr <- 1
gwas <- read_delim('largedata/top_12_traits/T37.assoc.txt', delim="\t")
hmp <- read_delim('largedata/Hapmaps/hapmap_chr1.txt', delim="\t", na = c("", "NA", "NN"))

zoom_taxon = "Stenotrophomonas maltophilia"
zoom_bin = 10736
zoom_chr <- 2
gwas <- read_delim('largedata/top_12_traits/T16.assoc.txt', delim="\t")
hmp <- read_delim('largedata/Hapmaps/hapmap_chr2.txt', delim="\t", na = c("", "NA", "NN"))

zoom_taxon = "Massilia putida"
zoom_bin = 11925
zoom_chr <- 3
gwas <- read_delim('largedata/top_12_traits/T129.assoc.txt', delim="\t")
hmp <- read_delim('largedata/Hapmaps/hapmap_chr3.txt', delim="\t", na = c("", "NA", "NN"))

zoom_taxon = "Acinetobacter nosocomialis"
zoom_bin = 18736
zoom_chr <- 3
gwas <- read_delim('largedata/top_12_traits/T6.assoc.txt', delim="\t")
hmp <- read_delim('largedata/Hapmaps/hapmap_chr3.txt', delim="\t", na = c("", "NA", "NN"))

zoom_taxon = "Pseudomonas umsongensis"
zoom_bin = 18001
zoom_chr <- 7
gwas <- read_delim('largedata/top_12_traits/T118.assoc.txt', delim="\t")
hmp <- read_delim('largedata/Hapmaps/hapmap_chr7.txt', delim="\t", na = c("", "NA", "NN"))

zoom_taxon = "f_A21b"
zoom_bin = 9275
zoom_chr <- 9
gwas <- read_delim('largedata/top_12_traits/T128.assoc.txt', delim="\t")
hmp <- read_delim('largedata/Hapmaps/hapmap_chr9.txt', delim="\t", na = c("", "NA", "NN"))




#### SET VALUES


zoom_from <- zoom_bin-15
zoom_to <- zoom_bin+15

zoom_dat <- gwas %>%
  filter(chr == zoom_chr) %>%
  mutate(log10p = -log10(p_wald)) %>%
  left_join(bpadd, by = "chr") %>%
  mutate(psabs = chindex + ps) %>%
  mutate(psabs = as.numeric(gsub('.{4}$', '', psabs))) %>% ## pos. in 10kb units
  mutate(bin = as.numeric(gsub(paste0(".{",4,"}$"), '', ps)) + 1) %>% # absolute bin position
  filter(bin >= zoom_from) %>%
  filter(bin <= zoom_to) %>%
  mutate(snp_id = paste0("S", chr, "_", ps)) %>%
  mutate(chr = as.factor(chr))


zoom_range <- seq( min(zoom_dat$bin), max(zoom_dat$bin))
zoom_ps_min <- as.numeric(gsub(paste0(".{",4,"}$"), '', min(zoom_dat$ps)))*10000
zoom_ps_max <- (as.numeric(gsub(paste0(".{",4,"}$"), '', max(zoom_dat$ps)))+1)*10000
axis_range <- seq(zoom_ps_min, zoom_ps_max, 10000)

#### SET VALUES


plotZoom()



##### DRAW LD PLOT

colnames(hmp)[1] <- "snp_id"

## randomly select up to n snps per bin for LD plot
set.seed(2021)
LD_dat <- zoom_dat %>%
  sample_n(200) #%>%
  #group_by(bin) %>%
  #sample_n(zoom_sample_n, replace = TRUE) %>%
  #unique()
bins <- zoom_range
#print(LD_snps$snp_id)



LD_haplotypes <- hmp %>%
  filter(snp_id %in% LD_dat$snp_id)

#selected_haplotypes <- LD_haplotypes %>%
  #filter(snp_id == 0)

#for (bin in bins) {
  #start <- as.numeric(paste0(as.character(bin-1), "0000"))
  #end <- as.numeric(paste0(as.character(bin), "0000"))
  #add_haplotypes <- LD_haplotypes %>%
    #filter(pos > start & pos <= end)
  #selected_haplotypes <- rbind(selected_haplotypes, add_haplotypes)
#}

#colnames(selected_haplotypes)


## make filename
binrange <- ifelse(length(bins) > 1, paste0(min(bins), "-", max(bins)), max(bins))
ntreat <- ifelse(nitr == "+N", "stdN", "lowN")
fname <- paste(str_replace(zoom_taxon," ","_"), ntreat, binrange, sep = "_")
fpath <- paste0("figures/LD_plots/", fname, ".png")

#LD_haplotypes <- data.frame(LD_haplotypes)
#rownames(LD_haplotypes) <- LD_haplotypes$snp_id
gene <- data.frame(t(LD_haplotypes[,12:ncol(LD_haplotypes)]))
gty <- makeGenotypes(gene,sep="")

rgb.palette <- colorRampPalette(rev(c("blue","orange" ,"red")), space = "rgb")

png(file=fpath,res=300,width=1000,height=1000)
myld <- LDheatmap(gty,genetic.distances=LD_haplotypes$pos,flip=TRUE, text=FALSE, color=rgb.palette(20), title=fname)
dev.off()




##### DRAW LD PLOT



repSNP()

### find representative snp





