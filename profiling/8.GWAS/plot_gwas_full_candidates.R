#plot_gwas_full



library("tidyverse")
library("stringr")
library("gridExtra")

tgr <- "Sphingoaurantiacus"
asv <- "asv_000281"
gwas_dat <- read_delim("largedata/GWAS2/top13/HN/asv_000281.assoc.txt", delim = "\t")
nit <- "HN"
#2-230389923 pos 236826599

#test <- strsplit(gwas_dat$rs, "-")

#test <- sapply(strsplit(gwas_dat$rs, split='-', fixed=TRUE), `[`, 2)



tgr <- "Acinetobacter nosocomialis"
asv <- "asv_000095"
gwas_dat <- read_delim("largedata/GWAS2/top13/HN/asv_000095.assoc.txt", delim = "\t")
nit <- "HN"
#3-184434608



tgr <- "Steroidobacter"
asv <- "asv_001993"
gwas_dat <- read_delim("largedata/GWAS2/top13/LN/asv_001993.assoc.txt", delim = "\t")
nit <- "LN"
#5-183779458 pos 188470211


tgr <- "Conexibacter"
asv <- "asv_001234"
gwas_dat <- read_delim("largedata/GWAS2/top13/LN/asv_001234.assoc.txt", delim = "\t")
nit <- "LN"
#6-167618482 pos 171561709



### fix snp position ### not necessary, believe the position, snp naming is different!!
#gwas_dat$ps <- as.numeric(sapply(strsplit(gwas_dat$rs, split='-', fixed=TRUE), `[`, 2))



load("data/data_summary_150_traits.rda")



chr_len_10k <- c(308460000, 243680000, 238020000, 250340000, 226360000, 181360000, 185810000, 182420000, 163010000, 152440000)




nCHR <- 10
chindex <- c(0) ## this number needs to be added to position within each chromosome
for (i in c(1:nCHR)){
  #print(i)
  chindex[i+1] <- sum(chr_len_10k[1:i])
}
chindex <- chindex[1:10]
bpadd <- data.frame("chr" = c(1:nCHR), "chindex" = chindex)



options(scipen = 999)








draw_manhattan <- function(gwas_dat, taxgrp, nitr="HN", threshold_line = 5) {
  
  #snp_dat <- read_tsv("largedata/top_10_traits/T75.assoc.txt")
  #snp_dat <- read_tsv("largedata/top_10_traits/T6.assoc.txt")
  #taxgrp <- "Acinetobacter nosocomialis"
  #taxgrp <- "Acinetobacter johnsonii"
  
  #chromosomes <- c(1:10)
  #snp_dat <- snps
  threshold_line <- 5
  #taxgrp <- "Pseudomonas chlororaphis"
  #nitr <- "HN"
  
  taxgrp <- tgr
  nitr <- nit
  title <- paste(taxgrp, "|", nitr)
  
  #chr_len <- c(308452471, 243675191, 238017767, 250330460, 226353449, 181357234, 185808916, 182411202, 163004744, 152435371)
  
  
  
  # add bins
  plot_obj <- gwas_dat %>%
    mutate(log10p = -log10(p_wald)) %>%
    #filter(log10p >= 5) %>%
    left_join(bpadd, by = "chr") %>%
    mutate(psabs = chindex + ps) %>%
    #mutate(psabs = as.numeric(gsub('.{4}$', '', psabs))) %>% ## pos. in 10kb units
    mutate(ps = as.character(ps)) %>%
    mutate(bin = as.numeric(gsub('.{4}$', '', ps)) + 1) %>% # absolute bin position in 10kb units
    group_by(chr, bin) %>%
    add_count(count = sum(-log10(p_wald) >= 5), name="total_sign_snps") %>%
    ##mutate(total_mean_p = mean(log10p)) %>%
    group_by(chr, bin) %>%
    add_count(count = sum(-log10(p_wald) >= 5), name="sign_snps") %>%
    ##mutate(mean_p = mean(log10p)) %>%
    mutate(ASV=asv) %>%
    left_join(dplyr::select(data_summary_150_traits, ASV, tax_group)) %>%
    #mutate(color = ifelse(chr %in% c(1,3,5,7,9), "#276FBF", "#183059")) %>%
    mutate(chr = as.character(chr)) 
  
  test <- plot_obj %>%
    filter(chr =="6")
  
  
  
  #vlines <- c(0,cumsum(chr_len)) ## draw vertical lines at chomosome ends
  
  ## find position for labels
  axis_set <- plot_obj %>% 
    group_by(chr) %>% 
    summarize(center = (max(psabs) + min(psabs)) / 2) # find middle of chromosome (to place labels)
  
  
  ## find range for y axis
  ylim <- floor(max(plot_obj$log10p, na.rm = TRUE)) + 1
  
  
  colors = c("1" = "#58CCED", "3" = "#58CCED", "5" = "#58CCED", "7" = "#58CCED", "9" = "#58CCED",
             "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059", "sig" = "#ff0000")
  
  
  
  
  
  chdx <- chindex[6]
  
  #171560000

  xmin <- 171410000 + chdx
  xmax <- 171710000 + chdx
  breaks <- seq(from=xmin, to=xmax, by=10000)
  
  labels <- seq(xmin-chdx, xmax-chdx, 10000 ) /1000
  

  manhplot <- ggplot(plot_obj, aes(x = psabs, y = log10p, color=chr)) +
    geom_point(alpha = 0.75) +
    #geom_hline(yintercept = -log10(5), color = "grey40", linetype = "dashed") + 
    scale_x_continuous(limits = c(xmin, xmax), breaks = breaks, labels=labels) +
    #scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(1, ylim), breaks = c(1:ylim)) +
    geom_hline(yintercept=threshold_line, linetype="dashed", color = "red") +
    #geom_vline(xintercept = vlines) +
    scale_color_manual(values = colors) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, y = "-log10(p)") + 
    ggtitle(title) +
    theme_bw() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
    )
  
  manhplot
  
  
  return(manhplot)
  
}

p <- draw_manhattan(plotobj, tgr, nit)



