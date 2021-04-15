library("tidyverse")
library("stringr")

# plot GWAS of all microbial traits

load("largedata/combined_SNP_data.rda")

load("data/group_data.rda")

# find bins with at least 27 sig snps
bin_regex <- paste0(".{",4,"}$")

# add bins
plot_dat <- gwas_dat %>%
  select(-X1) %>%
  #filter(trait %in% top_traits$trait) %>%
  #filter(-log10(p_wald) >= 1) %>%
  mutate(bin = as.numeric(gsub(bin_regex, '', ps)) + 1) %>%
  group_by(chr, bin, nitrogen, trait) %>%
  add_count(count = sum(-log10(p_wald) >= 5)) %>%
  left_join(select(group_data, trait, tax_group)) %>%
  filter(count >= 10) %>%
  mutate(log10p = -log10(p_wald)) %>%
  mutate(chrbin = paste0(chr, "_", bin))


draw_manhattan <- function(snp_dat, taxgrp, nitr="+N", chromosomes = c(1:10), threshold_line = 5) {
  
  #snp_dat <- read_tsv("largedata/top_10_traits/T75.assoc.txt")
  #snp_dat <- read_tsv("largedata/top_10_traits/T6.assoc.txt")
  #taxgrp <- "Acinetobacter nosocomialis"
  #chromosomes <- c(1:10)
  #threshold_line <- 5
  
  #snp_dat <- snps
  #taxgrp <- tgr
  
  
  sign_snps <- plot_dat %>%
    filter(tax_group == taxgrp) %>%
    filter(nitrogen == nitr) %>%
    filter(n >= 27)
  
  title <- paste(taxgrp, "|", nitr)
  
  
  plot_obj <- snp_dat %>%
    select(rs, chr, ps, p_wald) %>%
    mutate(log10p = -log10(p_wald)) %>%
    filter(log10p >= 2) %>%
    select(-p_wald)
  
  
  chr_len <- c(306970668, 244420837, 235651057, 246967222, 223706090, 173536910, 181718119, 181046068, 159686117, 150929986)
  #vlines <- c(0,cumsum(chr_len)) ## draw vertical lines at chomosome ends
  
  nCHR <- length(unique(snp_dat$chr))
  
  chindex <- c(0) ## this number needs to be added to position within each chromosome
  for (i in c(1:nCHR)){
    #print(i)
    chindex[i+1] <- sum(chr_len[1:i])
  }
  chindex <- chindex[1:10]
  bpadd <- data.frame("chr" = c(1:nCHR), "chindex" = chindex)
  
  
  plot_obj <- plot_obj %>%
    left_join(bpadd, by = "chr") %>%
    mutate(psabs = ps + chindex) # calculate absolute snp position
  
  
  ## find position for labels
  axis_set <- plot_obj %>% 
    group_by(chr) %>% 
    summarize(center = (max(psabs) + min(psabs)) / 2) # find middle of chromosome (to place labels)
  
  
  ## subset chromosomes
  plot_obj <- plot_obj %>%
    filter(chr %in% chromosomes)
  
  # find bins with at least 27 sig snps
  bin_regex <- paste0(".{",4,"}$")
  
  plot_obj <- plot_obj %>%
    mutate(bin = as.numeric(gsub(bin_regex, '', ps)) + 1) %>%
    mutate(chrbin = paste0(chr, "_", bin)) %>%
    mutate(color_group = ifelse(chrbin %in% unique(sign_snps$chrbin), "sig", chr))
  
  
  ## find range for y axis
  ylim <- floor(max(plot_obj$log10p, na.rm = TRUE)) + 1
  
  
  colors = c("1" = "#58CCED", "3" = "#58CCED", "5" = "#58CCED", "7" = "#58CCED", "9" = "#58CCED",
             "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059", "sig" = "#ff0000")
  
 
  manhplot <- ggplot(plot_obj, aes(x = psabs, y = log10p, color = color_group)) +
    geom_point(alpha = 0.75) +
    #geom_hline(yintercept = -log10(5), color = "grey40", linetype = "dashed") + 
    #scale_x_continuous(label = axis_set$chr, breaks = axis_set$center, limits = c(1273000000, 1273300000)) +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(2, ylim), breaks = c(2:ylim)) +
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
      panel.grid.minor.x = element_blank()#,
      #axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
    )
    
  
  #manhplot
  
  return(manhplot)
  
}



plot_manh <- function(tgr, N){
  #tgr <- "Acinetobacter nosocomialis"
  trait <- filter(group_data, tax_group == tgr)$trait
  infile <- paste0("largedata/GWAS/",ifelse(N == "+N", "short_1_out_traits_150_stdN_201109-103909", "short_1_out_traits_150_lowN_201109-100251"),"/", trait, ".assoc.txt")
  print(paste(N, "plotting", tgr))
  if(file.exists(infile)){
    snps <- read_tsv(infile)
    mplot <- draw_manhattan(snps, tgr, nitr=N, chromosomes = c(1:10))
    outfile <- paste0("largedata/manhattan_plots/",ifelse(N == "+N", "stdN", "lowN"),"/", str_replace_all(tgr," ","_"),".png")
    ggsave(filename=outfile, plot=mplot, width = 8, height = 4)
  } else {
    print(paste(infile, "not found"))
  }
}


## draw all 10 Manhattan plots

for (N in c("+N", "-N")){
  for (tgr in group_data$tax_group){
    plot_manh(tgr, N)
  }
}










