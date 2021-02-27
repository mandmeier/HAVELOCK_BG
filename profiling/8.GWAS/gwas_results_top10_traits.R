# plot GWAS of top 10 microbial traits

# bins were selected if there were at least 27 significant snps above p value 5.
# as determined by looking at bins that map to known genes.

library("tidyverse")


## select top 10 microbial traits

microbial_traits <- read_csv("data/microbial_traits.csv")

top_traits <- microbial_traits %>%
  filter(GWAS_signal == "strong") %>%
  mutate(N_treatment = ifelse(N_treatment == "stdN", "+N", "-N")) %>%
  rename(nitrogen = N_treatment)


## get top 10kb bins of top 10 microbes


load("largedata/combined_SNP_data.rda")

# find bins with at least 27 sig snps
bin_regex <- paste0(".{",4,"}$")

# add bins
plot_dat <- gwas_dat %>%
  select(-X1) %>%
  filter(trait %in% top_traits$trait) %>%
  #filter(-log10(p_wald) >= 1) %>%
  mutate(bin = as.numeric(gsub(bin_regex, '', ps)) + 1) %>%
  group_by(chr, bin, nitrogen, trait) %>%
  add_count(count = sum(-log10(p_wald) >= 5)) %>%
  filter(count >= 10) %>%
  left_join(select(top_traits, trait, tax_group, nitrogen, GWAS_signal)) %>%
  filter(GWAS_signal == "strong") %>%
  mutate(log10p = -log10(p_wald)) %>%
  mutate(chrbin = paste0(chr, "_", bin))
  
#length(unique(plot_dat$rs))




### manhattan plots





draw_manhattan <- function(snp_dat, taxgrp, chromosomes = c(1:10), threshold_line = 5, zoom = "none") {
  
  #snp_dat <- read_tsv("largedata/top_10_traits/T75.assoc.txt")
  
  #snp_dat <- read_tsv("largedata/top_10_traits/T6.assoc.txt")
  #taxgrp <- "Acinetobacter nosocomialis"
  #chromosomes <- c(1:10)
  #threshold_line <- 5
  
  sign_snps <- plot_dat %>%
    filter(tax_group == taxgrp)
  
  title <- paste(taxgrp, "|", sign_snps$nitrogen[1])
  
  
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
  
  if(is.list(zoom)) {
    plot_obj <- plot_obj %>%
      filter(chr == zoom[1]) %>%
      filter(ps >= zoom[2] & ps <= zoom[3])
  }
  
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
  
  if(is.list(zoom)) {
    
    chindex <- plot_obj$chindex[1]
    lbl <- paste0(c(seq(unlist(zoom[2])/1000, unlist(zoom[3])/1000, 10)), "k")
    brks <- c(seq(unlist(zoom[2])+chindex, unlist(zoom[3])+chindex, 10000))
    lms <- c(unlist(zoom[2])+chindex, unlist(zoom[3])+chindex)
    
    manhplot <- ggplot(plot_obj, aes(x = psabs, y = log10p, color = color_group)) +
      geom_point(alpha = 0.75) +
      scale_x_continuous(label = lbl,  breaks = brks,  limits = lms) +
      scale_y_continuous(expand = c(0,0), limits = c(2, ylim), breaks = c(2:ylim)) +
      scale_color_manual(values = colors) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = paste("chr", zoom[1]), y = "-log10(p)") + 
      theme_bw() +
      theme( 
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
      )
    
  } else {
    
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
  
  }
  
  return(manhplot)
  
}

library("stringr")


plot_manh <- function(tgr){
  #tgr <- "Acinetobacter nosocomialis"
  trait <- filter(top_traits, tax_group == tgr)$trait
  infile <- paste0("largedata/top_10_traits/", trait, ".assoc.txt")
  snps <- read_tsv(infile)
  mplot <- draw_manhattan(snps, tgr, chromosomes = c(1:10))
  outfile <- paste0("figures/Manhattan_plots_top10/", str_replace_all(tgr," ","_"),".png")
  ggsave(filename=outfile, plot=mplot, width = 8, height = 4)
}


## draw all 10 Manhattan plots

for (tgr in top_traits$tax_group){
  print(paste("plotting", tgr))
  plot_manh(tgr)
}





### zoom in to region around SNPS







plot_zoom <- function(taxgroup, snps, bins, interval = 10, save = FALSE){
  
  ## get chrom
  chrm <- str_split(bins[1], "_")[[1]][1]
  ## get bin positions
  binpos <- bins %>%
    map(function(x) str_split(x, "_")[[1]][2]) %>%
    unlist(use.names=FALSE) %>%
    as.numeric()
  
  ## get middle position
  middle_pos <- (binpos[1]-1) + (binpos[length(binpos)]-(binpos[1]-1))/2
  lower <- (middle_pos-interval)*10000
  upper <- (middle_pos+interval)*10000
  zoom <- list(chrm, lower, upper)
  
  mplot <- draw_manhattan(snps, taxgroup, zoom = zoom)
  
  if (save == TRUE) {
    outfile <- paste0("figures/Manhattan_plots_top10/", str_replace_all(taxgroup," ","_"),"_chr",chrm,"_",lower,"-",upper,".png")
    ggsave(filename=outfile, plot=mplot, width = 8, height = 2)
  }

  return(mplot)
  
}




tgr <- "Acinetobacter nosocomialis"
tgr_snps <- read_tsv("largedata/top_10_traits/T6.assoc.txt")
tgr_top_snps <- plot_dat %>%
  filter(tax_group == tgr)
## bins
tgr_top_bins <- unique(tgr_top_snps$chrbin)

mplot <- plot_zoom(tgr, tgr_snps, tgr_top_bins, save = TRUE)
mplot



tgr <- "Candidatus Udaeobacter copiosus"
tgr_snps <- read_tsv("largedata/top_10_traits/T54.assoc.txt")
tgr_top_snps <- plot_dat %>%
  filter(tax_group == tgr)
## bins
tgr_top_bins <- unique(tgr_top_snps$chrbin)

mplot <- plot_zoom(tgr, tgr_snps, tgr_top_bins, save = TRUE)
mplot



tgr <- "f_A21b"
tgr_snps <- read_tsv("largedata/top_10_traits/T128.assoc.txt")
tgr_top_snps <- plot_dat %>%
  filter(tax_group == tgr)
## bins
tgr_top_bins <- unique(tgr_top_snps$chrbin)

chr3_loc1 <- c("3_2692", "3_2693")
mplot <- plot_zoom(tgr, tgr_snps, chr3_loc1, save = TRUE)
mplot

chr3_loc2 <- c("3_17437", "3_17438", "3_17441", "3_17442")
mplot <- plot_zoom(tgr, tgr_snps, chr3_loc2, save = TRUE)
mplot


chr9_loc1 <- c("9_9224", "9_9225", "9_9229", "9_9231")
mplot <- plot_zoom(tgr, tgr_snps, chr9_loc1, save = TRUE)
mplot

chr9_loc2 <- c("9_9253", "9_9256")
mplot <- plot_zoom(tgr, tgr_snps, chr9_loc2, save = TRUE)
mplot

chr9_loc3 <- c("9_9274", "9_9277")
mplot <- plot_zoom(tgr, tgr_snps, chr9_loc3, save = TRUE)
mplot

chr9_loc4 <- c("9_13916", "9_13917")
mplot <- plot_zoom(tgr, tgr_snps, chr9_loc4, save = TRUE)
mplot


## larger window chr 9
chr9 <- c("9_9220", "9_9285")
mplot <- plot_zoom(tgr, tgr_snps, chr9, interval=35, save = TRUE)
mplot



tgr <- "f_Comamonadaceae Unknown Genus"
tgr_snps <- read_tsv("largedata/top_10_traits/T145.assoc.txt")
tgr_top_snps <- plot_dat %>%
  filter(tax_group == tgr)
## bins
tgr_top_bins <- unique(tgr_top_snps$chrbin)


chr2_loc1 <- c("2_21744", "2_21745")
mplot <- plot_zoom(tgr, tgr_snps, chr2_loc1, save = TRUE)
mplot

chr5_loc1 <- c("5_14749", "5_14750", "5_14753", "5_14754", "5_14755", "5_14756")
mplot <- plot_zoom(tgr, tgr_snps, chr5_loc1, save = TRUE)
mplot




tgr <- "Filimonas sp 2"
tgr_snps <- read_tsv("largedata/top_10_traits/T75.assoc.txt")
tgr_top_snps <- plot_dat %>%
  filter(tax_group == tgr)
## bins
tgr_top_bins <- unique(tgr_top_snps$chrbin)


chr1_loc1 <- c("1_1277", "1_1278", "1_1280", "1_1282", "1_1286")
mplot <- plot_zoom(tgr, tgr_snps, chr1_loc1, save = TRUE)
mplot

chr1_loc2 <- c("1_21325", "1_21326", "1_21327", "1_21331", "1_21332", "1_21333", "1_21334", "1_21335", "1_21338", "1_21340", "1_21344")
mplot <- plot_zoom(tgr, tgr_snps, chr1_loc2, save = TRUE)
mplot

chr4_loc1 <- c("4_16096", "4_16097")
mplot <- plot_zoom(tgr, tgr_snps, chr4_loc1, save = TRUE)
mplot

## left some out "5_2686"  "5_2691" "5_3122"

chr5_loc3 <- c("5_3579", "5_3580", "5_3585", "5_3592", "5_3600", "5_3608")
mplot <- plot_zoom(tgr, tgr_snps, chr5_loc3, save = TRUE)
mplot

chr5_loc4 <- c("5_3757", "5_3769", "5_3771", "5_3772", "5_3773", "5_3774")
mplot <- plot_zoom(tgr, tgr_snps, chr5_loc4, save = TRUE)
mplot




tgr <- "Ilumatobacter"
tgr_snps <- read_tsv("largedata/top_10_traits/T97.assoc.txt")
tgr_top_snps <- plot_dat %>%
  filter(tax_group == tgr)
## bins
tgr_top_bins <- unique(tgr_top_snps$chrbin)


chr2_loc1 <- c("2_10610", "2_10649")
mplot <- plot_zoom(tgr, tgr_snps, chr2_loc1, save = TRUE)
mplot

chr5_loc1 <- c("5_8326", "5_8333")
mplot <- plot_zoom(tgr, tgr_snps, chr5_loc1, save = TRUE)
mplot

chr5_loc2 <- c("5_8343", "5_8346")
mplot <- plot_zoom(tgr, tgr_snps, chr5_loc2, save = TRUE)
mplot

chr5_loc3 <- c("5_8363", "5_8390")
mplot <- plot_zoom(tgr, tgr_snps, chr5_loc3, save = TRUE)
mplot


## larger window chr 5
chr5 <- c("5_8323", "5_8386")
mplot <- plot_zoom(tgr, tgr_snps, chr5, interval=35, save = TRUE)
mplot



tgr <- "Massilia niabensis"
tgr_snps <- read_tsv("largedata/top_10_traits/T130.assoc.txt")
tgr_top_snps <- plot_dat %>%
  filter(tax_group == tgr)
## bins
tgr_top_bins <- unique(tgr_top_snps$chrbin)

mplot <- plot_zoom(tgr, tgr_snps, tgr_top_bins, save = TRUE)
mplot


tgr <- "Niabella yanshanensis"
tgr_snps <- read_tsv("largedata/top_10_traits/T83.assoc.txt")
tgr_top_snps <- plot_dat %>%
  filter(tax_group == tgr)
## bins
tgr_top_bins <- unique(tgr_top_snps$chrbin)


chr8_loc1 <- c("8_11943", "8_11970")
mplot <- plot_zoom(tgr, tgr_snps, chr8_loc1, interval = 15, save = TRUE)
mplot





tgr <- "Rhizobium daejeonense"
tgr_snps <- read_tsv("largedata/top_10_traits/T34.assoc.txt")
tgr_top_snps <- plot_dat %>%
  filter(tax_group == tgr)
## bins
tgr_top_bins <- unique(tgr_top_snps$chrbin)


chr4_loc1 <- c("4_21403", "4_21404")
mplot <- plot_zoom(tgr, tgr_snps, chr4_loc1, save = TRUE)
mplot

chr6_loc1 <- c("6_1818", "6_1835")
mplot <- plot_zoom(tgr, tgr_snps, chr6_loc1, save = TRUE)
mplot

chr6_loc2 <- c("6_1839", "6_1855")
mplot <- plot_zoom(tgr, tgr_snps, chr6_loc2, save = TRUE)
mplot

chr6_loc3 <- c("6_1863", "6_1868")
mplot <- plot_zoom(tgr, tgr_snps, chr6_loc3, save = TRUE)
mplot

## omitted chr 7 and 10 "7_10195" "7_10217" "10_6676" "10_8723"


## larger window chr 6
chr6 <- c("6_1816", "6_1875")
mplot <- plot_zoom(tgr, tgr_snps, chr6, interval=35, save = TRUE)
mplot



tgr <- "Sphingobium herbicidovorans 1"
tgr_snps <- read_tsv("largedata/top_10_traits/T28.assoc.txt")
tgr_top_snps <- plot_dat %>%
  filter(tax_group == tgr)
## bins
tgr_top_bins <- unique(tgr_top_snps$chrbin)


chr2_loc1 <- c("2_6979", "2_6980")
mplot <- plot_zoom(tgr, tgr_snps, chr2_loc1, save = TRUE)
mplot

chr3_loc1 <- c("3_310", "3_311")
mplot <- plot_zoom(tgr, tgr_snps, chr3_loc1, save = TRUE)
mplot

chr8_loc1 <- c("8_14527", "8_14530")
mplot <- plot_zoom(tgr, tgr_snps, chr8_loc1, save = TRUE)
mplot

chr9_loc1 <- c("9_1163", "9_1164")
mplot <- plot_zoom(tgr, tgr_snps, chr9_loc1, save = TRUE)
mplot










### find intervals for LD plots for all candidate loci

loci <- list(
  list(c(3, 187340, 187420)),
  list(c(9, 103530, 103590))
)


names(loci) <- c("Acinetobacter nosocomialis", "Candidatus Udaeobacter copiosus")

loci




loci$`Acinetobacter nosocomialis`[[1]][2]


