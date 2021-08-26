# plot_gwas_v5


library("tidyverse")
library("stringr")
library("gridExtra")


load("largedata/combined_SNP_data.rda")

load("largedata/GWAS2/plot_dat_v5_gt5.rda")

load("data/group_data.rda")


sort(unique(plot_dat$tax_group))


draw_manhattan <- function(taxgrp, nitr="HN", threshold_line = 5) {
  
  #snp_dat <- read_tsv("largedata/top_10_traits/T75.assoc.txt")
  #snp_dat <- read_tsv("largedata/top_10_traits/T6.assoc.txt")
  #taxgrp <- "Acinetobacter nosocomialis"
  #taxgrp <- "Acinetobacter johnsonii"
  
  #chromosomes <- c(1:10)
  #snp_dat <- snps
  #threshold_line <- 5
  #taxgrp <- "Pseudomonas chlororaphis"
  #nitr <- "HN"
  print(nitr)
  title <- paste(taxgrp, "|", nitr)
  
  
  plot_obj <- plot_dat %>%
    filter(tax_group == taxgrp) #%>%
    #filter(nitrogen == nitr) 
  
  chr_len <- c(308452471, 243675191, 238017767, 250330460, 226353449, 181357234, 185808916, 182411202, 163004744, 152435371)
  
  #vlines <- c(0,cumsum(chr_len)) ## draw vertical lines at chomosome ends
  
  ## find position for labels
  axis_set <- plot_obj %>% 
    group_by(chr) %>% 
    summarize(center = (max(psabs) + min(psabs)) / 2) # find middle of chromosome (to place labels)
  
  
  ## find range for y axis
  ylim <- floor(max(plot_obj$log10p, na.rm = TRUE)) + 1
  
  
  colors = c("1" = "#58CCED", "3" = "#58CCED", "5" = "#58CCED", "7" = "#58CCED", "9" = "#58CCED",
             "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059", "sig" = "#ff0000")
  
  
  manhplot <- ggplot(plot_obj, aes(x = psabs, y = log10p, color=chr)) +
    geom_point(alpha = 0.75) +
    #geom_hline(yintercept = -log10(5), color = "grey40", linetype = "dashed") + 
    #scale_x_continuous(label = axis_set$chr, breaks = axis_set$center, limits = c(1273000000, 1273300000)) +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(2, ylim), breaks = c(2:ylim)) +
    geom_hline(yintercept=threshold_line, linetype="dashed", color = "red") +
    #geom_vline(xintercept = vlines) +
    scale_color_manual(values = colors) +
    scale_size_continuous(range = c(0.5,3)) +
    facet_wrap(~nitrogen, nrow=2, scales = "free_y") +
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
  
  return(manhplot)
  
}




##draw_manhattan(taxgrp = "Chryseobacterium sp 2", nitr <- "LN") this microbe was excluded


sort(unique(plot_dat$tax_group))


plots <- list()

taxa <- sort(total)
taxa <- taxa[taxa %in% unique(plot_dat$tax_group)]


for (i in c(1:length(taxa))){
  tax = taxa[i]
  print(tax)
  plot <- draw_manhattan(taxgrp = tax)
  plots[[i]] <- list(plot)
}



pdf("largedata/GWAS2/p5_manhattan_CC.pdf", onefile = TRUE)

for (i in seq(length(plots))) {
  do.call("grid.arrange", plots[[i]])  
}
dev.off()





draw_manhattan(taxgrp = "Ellin6067", nitr = "LN")
draw_manhattan(taxgrp = "Pseudomonas chlororaphis", nitr = "LN")
draw_manhattan(taxgrp = "Massilia niabensis", nitr = "LN") #### NICE
draw_manhattan(taxgrp = "Niabella yanshanensis", nitr = "LN")
draw_manhattan(taxgrp = "Pelomonas saccharophila", nitr = "LN") #### NICE
draw_manhattan(taxgrp = "Steroidobacter", nitr = "LN") #### NICE
draw_manhattan(taxgrp = "Conexibacter", nitr = "LN") #### NICE
draw_manhattan(taxgrp = "f_Comamonadaceae Unknown Genus", nitr = "LN")


draw_manhattan(taxgrp = "Archangium gephyra", nitr = "HN") #### NICE
draw_manhattan(taxgrp = "Variovorax paradoxus", nitr = "HN")
draw_manhattan(taxgrp = "Sphingoaurantiacus", nitr = "HN") #### NICE
draw_manhattan(taxgrp = "Massilia putida", nitr = "HN")
draw_manhattan(taxgrp = "Nitrospira japonica", nitr = "HN")
draw_manhattan(taxgrp = "Acinetobacter nosocomialis", nitr = "HN") #### NICE
draw_manhattan(taxgrp = "Acidothermus", nitr = "HN") #### NICE
draw_manhattan(taxgrp = "Salmonella bongori", nitr = "HN")
draw_manhattan(taxgrp = "Rhizobium phaseoli-multihospitium-favelukesii", nitr = "HN")
draw_manhattan(taxgrp = "f_A21b", nitr = "HN") #### NICE










draw_manhattan(taxgrp = "Rhizobium daejeonense", nitr = "HN")


v





draw_manhattan(taxgrp = "Acinetobacter nosocomialis", nitr <- "HN")









