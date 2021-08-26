#Yield_traits_GWAS_zoom.R

library("LDheatmap")
library("genetics")




# chr 4 bin 20790



### plot zoom

plotZoom <- function(){
  
  
  colors = c("1" = "#276FBF", "3" = "#276FBF", "5" = "#276FBF", "7" = "#276FBF", "9" = "#276FBF",
             "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059", "sig" = "#ff0000")
  
  
  
  colors = c("1" = "#fcae1e", "3" = "#fcae1e", "5" = "#fcae1e", "7" = "#fcae1e", "9" = "#fcae1e",
             "2" = "#dd571c", "4" = "#dd571c", "6" = "#dd571c", "8" = "#dd571c", "10" = "#dd571c", "sig" = "#ff0000")
  

  
  zoom_plot <-  ggplot(zoom_dat, aes(x = ps, y = log10p, color = chr)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = 3, color="red", linetype="dashed") +
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




nitr = "+N"
zoom_taxon <- "Pajaroellobacter"
#zoom_bin <- 20732 # "Pajaroellobacter"
zoom_bin <- 20790 # "Pajaroellobacter"
#zoom_bin <- 20791 # "Pajaroellobacter"
zoom_chr <- 4
gwas <- read_delim('largedata/GWAS/short_1_out_traits_150_stdN_201109-103909/T38.assoc.txt', delim="\t")

hmp <- read_delim('largedata/Hapmaps/hapmap_chr4.txt', delim="\t", na = c("", "NA", "NN"))


## ExG, GLI, RGB
gwas <- read_delim('largedata/Yield_Traits_GWAS/HN-ExG.assoc.txt', delim="\t")





load(file="cache/Agro_traits_GWAS/microbe_gwas_overview.rda")


load(file="cache/Yield_traits_GWAS/yield_gwas_overview.rda")




test <- filter(microbe_gwas_overview, bin == zoom_bin)

unique(test$psabs)

test2 <- filter(yield_gwas_overview, psabs %in% unique(test$psabs))


table(test2$trait)


#### SET VALUES


zoom_from <- zoom_bin-10
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




plotZoom()