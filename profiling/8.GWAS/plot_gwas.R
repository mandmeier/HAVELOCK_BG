library("tidyverse")
library("ape")

load("data/group_data.rda")
gwas_dat_stdN <- read_csv("data/merged_short_5_out_traits_150_stdN_201109-103909.csv")
gwas_dat_lowN <- read_csv("data/merged_short_5_out_traits_150_lowN_201109-100251.csv")


gwas_dat_stdN <- read_tsv("temp/T118.assoc.stdN.txt")
gwas_dat_stdN$trait <- "T118"
gwas_dat_lowN <- read_tsv("temp/T118.assoc.lowN.txt")
gwas_dat_lowN$trait <- "T118"


gwas_dat_stdN <- read_tsv("temp/T34.assoc.stdN.txt")
gwas_dat_stdN$trait <- "T34"
gwas_dat_lowN <- read_tsv("temp/T34.assoc.lowN.txt")
gwas_dat_lowN$trait <- "T34"



### combine both N treatments
gwas_dat_stdN$nitrogen <- "+N"
gwas_dat_lowN$nitrogen <- "-N"

gwas_dat <- rbind(gwas_dat_stdN, gwas_dat_lowN)


#### prepare gwas table for plotting ####

### add annotation data
gwas_dat <- left_join(gwas_dat, group_data, by = "trait")




nCHR <- length(unique(gwas_dat$chr))

## get chromosome lengths (counted with grep)
## grep "1-" T88.assoc.txt | sort -n -k 3| tail
#chr_len <- c(306971061, 289663740, 235651490, 288208778, 284473431, 218627434, 186522559, 253629336, 279123320, 296807830)

## get chromosome lengths (counted from range of data)
#test <- gwas_dat %>%
#  filter(chr == 10)
#range(test$ps)

chr_len <- c(306970668, 244420837, 235651057, 246967222, 223706090, 173536910, 181718119, 181046068, 159686117, 150929986)
vlines <- c(0,cumsum(chr_len)) ## draw vertical lines at chomosome ends



chindex <- c(0) ## this number needs to be added to position within each chromosome
for (i in c(1:nCHR)){
  print(i)
  chindex[i+1] <- sum(chr_len[1:i])
}
chindex <- chindex[1:10]
bpadd <- data.frame("chr" = c(1:nCHR), "chindex" = chindex)
gwas_dat <- left_join(gwas_dat, bpadd, by = "chr")


# calculate absolute snp position
gwas_dat$psabs <- gwas_dat$ps + gwas_dat$chindex

## find position for labels
axis_set <- gwas_dat %>% 
  group_by(chr) %>% 
  summarize(center = (max(psabs) + min(psabs)) / 2) # find middle of chromosome (to place labels)



#save(gwas_dat, file = "largedata/combined_SNP_data.rda")


#### annotate significant peaks ####

## add 10kb bins
#gwas_dat <- gwas_dat %>%
  #mutate(bin = as.numeric(gsub('.{4}$', '', ps)) + 1)

#### function to draw manhattan plot ####

draw_manhattan <- function(gwas_dat, tax_grp, binsize = 4, count = 10,chromosomes = c(1:10)){
  
  binsize <- 4
  chromosomes <- c(1:10)
  count <- 10
  #tax_grp <- "Streptomyces sp 1"
  #tax_grp <- "Sphingobium herbicidovorans 2"
  #tax_grp <- "Pseudomonas umsongensis"
  #tax_grp <- "Rhizobium daejeonense"
  bin_regex <- paste0(".{",binsize,"}$")
  plot_dat <- gwas_dat %>%
    filter(tax_group == tax_grp) %>%
    filter(-log10(p_wald) >= 1) %>%
    mutate(bin = as.numeric(gsub(bin_regex, '', ps)) + 1) %>%
    filter(chr %in% chromosomes) %>%
    group_by(chr, bin, nitrogen, trait) %>%
    add_count(count = sum(-log10(p_wald) >= 5))
    #add_tally(name="count")
  

  
  ## mark significant snips red if counts are > count threshold
  plot_dat$color_group <- as.character(plot_dat$chr)
  plot_dat$color_group <- ifelse(plot_dat$count >= count, "sig", plot_dat$color_group)
  
  
  
  if (length(unique(plot_dat$nitrogen)) == 0){
    dummy <- gwas_dat[c(1,2),]
    dummy$nitrogen <- c("+N", "-N")
    dummy$p_wald <- NA
    plot_dat$color_group <- as.character( plot_dat$color_group)
    plot_dat <- rbind(plot_dat, dummy)
  }
  
  ### if one N treatment is missing add dummy value
  if (length(unique(plot_dat$nitrogen)) == 1){
    ntreat <- c("+N", "-N")
    # add dummy row
    dummy <- plot_dat[1,]
    dummy$nitrogen <- ntreat[unique(plot_dat$nitrogen) != ntreat]
    dummy$p_wald <- NA
    plot_dat <- rbind(plot_dat, dummy)
  }
  
  ## find range for y axis
  ylim <- abs(floor(log10(min(plot_dat$p_wald, na.rm = TRUE)))) + 1
  #print(ylim)
  
  colors = c("1" = "#276FBF", "3" = "#276FBF", "5" = "#276FBF", "7" = "#276FBF", "9" = "#276FBF",
    "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059", "sig" = "#ff0000")
  
  
  manhplot <- ggplot(plot_dat, aes(x = psabs, y = -log10(p_wald), color = color_group)) +
    geom_point(alpha = 0.75) +
    #geom_hline(yintercept = -log10(5), color = "grey40", linetype = "dashed") + 
    #scale_x_continuous(label = axis_set$chr, breaks = axis_set$center, limits = c(1273000000, 1273300000)) +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(1, ylim)) +
    geom_hline(yintercept=5, linetype="dashed", color = "red") +
    #geom_vline(xintercept = vlines) +
    scale_color_manual(values = colors) +
    scale_size_continuous(range = c(0.5,3)) +
    facet_wrap(~nitrogen, nrow = 2) +
    labs(x = NULL, 
      y = "-log10(p)") + 
    ggtitle(paste0(tax_grp, " | binsize = ", 10^binsize/1000, "kb | marked if >= ", count, " sig. snps per bin")) +
    theme_bw() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
    )
  
  return(manhplot)
  
}


draw_manhattan(gwas_dat, "Pseudomonas umsongensis")



draw_manhattan(gwas_dat, "Bryobacter")
draw_manhattan(gwas_dat, "Niabella yanshanensis")
draw_manhattan(gwas_dat, "Rhizobium daejeonense")
draw_manhattan(gwas_dat, "Pseudomonas umsongensis")
draw_manhattan(gwas_dat, "Acidothermus")
draw_manhattan(gwas_dat, "Acinetobacter nosocomialis")
draw_manhattan(gwas_dat, "Candidatus Udaeobacter copiosus")
draw_manhattan(gwas_dat, "Filimonas sp 2")
draw_manhattan(gwas_dat, "Massilia niabensis")
draw_manhattan(gwas_dat, "Microbacterium testaceum")
draw_manhattan(gwas_dat, "Sphingoaurantiacus")
draw_manhattan(gwas_dat, "Sphingobium herbicidovorans 1")
draw_manhattan(gwas_dat, "Taibaiella chishuiensis")
draw_manhattan(gwas_dat, "TM7a sp 1")



sort(unique(group_data$tax_group))
sort(unique(gwas_dat$tax_group))



sort(unique(group_data$tax_group))

#### print plots to pdf for each tax group ####

# Make list of variable names to loop over.
var_list = sort(unique(group_data$tax_group))

# Make plots.
plot_list = list()
for (i in 1:length(var_list)) {
  p = draw_manhattan(gwas_dat, var_list[i])
  plot_list[[i]] = p
}

# Another option: create pdf where each page is a separate plot.
pdf("largedata/gwas_10kb_10counts_test.pdf")
for (i in 1:length(plot_list)) {
  print(var_list[i])
  print(plot_list[[i]])
}
dev.off()





#### print plots to pdf for each of top 28 tax groups ####

load("cache/top_taxa_bins_genes.rda")

top_28_taxa <- sort(unique(top_bins$tax_group))

# Make list of variable names to loop over.
var_list = sort(unique(group_data$tax_group))


var_list = var_list[var_list %in% top_28_taxa]

# Make plots.
plot_list = list()
for (i in 1:length(var_list)) {
  p = draw_manhattan(gwas_dat, var_list[i])
  plot_list[[i]] = p
}

# Another option: create pdf where each page is a separate plot.
pdf("largedata/gwas_10kb_10counts_top28.pdf")
for (i in 1:length(plot_list)) {
  print(var_list[i])
  print(plot_list[[i]])
}
dev.off()






#### plot peaks associated with multiple traits ####




draw_group_manhattan <- function(gwas_dat, tax_groups, nitr = "+N", binsize = 4, count = 10,chromosomes = c(1:10), zoom = NA){
  
  binsize <- 4
  chromosomes <- c(6)
  count <- 10
  #tax_groups <- c("Oryzihumus terrae", "f_Vicinamibacteraceae", "Candidatus Nitrocosmicus oleophilus", "Parafilimonas")
  tax_groups <- c("Rhizobium daejeonense")
  nitr = "-N"
  #min = 0
  #max = 20000000
  #tax_grp <- "Streptomyces sp 1"
  #tax_grp <- "Sphingobium herbicidovorans 2"
  #tax_grp <- "Niabella yanshanensis"
  bin_regex <- paste0(".{",binsize,"}$")
  plot_dat <- gwas_dat %>%
    filter(tax_group %in% tax_groups) %>%
    filter(nitrogen == nitr) %>%
    #filter(ps > min & ps < max) %>%
    mutate(bin = as.numeric(gsub(bin_regex, '', ps)) + 1) %>%
    filter(chr %in% chromosomes) %>%
    group_by(chr, bin, nitrogen, trait) %>%
    add_tally(name="count")
  
  ### if zoom
  
  
  
  zoom <- c(6, 1868)
  
  
  
  ## mark significant snips red if counts are > count threshold
  plot_dat$color_group <- as.character(plot_dat$chr)
  plot_dat$color_group <- ifelse(plot_dat$count >= count, "sig", plot_dat$color_group)
  
  
  
  if (length(unique(plot_dat$nitrogen)) == 0){
    dummy <- gwas_dat[c(1,2),]
    dummy$nitrogen <- c("+N", "-N")
    dummy$p_wald <- NA
    plot_dat$color_group <- as.character( plot_dat$color_group)
    plot_dat <- rbind(plot_dat, dummy)
  }
  
  ### if one N treatment is missing add dummy value
  if (length(unique(plot_dat$nitrogen)) == 1){
    ntreat <- c("+N", "-N")
    # add dummy row
    dummy <- plot_dat[1,]
    dummy$nitrogen <- ntreat[unique(plot_dat$nitrogen) != ntreat]
    dummy$p_wald <- NA
    plot_dat <- rbind(plot_dat, dummy)
  }
  
  ## find range for y axis
  ylim <- abs(floor(log10(min(plot_dat$p_wald, na.rm = TRUE)))) + 1
  #print(ylim)
  
  colors = c("1" = "#276FBF", "3" = "#276FBF", "5" = "#276FBF", "7" = "#276FBF", "9" = "#276FBF",
    "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059", "sig" = "#ff0000")
  
  
  manhplot <- ggplot(plot_dat, aes(x = psabs, y = -log10(p_wald), color = color_group)) +
    geom_point(alpha = 0.75) +
    #geom_hline(yintercept = -log10(5), color = "grey40", linetype = "dashed") + 
    #scale_x_continuous(label = axis_set$chr, breaks = axis_set$center, limits = c(1273000000, 1273300000)) +
    #scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(5, ylim)) +
    scale_color_manual(values = colors) +
    scale_size_continuous(range = c(0.5,3)) +
    facet_wrap(~tax_group, nrow = length(tax_groups)) +
    labs(x = NULL, 
      y = "-log10(p)") + 
    ggtitle(paste0(nitr, " | binsize = ", 10^binsize/1000, "kb | marked if >= ", count, " sig. snps per bin")) +
    theme_bw() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
    )
  
  
  if(!is.na(zoom)){
    abs <- round(mean(filter(plot_dat, chr == as.character(zoom[1]) & bin == zoom[2])$psabs), -binsize)
    bin_len = 10^binsize
    min <- (abs) - 6*bin_len
    if(min < 0){min <- 0}
    max <- (abs) + 1*bin_len
    
    manhplot <- manhplot  +
      scale_x_continuous(limits = c(min, max), breaks=seq(min,max,bin_len))
      
  }
  
  
  
  return(manhplot)
  
}







