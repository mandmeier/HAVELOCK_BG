##mediation_analysis


library("tidyverse")



mediators_stdN <- read_delim("data/Mediation_analysis/mediators_microbiome_stdN.txt", delim="\t")
mediators_lowN <- read_delim("data/Mediation_analysis/mediators_microbiome_lowN.txt", delim="\t")
load("data/data_summary_150_traits.rda")
names(mediators_stdN) <- names(mediators_lowN)

mediators <- rbind(mediators_stdN, mediators_lowN)


## find index, make compatible with 10kbin data


#chr_len <- c(306970668, 244420837, 235651057, 246967222, 223706090, 173536910, 181718119, 181046068, 159686117, 150929986)
chr_len_10k <- c(306980000, 244430000, 235660000, 246970000, 223710000, 173540000, 181720000, 181050000, 159690000, 150930000)

nCHR <- length(unique(gwas_dat$chr))

chindex <- c(0) ## this number needs to be added to position within each chromosome
for (i in c(1:nCHR)){
  #print(i)
  chindex[i+1] <- sum(chr_len_10k[1:i])
}



chindex <- chindex[1:10]
bpadd <- data.frame("chr" = as.character(c(1:nCHR)), "chindex" = chindex)


root_mediators <- mediators %>%
  mutate(trait = strsplit(trait, "_")[[1]][1]) %>%
  filter(tissue == "GRoot") %>%
  left_join(bpadd, by = "chr") %>%
  mutate(psabs = pos + chindex) %>%
  left_join(select(data_summary_150_traits, trait, tax_group)) %>%
  mutate(nitrogen = ifelse(treatment == "stdN", "+N", "-N")) %>%
  select(-Freq, -treatment) %>%
  group_by(nitrogen, id_v4) %>%
  add_tally(name="Freq")
  




plot_dat <- root_mediators %>%
  select(nitrogen,Freq, chr, psabs, tax_group, id_v4) %>%
  unique()



colors = c("1" = "#276FBF", "3" = "#276FBF", "5" = "#276FBF", "7" = "#276FBF", "9" = "#276FBF",
           "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059", "sig" = "#ff0000")


## find position for labels
axis_set <- plot_dat %>% 
  group_by(chr) %>% 
  summarize(center = (max(psabs) + min(psabs)) / 2) # find middle of chromosome (to place labels)



=
mediation_plot <-  ggplot(plot_dat, aes(x = psabs, y = Freq, color = chr)) +
  #geom_point() +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~nitrogen, nrow = 2) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_color_manual(values = colors) +
  labs(x = NULL, y = "Number of mediated traits") + 
  geom_vline(xintercept =  c(0, 2104600000), color = "red", size=0) + # to combine plots
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

mediation_plot



  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_color_manual(values = colors) +
  scale_size_continuous(breaks = c(1, 5, 10)) +
  #scale_size_continuous(range = c(0.5,3)) +
  facet_wrap(~nitrogen, nrow = 2) +
  labs(x = NULL, y = "Total sign. SNPs per 10kb bin") + 
  #ggtitle("total signif. SNPs for 399 10kb bins with nearby annotated genes") +
  theme_bw() +
  theme( 
    #legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

overview_plot







