library("tidyverse")



load("data/group_data.rda")



data_blup_HN <- read_delim(file="data/Blup_and_Heritability/data_blup_HN_result.txt", delim = "\t")

unique(data_blup_HN$id) # 210

data_blup_HN <- data_blup_HN %>%
  drop_na() %>%
  pivot_longer(cols=c(-id),names_to="trait", values_to = "abundance") %>%
  group_by(trait) %>%
  summarize(mean=mean(abundance), sd=sd(abundance)) %>%
  rename(ASV = trait, mean_gx = mean, sd_gx = sd)

HN <- left_join(data_blup_HN, dplyr::select(group_data, ASV, trait))



blup_stdN <- read_csv("data/blup_stdN_150_tax_groups.csv")

unique(blup_stdN$genotype) # 226

blup_stdN <- blup_stdN %>%
  pivot_longer(cols=c(-genotype),names_to="trait", values_to = "abundance")%>%
  group_by(trait) %>%
  summarize(mean=mean(abundance), sd=sd(abundance)) %>%
  rename(mean_mm = mean, sd_mm = sd)

HN <- left_join(HN, blup_stdN)
HN$nitrogen <- "HN"





data_blup_LN <- read_delim(file="data/Blup_and_Heritability/data_blup_LN_result.txt", delim = "\t")


unique(data_blup_LN$id) # 209


data_blup_LN <- data_blup_LN %>%
  drop_na() %>%
  pivot_longer(cols=c(-id),names_to="trait", values_to = "abundance") %>%
  group_by(trait) %>%
  summarize(mean=mean(abundance), sd=sd(abundance)) %>%
  rename(ASV = trait, mean_gx = mean, sd_gx = sd)

LN <- left_join(data_blup_LN, dplyr::select(group_data, ASV, trait))

blup_lowN <- read_csv("data/blup_lowN_150_tax_groups.csv")

unique(blup_lowN$genotype) #221

blup_lowN <- blup_lowN %>%
  drop_na() %>%
  pivot_longer(cols=c(-genotype),names_to="trait", values_to = "abundance")%>%
  group_by(trait) %>%
  summarize(mean=mean(abundance), sd=sd(abundance)) %>%
  rename(mean_mm = mean, sd_mm = sd)

LN <- left_join(LN, blup_lowN)
LN$nitrogen <- "LN"

BLUPS <- rbind(HN, LN)


ggplot(BLUPS, aes(x=mean_mm, y=mean_gx)) +
  facet_wrap(~nitrogen, nrow = 1) +
  ggtitle("mean BLUPs for each trait") +
  geom_point()













hist(data_blup_HN)
