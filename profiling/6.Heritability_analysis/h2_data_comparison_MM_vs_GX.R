#h2_data_comparison_MM_vs_GX.R




load("data/group_data.rda")



data_h2_HN <- read_delim(file="data/Blup_and_Heritability/Heritability_HN_result.txt", delim = "\t")

data_h2_HN <- data_h2_HN %>%
  rename(ASV=Trait, h2_gx=Heritability) %>%
  mutate(nitrogen="HN")


data_h2_LN <- read_delim(file="data/Blup_and_Heritability/Heritability_LN_result.txt", delim = "\t")
data_h2_LN <- data_h2_LN %>%
  rename(ASV=Trait, h2_gx=Heritability) %>%
  mutate(nitrogen="LN")

gx <- rbind(data_h2_HN, data_h2_LN)



h2_HN <- group_data %>%
  dplyr::select(ASV, tax_group, H2_stdN_19) %>%
  rename(h2_mm=H2_stdN_19) %>%
  mutate(nitrogen="HN")


h2_LN <- group_data %>%
  dplyr::select(ASV, tax_group, H2_lowN_19) %>%
  rename(h2_mm=H2_lowN_19) %>%
  mutate(nitrogen="LN")


mm <- rbind(h2_HN, h2_LN)


H2_data <- left_join(mm, gx)




ggplot(H2_data, aes(x=h2_mm, y=h2_gx)) +
  facet_wrap(~nitrogen, nrow = 1) +
  ggtitle("Heritability for each microbe trait") +
  geom_point()



