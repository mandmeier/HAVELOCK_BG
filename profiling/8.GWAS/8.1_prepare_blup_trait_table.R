#### prepare trait matrix for GWAS ####

library("phyloseq")
library("tidyverse")


## add unused genotypes as NA
names <- read.csv("data/BG_MM_Gen_names.csv")
## get in right order. this matches the order in data/PCA_277.csv
order <- c("33-16", "38-11Goodman-Buckler", "4226", "4722", "A6", "A188", "A214NGoodman-Buckler", "A239", "A441-5", "A554", "A556", "A619", "A632", "A634", "A635", "A641", "A654", "A659", "A661", "A679", "A680", "A682", "Ab28A", "B2", "B10", "B14A", "B37", "B46", "B52", "B57", "B64", "B68", "B73", "B73Htrhm", "B75", "B76", "B77", "B79", "B84", "B97", "B103", "B104", "B105", "B109", "B164", "C49A", "C103", "C123", "CH9", "CH701-30", "CI3A", "Ci7Goodman-Buckler", "CI21E", "CI28AGoodman-Buckler", "CI31A", "CI64", "CI66", "Ci90CGoodman-Buckler", "Ci91BGoodman-Buckler", "CI187-2", "CM7", "CM37", "CM105", "CM174", "CML5", "CML10", "CML11", "CML14", "CML38", "CML45", "CML52", "CML61", "CML69", "CML77", "CML91", "CML92", "CML103", "CML108", "CML154Q", "CML157Q", "CML158Q", "CML218", "CML220", "CML228", "CML238", "CML247", "CML254", "CML258", "CML261", "CML264", "CML277", "CML281", "CML287", "CML311", "CML314", "CML321", "CML322", "CML323", "CML328", "CML331", "CML332", "CML333", "CML341", "CMV3", "CO106", "CO125", "CO255", "D940Y", "DE-2", "DE1", "DE811", "E2558W", "EP1", "F6", "F7", "F44", "F2834T", "GA209", "GT112", "H49", "H84", "H91", "H95", "H99", "H105W", "Hi27Goodman-Buckler", "HP301", "I29", "I137TN", "I205", "i1677a", "IA2132Goodman-Buckler", "Ia5125", "IDS28", "IDS69", "IDS91", "Il14H", "Il101", "ILLHy", "K4", "K55", "K64", "K148", "Ki3", "Ki11", "Ki14", "Ki21", "Ki43", "Ki44", "Ki2021", "Ky21", "KY226", "KY228", "L317", "L578", "M14", "M37W", "M162W", "MEF156-55-2", "Mo1W", "Mo17", "Mo18W", "Mo24W", "Mo44", "Mo45", "Mo46", "Mo47", "MoG", "Mp339", "MS71", "MS153", "MS1334", "Mt42", "N6", "N7A", "N28Ht", "N192", "NC33", "NC222", "NC230", "NC232", "NC236", "NC238", "NC250", "NC258", "NC260", "NC262", "NC264", "NC290A", "NC294", "NC296", "NC296A", "NC298", "NC300", "NC302", "NC304", "NC306", "NC310", "NC314", "NC318", "NC320", "NC324", "NC326", "NC328", "NC336", "NC338", "NC340", "NC342", "NC344", "NC346", "NC348", "NC350", "NC352", "NC354", "NC356", "NC358", "NC360", "NC362", "NC364", "NC366", "NC368", "ND246", "Oh7B", "Oh40B", "Oh43", "Oh43E", "Oh603", "Os420", "P39Goodman-Buckler", "Pa91", "Pa762", "Pa875", "Pa880", "R4", "R109B", "R168", "R177", "R229", "SA24", "SC55", "SC213R", "SC357", "SD40", "SD44", "Sg18", "Sg1533", "T8", "T232", "T234", "Tx303", "Tx601", "Tzi8", "Tzi9", "Tzi10", "Tzi11", "Tzi16", "Tzi18", "Tzi25", "U267Y", "Va14", "Va17", "Va22", "Va26", "Va35", "Va59", "Va85", "Va99", "VA102", "VaW6", "W22", "W64A", "W117Ht", "W153R", "W182B", "WD", "Wf9", "Yu796")

### 150 tax groups, stdN
blup <- read_csv("data/blup_stdN_150_tax_groups.csv") %>%
  rename(MM_name = genotype)

### use Gen Xu's nomenclature of maize genotypes
blup2 <-  names %>%
  dplyr::select(MM_name, GX_name) %>%
  filter(GX_name != "" & MM_name != "") %>%
  left_join(blup) %>%
  dplyr::select(-MM_name)

## get in same format as PCA file, i.e. 277 genotypes
blup2 <- blup2[match(order, blup2$GX_name),]
blup2 <- blup2[, -1]


## remember which trait is which ASV
#traits <- data.frame("ASV" = colnames(blup2), "trait" = paste0("T", c(1:ncol(blup2))))
## add to group data
#load("data/group_data.rda")
#group_data <- left_join(group_data, traits)
#save(group_data, file = "data/group_data.rda")

### save trait data
write.table(blup2, file = "cache/GWAS/traits_150_stdN.txt", sep = "\t", col.names = FALSE, row.names = FALSE)



### 150 tax groups, lowN
blup <- read_csv("data/blup_lowN_150_tax_groups.csv") %>%
  rename(MM_name = genotype)

### use Gen Xu's nomenclature of maize genotypes
blup2 <-  names %>%
  dplyr::select(MM_name, GX_name) %>%
  filter(GX_name != "" & MM_name != "") %>%
  left_join(blup) %>%
  dplyr::select(-MM_name)

## get in same format as PCA file, i.e. 277 genotypes
blup2 <- blup2[match(order, blup2$GX_name),]
blup2 <- blup2[, -1]

### save trait data
write.table(blup2, file = "cache/GWAS/traits_150_lowN.txt", sep = "\t", col.names = FALSE, row.names = FALSE)


test <- read.table("cache/GWAS/traits_150_lowN.txt")




