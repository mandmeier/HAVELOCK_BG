library("phyloseq")
library("tidyverse")

load("cache/ps_core5t.rda")

### agglomerate ASV counts for each taxonomic group

ps_glom <- tax_glom(ps_core5t, taxrank="tax_group")

save(ps_glom, file = "cache/ps_glom.rda")
