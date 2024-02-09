getwd()
library("dplyr")


## loading data
data_meta <- read.csv("2_metadata_species_list.csv")
data_comp <- read.csv("2_allsites_sp.comp_clean.csv")

## site comp, and metadata colnames, verify which colnames to merge by
colnames(data_comp)
colnames(data_meta)

## issue with title. data_meta, "code_field_updated" should be "code"
colnames(data_meta)[5] <- "code"


## merging dataset by code, site
full_sp.comp <- full_join(data_meta, data_comp, 
                  by = c("site", "taxon", "code"), 
                  keep = TRUE,
                  relationship = "many-to-many")

## looking at data
write.csv(full_sp.comp, file = "03_spcomp_full.csv")

