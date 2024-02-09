
## 1: adding in species comp data per site, and meta data. 
comp_metadata <- read.csv("form-1__mk-nutnet2023-5dominantplants(edits).csv")
comp_sevi <- read.csv("branch-4__species-composition-sevi(edits).csv")
comp_temple <- read.csv("branch-1__species-composition-temple(edits).csv")
comp_lubb <- read.csv("branch-2__species-composition-lubb(edits).csv")
comp_arch <- read.csv("branch-3__species-composition-arch(edits).csv")


## checking col. names
colnames(comp_arch)
colnames(comp_lubb)
colnames(comp_sevi)
colnames(comp_temple)



## 2: Merging datasets
# loading the packages. 
install.packages("dplyr")
install.packages("tidyverse")
library("dplyr")
library("tidyverse")

# all sites together, using rbind
comp_all <- rbind(comp_arch, comp_lubb, comp_sevi, comp_temple)

# Checking that observation #'s match the comp_all observaitons!!
71 + 62 +84 +144





# renaming the metadata similar colname
# ec5_branch_uuid" -> "ec5_uuid"

full_join(comp_metadata, comp_all_sites, by = "ec5_uuid",)
full_join()
? "full_join"



## trying to merge the datasets again 
## following epicollect guide
# https://docs.epicollect.net/common-use-cases/consolidate-data#using-r

test <- full_join(comp_lubb, comp_arch, by = "ec5_branch_uuid")


colnames(comp_arch)
