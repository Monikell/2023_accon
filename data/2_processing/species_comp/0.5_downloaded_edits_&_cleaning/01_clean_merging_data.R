
## 1: adding in species comp data per site, and meta data. 
comp_metadata <- read.csv("form-1__mk-nutnet2023-5dominantplants.csv")
comp_sevi <- read.csv("branch-4__species-composition-sevi.csv")
comp_temple <- read.csv("branch-1__species-composition-temple.csv")
comp_lubb <- read.csv("branch-2__species-composition-lubb.csv")
comp_arch <- read.csv("branch-3__species-composition-arch.csv")


# merging meta data and species comp. 
? "merge"



## 2: Merging datasets
# packages
install.packages("dplyr")
install.packages("tidyverse")
library("dplyr")
library("tidyverse")

# all sites together.. missing like 
comp_all_sites <- bind_rows(comp_arch, comp_lubb, comp_sevi, comp_lubb)

# observation of sites math, #shows 279 observations
71+62+84+144
361 - 279
# missing like 82 observaitons not sure why??

colnames(comp_all_sites)
colnames(comp_metadata)


# renaming the metadata similar colname
# ec5_branch_uuid" -> "ec5_uuid"

full_join(comp_metadata, comp_all_sites, by = "ec5_uuid",)
full_join()
? "full_join"



## trying to merge the datasets again 
## following epicollect guide
# https://docs.epicollect.net/common-use-cases/consolidate-data#using-r


