getwd()
setwd("C:/Users/monik/Documents/reserach/fieldwork_summer_2023/data/
      2_processing/isotope_weight_values/3_isotope_rolling_20240103")


# Version 2 = made the colnames more R friendly
data_multi <- read.csv("multi_sample_library_v2.0.csv")
data_meta <- read.csv("species_list.csv")

colnames(data_multi)

?"read.csv"


# Installing packages
install.packages("dplyr")
install.packages("tidyverse")
library("dplyr")
library("tidyverse")

test <- merge.data.frame(data_meta, data_multi, by.x = "taxon", all = TRUE)

? "merge.data.frame"
test <- full_join(data_meta, data_multi, by = "code")

          
write.csv(test, file = "test.csv")
