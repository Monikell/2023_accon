getwd()
setwd("C:/Users/monik/Documents/reserach/fieldwork_summer_2023/data/
      2_processing/isotope_weight_values/evan_tim/tim_sevi")


## Installing packages
install.packages("plyr")
library("plyr")
install.packages("dplyr")
library("dplyr")


# loading sevi trait data, thank you Tim! 
data_raw <- read.csv(file = "sevi_plant_traits.csv")


# Keeping the columns that I want.All the meta data and the c and n values. 
## columns 1 - 25
data_relevant <- data_raw[c(1:25)]
colnames(data_relevant)


# removing columns with NA values
data_noNA <- na.omit(data_relevant)


# Added a new column called "taxon" that combines genus and species 
data_noNA$taxon <- paste(data_noNA$genus_new, data_noNA$species_new, sep = " ")

colnames(data_noNA)

# Getting the mean C per species
carbon_mean_sp <- ddply(data_noNA, .(genus_new, species_new), 
                       transform, mean.c.percent_taxon = 
                         mean(lf_iso_pC))

# Adding the mean N to the C mean per species to make one document
cn_means_sp <- ddply(carbon_mean_sp, .(genus_new, species_new),
                        transform, mean.n.value_taxon = 
                       mean(lf_iso_pN))


# seems to have acheived the goal! Writing to a csv. 
write.csv(cn_means_sp, file = "sevi_cleaned_data.csv")






