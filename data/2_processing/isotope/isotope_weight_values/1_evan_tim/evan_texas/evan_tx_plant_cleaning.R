getwd()

install.packages("plyr")
install.packages("dplyr")
install.packages("ggplot2")

library(plyr)
library(dplyr)
library(ggplot2)

# data enter
evan_data <- read.csv(file = "TXeco_data_CN rates for plants.csv")

# see column names
colnames(evan_data)


# Dataset ID, Observation and AccSp. ID and name, lat long and N and C values
data_clean <- evan_data[c(1:4,6,8:13,15:17)]
colnames(data_clean)


# inserting taxon information
taxon_info <- read.csv(file = "species_code_taxon.csv")



# Adding taxon informaiton? 
test <- merge(taxon_info, data_clean, all = TRUE)

?"merge"


# Getting the mean N and C for each species and site:
c_means_site <- ddply(test, .(site, NCRS.code), 
                 transform, mean.c.value_persite = mean(c.leaf))

# Getting means for N, and combining it with "c_means_site"
cn_means_site <- ddply(c_means_site, .(site, NCRS.code),
                  transform, mean.n.value_persite = mean(n.leaf))

colnames(cn_means_site)

# Getting the mean C per species 
c_means_taxon <- ddply(cn_means_site, .(NCRS.code), 
                      transform, mean.c.value_taxon = mean(c.leaf))

# Adding the mean N to the C mean per species to make one document
cn_means_taxon <- ddply(c_means_taxon, .(NCRS.code),
                       transform, mean.n.value_taxon = mean(n.leaf))

colnames(cn_means_taxon)


# I think it worked, writing to csv. 
write.csv(cn_means_taxon, file = "evan_tx_data_cn_means_mkelley.csv")
