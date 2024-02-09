getwd()

try.data <- read.csv(file = "TryAccspecies.csv")


library(tidyverse)
library(dplyr)

# Isolating rows with lat, long, N & C values only.
trydata_ofintrest <- try.data %>% filter(
    DataName == "Longitude" | 
      DataName == "Latitude" | 
      DataName == "Leaf nitrogen content per dry mass (Nmass)" |
      DataName == "Leaf carbon content per dry mass"
  )



# Writing the isolated data frame to a csv
write.csv(trydata_ofintrest, file = "try_cleaned_dataofintrest.csv", 
          row.names = FALSE)



unique(try.data$SpeciesName)
length(unique(try.data$SpeciesName))
length(unique(trydata_ofintrest$SpeciesName))


pivot(trydata_ofintrest, trydata_ofintrest$AccSpeciesID, 
      trydata_ofintrest$TraitName)







## Chat GPT ask

library(tidyr)
library(dplyr)


# Identify and create unique IDs for each group of related rows
datatest <- trydata_ofintrest %>%
  mutate(DataName = cumsum(!is.na(Latitude) | !is.na(Longitude)))


# Reshape the data to have one row per observation
reshaped_data <- your_data %>%
  group_by(GroupID) %>%
  fill(c(Latitude, Longitude), .direction = "downup") %>%
  filter(row_number() == n() | is.na(Latitude)) %>%
  ungroup() %>%
  select(-GroupID)