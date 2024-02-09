
# following along with the rtry vignettes: 
# https://cran.r-project.org/web/packages/rtry/vignettes/rtry-introduction.html

####################
### Inserting data into environment, mergining datasets into 1 
####################

getwd()

# Loading packages
install.packages(c("data.table", "dplyr", "tidyr", "jsonlite", "curl"))
install.packages("rtry")
library("rtry")
library("dplyr")
library("tidyr")

# Adding the data from TRY into the environment
data1 <- rtry_import(input = '30770.txt')
data2 <- rtry_import(input = '30769.txt')
data3 <- rtry_import(input = '30771.txt')

# Verifying col. names match before merging the datasets
d1names <- colnames(data1)
d2names <- colnames(data2)
d3names <- colnames(data3)

# identical vectors? 
identical(d1names, d2names)
identical(d2names, d3names)


# Merging the data into 1 massive data set
data <- rbind(data1, data2, data3)
colnames(data)


####################
### Cleaning data using rtry
####################


## Following along with ?"rtry_trans_wider" example 

# 1. Select only the trait records that have standardized numeric values.
num_traits <- rtry_select_row(data, complete.cases(TraitID) 
                              & complete.cases(StdValue))

# 2. Select the relevant columns for transformation.
num_traits <- rtry_select_col (num_traits,
                               ObservationID, AccSpeciesID, AccSpeciesName, 
                               TraitID, TraitName,StdValue, UnitName)

# 3. Extract the values of georeferences and the corresponding ObservationID.
lat <- rtry_select_anc(data, 59)
lon <- rtry_select_anc(data, 60)

# 4. Merge the relevant data frames based on the ObservationID 
# using rtry_join_left().
num_traits_georef <- rtry_join_left (num_traits, lat, baseOn = ObservationID)
num_traits_georef <- rtry_join_left(num_traits_georef, lon, baseOn = ObservationID)


# 5. Perform wide table transformation of TraitID, TraitName and UnitName 
# based on ObservationID, AccSpeciesID and AccSpeciesName with cell 
# values from StdValue. If several records with StdValue were provided for 
# one trait with the same ObservationID, AccSpeciesID and 
# AccSpeciesName, calculate their mean.
num_traits_georef_wider <- rtry_trans_wider(num_traits_georef,
                                            names_from = c(TraitID, TraitName, UnitName),
                                            values_from = c(StdValue),
                                            values_fn = list(StdValue = mean))


# Worked and it is amazing, thank you rtry. writing to csv. 
write.csv(num_traits_georef_wider, file = "data_try_cleaned.csv")



####################
### Need nitrogen and carbon values in % adding column here
####################

