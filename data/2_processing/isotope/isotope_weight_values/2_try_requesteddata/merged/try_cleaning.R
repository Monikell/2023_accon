# checking directory 
getwd()

# reading in TRY database data
my_data <- read.delim("29103.txt")

# checking out col. names
colnames(my_data)


# isolating rows with N/C and lat/long info
# 14 = carbon, 15 = Nitrogen, 59 = lat, 60 = long
data_isolate <- subset(my_data, DataID %in% 
                         c (14, 15, 59, 60)) 


# checking info in the OriglName col.
unique(data_isolate$DataName)


## Trying nicks suggestion of subsetting 
# Data id notes: # 14 = carbon, 15 = Nitrogen, 59 = lat, 60 = long

# Making individual dataframes based on teh values I want. 
d_carbon <- subset(data_isolate, DataID %in% c (14))
d_nitrogen <- subset(data_isolate, DataID %in% c (15)) 
d_lat <- subset(data_isolate, DataID %in% c (59)) 
d_long <- subset(data_isolate, DataID %in% c (60)) 


# Chaning col. names, of "StandardValue"
colnames(d_carbon)[21] <- "%carbon"
colnames(d_nitrogen)[21] <- "%nitrogen"
colnames(d_lat)[21] <- "lat"
colnames(d_long)[21] <- "long"


# only columns I want, DataID (3), and StandardRate (21)
only_carbon <- d_carbon[,c(3,21)]
only_nitrogen <- d_nitrogen[,c(3,21)]
only_lat <- d_lat[,c(3,21)]
only_long <- d_long[,c(3,21)]


# Okay! Merging with original d_carbon data set, and
# the "only" datasets. That way have all the original info + 
# the bits and pieces we want. 


# merging with whole dataset = bad, puts values on everything. 
test <- merge(data_isolate, only_carbon)
head(test)
View(test)
? "merge"

# trying to just merge the only data
# all things get a value and I don't think that is right. 
# will write the data down into a file and go from there. 
test <- merge(only_carbon, only_nitrogen)
View(test)

# writing data to look at it. 
# seems okay, it looks like all carbon values have a nitrogen value
# make sense, getting carbon = getting nitrogen. try on lat and long now.
write.csv(data_isolate, file = "data_isolte.csv")
write.csv(test, file = "test.csv")


# merging only data now
# already merged N and C previous, just renaming them here. 
isotopes <- test
# mering long and lats
long_lat <- merge(only_lat, only_long)

# Mering long_lat and isotope
# Not working very well
nc_longlat <- merge(isotopes, long_lat)
View(nc_longlat)


# Loading packages
install.packages("dplyr")
library("dplyr")
install.packages("tidyverse")
library("tidyverse")

? "full_join"


# attempting to repeated rows 75 - 78, joining isotopes and long_lat
# worked but created a massive data set. Looks like I was sorting by 
# the wrong value. Should sort by "ObservationID" and not by "DataID". DataID 
# includes everything in that data set, which is what we don't want. 
# rework the code. 
View(isotopes)
View(long_lat)

merged_data <- isotopes %>% full_join(long_lat)

head(merged_data)





##########################################################################
## Restarting, using "ObservationID" column
##########################################################################
## Note: "ObservationID" column, seems to be unique per species, data set, and
# unique values collected (such as Nitrogen, and Carbon). 


# Loading packages
install.packages("dplyr")
install.packages("tidyverse")
library("dplyr")
library("tidyverse")

# Loading the TRY database informaiton
my_data <- read.delim("29103.txt")

# isolating rows with N/C and lat/long info
# "DataID" values: 14 = carbon, 15 = Nitrogen, 59 = lat, 60 = long
data_isolate <- subset(my_data, DataID %in% c (14, 15, 59, 60)) 

# double checking that subsetted correct data/ rows.
unique(data_isolate$DataName)


# Creating new vectors for each element fo interest to later be joined. 
d_carbon <- subset(data_isolate, DataID %in% c (14))
d_nitrogen <- subset(data_isolate, DataID %in% c (15)) 
d_lat <- subset(data_isolate, DataID %in% c (59)) 
d_long <- subset(data_isolate, DataID %in% c (60)) 


# Updating the names of the "StandardValue" column to make it easier to find
# these values when they join back with the whole dataset. 
colnames(d_carbon)[21] <- "%carbon"
colnames(d_nitrogen)[21] <- "%nitrogen"
colnames(d_lat)[21] <- "lat"
colnames(d_long)[21] <- "long"

# checking the column names to make sure things look good. 
View(d_carbon)


# Isolating columns we want "DatasetID" (3) "ObservationID" (8), 
# and "StandardRate" (21)
# previously used "DatasetID" (3) only = bad.
only_carbon <- d_carbon[,c(3,8,21)]
only_nitrogen <- d_nitrogen[,c(3,8,21)]
only_lat <- d_lat[,c(3,8,21)]
only_long <- d_long[,c(3,8,21)]

# checking data
View(only_carbon)

# Merging isotope data, and lat and long data
# make sure to add "all = TRUE" to add NA values where needed
isotopes <- merge(only_nitrogen, only_carbon, all = TRUE)
lat_long <- merge(only_lat, only_long, all = TRUE)
iso_locations <- merge(isotopes, lat_long, all = TRUE)

View(iso_locations)


# Merging 
merged_data <- merge(data_isolate, iso_locations, all = TRUE)

View(merged_data)

# It worked!! Isotope, lat. and long data all on one row. 
# writing to a csv. Only issue now have duplicate rows. 

write.csv(merged_data, "merged_data.csv")


## Took the above csv. and made edits in excel for ease. 
# removed unwanted rows, removed duplicates that resulted because 
# of this process. And locaitons with NA values. 
