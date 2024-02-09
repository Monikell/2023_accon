getwd()

data <- read.csv(file = "merged_data.csv")

install.packages("tidyverse")
library("tidyverse")

# removing duplicates values from the "ObservationID" column
# left over dirty data from the cleaning file. With the cleaning file 
# got the data into one row, but had duplicates of the extra rows. 
# now removed. 
data1 <- data[!duplicated(data$ObservationID), ]
colnames(data1)


# keeping the columns we want
# Dataset ID, Observation and AccSp. ID and name, lat long and N and C values
data_clean <- data1[c(2,3,8:10,31:34)]
colnames(data_clean)

# writing CSV file of the clean data 
write.csv(data_clean, file = "data_clean.csv")


## Tried to upload into google maps and it gave an error saying that it 
# couldn't do that so just looking at the data. 
# looking for unique values in latitude
unique_values_lat <- unique(data_clean$lat)
print(unique_values_lat)

# unique values in longitude
unique_values_long <- unique(data_clean$long)
print(unique_values_long)



## Trying to add countries and states based on lat and long 
install.packages("sp")
install.packages("rgdal")



