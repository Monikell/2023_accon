

## packages
library(stringr)
library(tidyverse)
library(lubridate)
library(dplyr)
library(stringi)
library(stringr)


# loading data
data_multi <- read.csv(file = "0_data_multi.csv")
soft_lubb <- read.csv(file = "field_soft_lubb.csv")
soft_sevi <- read.csv(file = "field_soft_sevi.csv")
soft_arch <- read.csv(file = "field_soft_arch.csv")
soft_temple <- read.csv(file = "field_soft_temple.csv")

# merging sites
sites_t.a <- rbind(soft_temple,soft_arch)
sites_l.s <- rbind(soft_lubb, soft_sevi)
data_soft_all <- rbind(sites_l.s, sites_t.a)

# lowering code names 
data_multi[[11]] <- tolower(data_multi[[11]])
data_multi[[12]] <- tolower(data_multi[[12]])
data_multi[[12]] <- tolower(data_multi[[13]])
data_soft_all[[6]] <- tolower(data_soft_all[[6]])

# modifying data_multi colnames to match field_soft
colnames(data_multi)[12] <- "trt"
colnames(data_multi)[11] <- "code"
colnames(data_multi)[9] <- "rep"

# sites
unique(data_multi$site)
unique(data_soft_all$site)

# temple not spelled the same repalcing tmpl -> temple
data_soft_all$site <- gsub("tmpl", "temple", data_soft_all$site)

# colname test
colnames(data_soft_all)
colnames(data_multi)

## try to merge with multi
test <- full_join(data_multi, data_soft_all, 
                  by = c("site", "plot", "trt", "code", "rep"))
getwd()

## manual clean
write.csv(data_multi, file = "data_multi_edits.csv")
write.csv(data_soft_all, file = "data_soft_all_edits.csv")

## checking things
unique(data_multi$plot)
unique(data_multi$rep)
unique(data_soft_all$plot)
unique(data_soft_all$rep)


# type conversion
unique(data_multi$plot)
unique(soft_all$plot)
data_multi$plot[is.na(data_multi$plot)] <- 0
soft_all$plot[is.na(soft_all$plot)] <- 0

# removing na values, by changing things to 0
data_multi$plot <- as.integer(data_multi$plot)
data_multi$rep <- as.integer(data_multi$rep)



########### Attempt to merge #2 after excel edits
## issues with mutlipe matches, add time maybe. 
data_multi2 <- read.csv("r_csv/data_multi_edits.csv")
data_soft2 <- read.csv("r_csv/data_soft_all_edits.csv")

## too many things are matching!! Try adding time range
data_merged <- full_join(data_multi2, data_soft2, 
                  by = c("site", "plot", "trt", "code", "rep", "time"))


## trying to merge by data.table 
install.packages("data.table")
library(data.table)



test <- data_soft2[data_multi2, 
                   .(time_10minus = i.time, end_time, value1, value2), 
                   on = .(start_time <= time, end_time >= time)]


