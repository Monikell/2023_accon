###############################################################################
## purpose: cleaning and fitting licor curves
## author: Monika Kelley
## date: 2024/04/05
###############################################################################

# 2024/04/05: skeletons of code from Evan, need to review/ updated


###############################################################################
## packages and libraries
###############################################################################

install.packages("devtools")
devtools::install_github("poales/readLicorData")
library(tidyverse)
library(readLicorData)

install.packages("rio")
library(rio)



###############################################################################
## loading and prepping data
###############################################################################

## need to mass convert the xsl. files to csv. , look into "rio" package


## load data from file
file.list <- list.files(path = "cleaned_licor_data/",
                        pattern = "*.csv")

## Merge licor1 and licor2

file.list <- list.files(path = "cleaned_licor_data/",
                        pattern = "*.csv")
file.list <- setNames(file.list, stringr::str_extract(basename(file.list), 
                                                      '.*(?=\\.csv)'))


merged_files <- lapply(file_path, read.csv) %>%
  reshape::merge_all()

write.csv(merged_files, "merged_file.csv", row.names = FALSE)

library(tidyverse)
library(plantecophys)

aci.df <- read.csv("merged_file.csv")

## Curve fit test 1 (low light)
curve1 <- aci.df %>% 
  filter(id == "brit.lc.3.11") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(curve1)
summary(curve1)
coef(curve1)

aci.coefs <- data.frame(id = "brit.lc.3.11", doy = "360", t(coef(curve1)))

## Curve fit test 2 (high light)
curve2 <- aci.df %>% 
  filter(id == "brit.hc.4.4") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(curve2)
summary(curve2)
coef(curve2)

aci.coefs[2,] <- c(id = "brit.hc.4.4", baseline = "n", t(coef(curve2)))
aci.coefs[3,] <- c(id = "brit.lc.3.11", day = "7", t(coef(curve1)))