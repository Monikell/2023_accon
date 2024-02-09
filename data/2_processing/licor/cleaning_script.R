###############################################################################
### Licor cleaning
###############################################################################

############
## Libraries
############
install.packages("devtools")
devtools::install_github("poales/readLicorData")
library(tidyverse)
library(readLicorData)


##############################################################

## Clean Licor data

licor1 <- licorData(location = "raw/monika_albert_2023/2023-05-09_test_albert.xlsx")
write.csv(licor1, "cleaned/01_2023-05-09_test_albert.csv", row.names = FALSE)

licor2 <- licorData(location = "raw/monika_gibson_2023/2023-05-09_test_gibson.xlsx")
write.csv(licor2, "cleaned/02_2023-05-09_test_gibson.csv", row.names = FALSE)

licor3 <- licorData(location = "raw/monika_gibson_2023/2023-06-07-0814_logdata.xlsx")
write.csv(licor3, "cleaned/03_2023-06-07-0814_logdata.csv", row.names = FALSE)

licor4 <- licorData(location = "raw/monika_gibson_2023/monika_gibson_2023/2023-05-09_test_gibson.xlsx")
write.csv(licor4, "cleaned/04_2023-05-09_test_gibson.csv", row.names = FALSE)

licor5 <- licorData(location = "raw/monika_gibson_2023/monika_gibson_2023/2023-06-07-0814_logdata.xlsx")
write.csv(licor5, "cleaned/05_2023-06-07-0814_logdata.csv", row.names = FALSE)

licor6 <- licorData(location = "raw/monika_gibson_2023/monika_gibson_2023/06_2023-06-16-0733_logdata_prgl2day.xlsx")
write.csv(licor6, "cleaned/06_2023-06-16-0733_logdata_prgl2day.csv", row.names = FALSE)

licor7 <- licorData(location = "raw/monika_yadi_2023/07_2023-06-07-1352_logdata.xlsx")
write.csv(licor7, "cleaned/07_2023-06-07-1352_logdata.csv", row.names = FALSE)

licor8 <- licorData(location = "raw/monika_yadi_2023/08_2023-05-24-1009_logdata.xlsx")
write.csv(licor8, "cleaned/08_2023-05-24-1009_logdata.csv", row.names = FALSE)

licor9 <- licorData(location = "raw/monika_yadi_2023/09_2023-05-25-0935_logdata.xlsx")
write.csv(licor9, "cleaned/09_2023-05-25-0935_logdata.csv", row.names = FALSE)

licor10 <- licorData(location = "raw/monika_yadi_2023/10_2023-05-26-0922_logdata.xlsx")
write.csv(licor10, "cleaned/10_2023-05-26-0922_logdata.csv", row.names = FALSE)

licor11 <- licorData(location = "raw/monika_yadi_2023/11_2023-05-26-1349_logdata.xlsx")
write.csv(licor11, "cleaned/11_2023-05-26-1349_logdata.csv", row.names = FALSE)

licor12 <- licorData(location = "raw/monika_yadi_2023/12_2023-05-27-0821_logdata.xlsx")
write.csv(licor12, "cleaned/12_2023-05-27-0821_logdata.csv", row.names = FALSE)

licor13 <- licorData(location = "raw/monika_yadi_2023/13_2023-06-06-0850_logdata.xlsx")
write.csv(licor12, "cleaned/13_2023-06-06-0850_logdata.csv", row.names = FALSE)










############################################################## Evan Below!

## Clean Licor data
licor1 <- licorData(location = "/Users/eaperkowski/Desktop/eve_test/2023-12-13-soybean-test")
write.csv(licor1, "/Users/eaperkowski/Desktop/eve_test/2023-12-13-soybean-test.csv", 
          row.names = FALSE)

licor2 <- licorData(location = "/Users/eaperkowski/Desktop/eve_test/2023-12-13-soybean-test_b")
write.csv(licor2, "/Users/eaperkowski/Desktop/eve_test/2023-12-13-soybean-test_b.csv")


## Merge licor1 and licor2

file.list <- list.files(path = "/Users/eaperkowski/Desktop/eve_test/",
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



