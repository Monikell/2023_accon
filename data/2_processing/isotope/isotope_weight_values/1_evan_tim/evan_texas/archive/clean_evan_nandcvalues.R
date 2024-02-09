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

# mean c values
c_means <- ddply(evan_data, .(site, NCRS.code), 
                      transform, mean.c.value = mean(c.leaf))

CN_means <- ddply(c_means, .(site, NCRS.code),
                  transform, mean.n.value = mean(n.leaf))

# non-drop list
drop <- c("site", "NCRS.code", "duration", "n.leaf", "c.leaf",
          "mean.c.value", "mean.n.value")


# keeping only the columns I want based on the names
test <- CN_means[,(names(CN_means) %in% drop)]

# check to see i've got it
colnames(test)

# writing this into a CV 
write.csv(test, "cn_mean_values_evan.csv")
