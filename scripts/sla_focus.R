# focusing on SLA


### load libraries -------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(lme4)
library(ggplot2)
library(car) # "Anova"
library(emmeans)
library(multcompView) # "cld" for emmeans
library(multcomp) # "cld" for emmeans 
library(corrr)
library(vegan)

## PCA
#install.packages("corrr")
library('corrr')
#install.packages("ggcorrplot")
library(ggcorrplot)
#install.packages("FactoMineR")
library(FactoMineR)
# install.packages("factoextra")
library(factoextra)

### loading data
## data created on 08/03, stat. has ranking of acquisition strat. based on SLA 
## values
data_spcomp_raw <- read.csv("../data/03_rproducts/spcomp_traits_sla.strat_rank.csv")


###############################################################################
## Calculating diversity indexes
###############################################################################

