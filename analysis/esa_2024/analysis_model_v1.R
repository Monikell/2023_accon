###############################################################################
## Project: ESA analysis, SLA and cover for modeling
## Author: Monika Kelley
## Date: 2024/02/27
###############################################################################

## load packages
library(tidyverse)
library(dplyr)
library(ggplot2)

## load data, sla, comp, weight data
sla.comp <- read.csv("data/00_data_sla.comp.csv")


##############################################################################
## community weighted means (internet)
##############################################################################
## resource: https://search.r-project.org/CRAN/refmans/BAT/html/cwm.html
## resource: https://rstudio-pubs-static.s3.amazonaws.com/502799_bf78c38d5f0a49a39207a0b0039fd5f5.html 
##############################################################################

summarize_sla.comp_cwm <- sla.comp %>%
  group_by(site, plot, trt, duration) %>%
  summarise(
    sla_cwm = weighted.mean(sla_cm2.g.1, percent_cover)
  )

lm(sla_cwm ~ trt * site, data = summarize_sla.comp_cwm)


##############################################################################
## community weighted means (Nick)
##############################################################################
## Goal: CWM SLA dif b/w trt

# calculate individual weighted mean (iwm)
sla.comp$iwm <- sla.comp$sla_cm2.g.1 * (sla.comp$percent_cover/100)

# calculate cwm 
# note: Σ 1/iwm/#ind/plot = cwm (take the sum of all the iwm values)
sla.comp$cwm <- sum(1/sla.comp$iwm)

# model
# note: lm (SLAcwm ~ trt *site, data = data) 




