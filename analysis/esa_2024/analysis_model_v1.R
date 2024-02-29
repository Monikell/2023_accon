###############################################################################
## Project: ESA analysis, SLA and cover for modeling
## Author: Monika Kelley
## Date: 2024/02/27
###############################################################################


## load packages
library(tidyverse)
library(car)
library(emmeans)
library(lme4)
library(nlme)

## load data, sla, comp, weight data
sla.comp <- read.csv("data/00_data_sla.comp.csv")

## Todo: need to add block as a random. 

### mod dataset 

sla.comp$block <- NA

# adding blocks
sla.comp$block[sla.comp$site == "arch" & sla.comp$plot < 11] <- "block_1"
sla.comp$block[sla.comp$site == "arch" & sla.comp$plot >= 11& sla.comp$plot < 21] <- "block_2"
sla.comp$block[sla.comp$site == "arch" & sla.comp$plot >= 21] <- "block_3"

sla.comp$block[sla.comp$site == "lubb" & sla.comp$plot < 15] <- "block_1"
sla.comp$block[sla.comp$site == "lubb" & sla.comp$plot >= 15& sla.comp$plot < 29] <- "block_2"
sla.comp$block[sla.comp$site == "lubb" & sla.comp$plot >= 29] <- "block_3"

sla.comp$block[sla.comp$site == "sevi"] <- "block_1"

sla.comp$block[sla.comp$site == "temple" & sla.comp$plot < 9] <- "block_1"
sla.comp$block[sla.comp$site == "temple" & sla.comp$plot >= 11& sla.comp$plot < 21] <- "block_2"
sla.comp$block[sla.comp$site == "temple" & sla.comp$plot >= 21] <- "block_3"


##############################################################################
## community weighted means (internet)
##############################################################################
## resource: https://search.r-project.org/CRAN/refmans/BAT/html/cwm.html
## resource: https://rstudio-pubs-static.s3.amazonaws.com/502799_bf78c38d5f0a49a39207a0b0039fd5f5.html 
##############################################################################

## Model (Nick and Evan)

# calculating cwm (Evan, ditch duration might be weighing down things)
summarize_sla.comp_cwm <- sla.comp %>%
  group_by(site, plot, trt, block) %>%
  summarise(
    sla_cwm = weighted.mean(sla_cm2.g.1, percent_cover)
  )

# Model
# model <- lm(sla_cwm ~ trt * site, data = summarize_sla.comp_cwm) # old model

colnames(summarize_sla.comp_cwm)

sla_lmer <- lme(sla_cwm ~ trt + site, random = ~(1|block), data = summarize_sla.comp_cwm) #nested rdnm 



# check model fit Evan
# Check model assumptions
plot(model)
qqnorm(residuals(model))
qqline(residuals(model))
densityPlot(residuals(model))
shapiro.test(residuals(model)) #test normality, p value greater than .05
outlierTest(model)

# Evan: Model output
summary(model)
Anova(model) # car package tpye to chisquare walds square test for mixed models. Can glean p values.

# post-hoc, emtrends = continous variable, (estimated marginal trend.)
# Evan: future replace site for climate predictors. Can extract MAT MAP (PRISM)
# Would change to linear mixed effect models.

# Individual effects
# tells us how means are arranged, but if we want to do the tukie test. 

#Evan: FL, major diff in VPD (not so much temp.). VPD greater effect on
# VPD. 
emmeans(model, pairwise~ site)






# ##############################################################################
# ## community weighted means (Nick)
# ##############################################################################
# ## Goal: CWM SLA dif b/w trt
# 
# # calculate individual weighted mean (iwm)
# sla.comp$iwm <- sla.comp$sla_cm2.g.1 * (sla.comp$percent_cover/100)
# 
# # calculate cwm 
# # note: Σ 1/iwm/#ind/plot = cwm (take the sum of all the iwm values)
# sla.comp$cwm <- sum(1/sla.comp$iwm)
# 
# # model
# # note: lm (SLAcwm ~ trt *site, data = data) 




