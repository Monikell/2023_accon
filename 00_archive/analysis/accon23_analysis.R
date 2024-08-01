################################################################################
## Author: Monika Kelley
## Date: 2024/05/16 to _______
## Purpose: Summer 2023, Acquisitive vs. Conservative Grasslands Communities
## Project code/ name: ACCON (AC-quisitive CON-servative)
################################################################################

### ctrl f notes 
# LEMON = last stopped
# HELP = issue/ not sure


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
library("FactoMineR")





### load data ------------------------------------------------------------------
## data
# soft data = field and lab data collected
raw_data_soft <- read.csv("../data/02_cleaned/field/accon_field.lab_2023.csv")
raw_data_species.comp <- read.csv("../data/02_cleaned/field/accon_species.comp_2023.csv")

## metadata
metadata_plots <- read.csv("../data/00_meta/plot-descriptions-02-August-2019.csv")
metadata_species.list <- read.csv("../data/00_meta/species_list.csv")
metadata_rank <- read.csv("../data/00_meta/ranking_taxon.csv")


################################################################################
## soft traits - calculating soft plant trait information
################################################################################

### cleaning data --------------------------------------------------------------

## renaming data set to be modified
data_soft_calc <- raw_data_soft

## lowering column names and the treatment column Values
names(data_soft_calc) <- tolower(colnames(data_soft_calc))

## adding percent cover to the soft dataset
# creating shorter version of plot metadata
metadata_plot_treatment <- 
  metadata_plots[,c("site_code","Treatment","plot","block")]

# easier to merge col. names
names(metadata_plot_treatment) <- tolower(colnames(metadata_plot_treatment))

## merging plot information with species comp information
data_spcomp.avg.cover <- 
  left_join(raw_data_species.comp, metadata_plot_treatment)


## getting sp.comp cover averages
# MRK: sp.comp average between plots
 df_grouped_spcomp_avg.cover <- data_spcomp.avg.cover %>%
   group_by(site_code, plot, taxon_code) %>%
   summarise(average.cover = mean(percent_cover, na.rm = TRUE))


 ### calculating soft traits ---------------------------------------------------
 
 # traits: leaf thickness, SLA, carbon concentration, nitrogen concentration, 
 # C:N ratio, SSD, plant height, δ13C, δ15N, LDMC 
 # Coming: phosphorous and potassium
 
### specific leaf area (SLA) ---------------------------------------------------
## calculating SLA: area/oven-dry mass
data_soft_calc$sla <- 
  data_soft_calc$leaf_area / data_soft_calc$leaf_dry_g

 
### leaf dry-matter content (LDMC) ---------------------------------------------
## calculating LDMC: oven-dry mass / water-saturated fresh mass
data_soft_calc$ldmc <- 
  data_soft_calc$leaf_dry_g / data_soft_calc$leaf_wet_g


### stem-specific density (SSD) ------------------------------------------------
## calculating SSD: oven-dry mass / volume of fresh stem

## notes ----------
# 2 methods for calculating stem volume: 1) dimensional, 2) water-displacement
# 01) dimensional method: V = (0.5D)^2 X pi X L
# 01) notes) D = diameter (ddh average), L = length of section of stem
# 02) water-displacement: physically calculated in the lab


## 01) dimensional method calculations ----------

## step 1: addressing hollow stems ----------

## step 1.a: renaming columns for clarity
# imagej specifically called out as no "set_scale" was used as we were only 
# looking for a percentage of hollow out of total stem cross-sectional area
data_soft_calc <- data_soft_calc %>%
  rename(
  stem_hollow.or.solid = stem_hollow, # clarifying if stem solid or hollow
  stem_imagej_area.whole = stem_area_whole,
  stem_imagej_area.hollow = stem_area_hollow,
  ddh_average_mm = ddh_average
  )


## step 1.b: calculating average % of hollow stem for each species by site
df_stem_hollow.average <- data_soft_calc %>%
  group_by(site, taxon_code) %>%
  summarise(stem_hollow.percent.average = mean(stem_area_percent_hollow, 
                                               na.rm = TRUE))

## step 1.c: merge the calculated averages back into the original data set
data_soft_calc <- data_soft_calc %>%
  left_join(df_stem_hollow.average, by = c("site", "taxon_code"))

## step 1.d: estimating hollow diameter
# calculate the solid and hollow area of the stem using the ddh average
data_soft_calc$stem_total.area <- 
  pi * (data_soft_calc$ddh_average_mm / 2)^2

# calculate the hollow area based on the calculated percentage
data_soft_calc$stem_hollow.area <- 
  (data_soft_calc$stem_hollow.percent.average / 100) * 
  data_soft_calc$stem_total.area

# calculate the estimated diameter of the hollow area
# purpose: used to check math against field ddh values (reasonable numbers)
data_soft_calc$stem_hollow.diameter_estimate <- 
  2 * sqrt(data_soft_calc$stem_hollow.area / pi)


## step 2: finalizing stem area ----------

# step 2.a: replacing NA values with 0 for the stem_hollow.area
# purpose: highlights solid vs. hollow stems + good for math in step 2.b
data_soft_calc <- data_soft_calc %>% 
  mutate(stem_hollow.area = replace_na(stem_hollow.area, 0))


## step 2.b: calculating solid cross-sectional area: whole.area - hollow.area
## last step for dealing with hollow stems! can move on to volume calculations 
data_soft_calc$stem_area.solid <- 
  data_soft_calc$stem_total.area - data_soft_calc$stem_hollow.area


## step 3: calculating volume  ----------
# formula V = (stem_area.solid) X length of stem

# converting the area from mm^2 to cm^2
data_soft_calc$stem_area.solid_cm2 <- data_soft_calc$stem_area.solid / 100

# calculating volume: V = (0.5D)^2 x pi x L
data_soft_calc$ stem_volume.dimensional <- 
  data_soft_calc$stem_area.solid_cm2 * data_soft_calc$stem_length_cm


## step 4: calculating SSD
data_soft_calc$ssd_dimensional <- 
  data_soft_calc$stem_dry_g / data_soft_calc$stem_volume.dimensional


## 02) water-displacement method calculations ----------

## step 1: updating column names for clarity
data_soft_calc$stem_volume.displacement <- data_soft_calc$stem_water.displaced.g


## step 2: calculating ssd 
# volume physically calculated in lab, no math required other than ssd calc. 
data_soft_calc$ssd_displacement <- 
  data_soft_calc$stem_dry_g / data_soft_calc$stem_volume.displacement


### carbon and nitrogen concentration ------------------------------------------
## total of sample (µg)/ sample weight (mg)

# carbon
data_soft_calc$c_concentration <- 
  data_soft_calc$c_total / data_soft_calc$cn_sample.weight

# nitrogen
data_soft_calc$n_concentration <- 
  data_soft_calc$n_total / data_soft_calc$cn_sample.weight


### carbon nitrogen ration (c:n) -----------------------------------------------
## carbon concentration / nitrogen concentration

data_soft_calc$cn_ratio <- 
  data_soft_calc$c_concentration / data_soft_calc$n_concentration


### soft plant data calculated creating new data frame  ------------------------
data_soft_full <- data_soft_calc

## looking over data for qc/qa
# write.csv(data_soft_full, "../data/03_rproducts/data_soft_full.csv")
# MRK: looks good, except sample "sevi_39_boer_1_491", the CN isotope values
# seemed to have an issue from the lab, removing data error. 

# removing CN isotope issue sample
data_soft_full <- data_soft_full[data_soft_full$id_full != "sevi_39_boer_1_491", ]


################################################################################
## soft traits - models
################################################################################




################################################################################
## community weighted means (CWM) - models
################################################################################

### merging dfs soft and species comp together for CWM -------------------------

## verifying photosynthesis pathways based on delta c 13 values
# Bender 1971; Smith & Epstein 1971: C3 = 13C (-24‰ to -34‰), c4 = (-6‰ to -19‰)
data_soft_full <- data_soft_full %>%
  mutate(photo_pathway_delta.c.13 = case_when(
    c_delta.13 >= -35 & c_delta.13 <= -23 ~ "c3",
    c_delta.13 >= -19 & c_delta.13 <= -6 ~ "c4"
  ))
# MRK: looks good, only issue is one "boer", that has a weird C13 value, lab
# noted that the sample was very low, and that likely explains the value.


## mergining soft and species comp. average cover
data_soft.spcomp <- left_join(data_soft_full, df_grouped_spcomp_avg.cover, 
                   by = c("site_code", "plot", "taxon_code"))

## mergining species list information into dataset 
data_soft.spcomp <- left_join(data_soft.spcomp, metadata_species.list, 
                         by = c("site_code", "taxon_code"))

## removing site.y, and changing site.x to just site, for simplicity
# remove column site.y
data_soft.spcomp$site.y <- NULL

# update name 
names(data_soft.spcomp)[names(data_soft.spcomp) == "site.x"] <- "site"


## sp.comp % cover update NA values from "NA" to "0", top 5 species data
# collected from all plots even if not observed in 1x1s
data_soft.spcomp$average.cover[is.na(data_soft.spcomp$average.cover)] <- 0


## qa/qc real quick
# write.csv(data_soft.spcomp, "../data/03_rproducts/data_soft.spcomp.csv")
# MRK: looks good/ okay! 


### CWM: SLA -------------------------------------------------------------------

## looking for NA values
any(is.na(data_soft.spcomp$sla))

## looking over raw data
hist(data_soft.spcomp$sla)
hist(log(data_soft.spcomp$sla))

qqnorm(data_soft.spcomp$sla) ; qqline(data_soft.spcomp$sla, col = "red")
qqnorm(log(data_soft.spcomp$sla)) ; qqline(log(data_soft.spcomp$sla))


## cwm calculation
cwm_sla_summarize <- data_soft.spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_sla = log(weighted.mean(sla, average.cover)))


## model (nested rdmn) _________________________________________________________ HELP, 'is singluar'
# overcomplicated random effects, potentailyl remove block?
lmer_sla <- 
  lmer(cwm_sla ~ treatment * site + (1 | block:site), 
       data = cwm_sla_summarize)

# # Simplified model with random effect for site only
# lmer_sla_simplified <- lmer(cwm_sla ~ treatment * site + (1 | site), 
#                             data = cwm_sla_summarize)


## plot ________________________________________________________________________ HELP, cone(ish)
plot(lmer_sla, which = 2)

## anova
Anova(lmer_sla) ## site ***, treatment:site * 

## emmeans
emmeans(lmer_sla, ~site)
emmeans(lmer_sla, ~treatment)
emmeans(lmer_sla, ~treatment*site)

# looking deeper, p-value (so close, p = 0.0659)
pairs(emmeans(lmer_sla, ~treatment, at = list(site = 'temple'))) # p 0.1155
pairs(emmeans(lmer_sla, ~treatment, at = list(site = 'arch'))) #p 0.1086
pairs(emmeans(lmer_sla, ~treatment, at = list(site = 'lubb')))
pairs(emmeans(lmer_sla, ~treatment, at = list(site = 'sevi')))


# visualize and interpret pairwise comparison (trt not sig. from each other)
cld(emmeans(lmer_sla, ~treatment, at = list (site = 'temple')))
cld(emmeans(lmer_sla, ~treatment, at = list (site = 'arch')))
cld(emmeans(lmer_sla, ~treatment, at = list (site = 'lubb')))
cld(emmeans(lmer_sla, ~treatment, at = list (site = 'sevi')))

## visualizing
# raw data
ggplot(data_soft.spcomp, aes (site, sla, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1)

# cwm data
# MRK: arch and temple (wet sites) the SLA is significant
ggplot(cwm_sla_summarize, aes (site, cwm_sla, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1)


## checking over data again to make sure it is okay
# q-q plot
residuals_lmer_sla <- residuals(lmer_sla)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_sla)
qqline(residuals_lmer_sla)


### CWM: Leaf thickness --------------------------------------------------------

## to log or not log, looking at data (log)
hist(data_soft_full$leaf_thickness)
hist(log(data_soft_full$leaf_thickness))
qqnorm(data_soft_full$leaf_thickness)
qqline(data_soft_full$leaf_thickness)
qqnorm(log(data_soft_full$leaf_thickness))
qqline(log(data_soft_full$leaf_thickness))

## looking for NA values
any(is.na(data_soft.spcomp$leaf_thickness))
unique(data_soft.spcomp$leaf_thickness)

## cwm calculation
cwm_leaf.thickness_summarize <- data_soft.spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_leaf.thickness = 
              log(weighted.mean(leaf_thickness, average.cover)))

## model (nested rdmn)
lmer_leaf.thickness <- 
  lmer(cwm_leaf.thickness ~ treatment * site + (1 | block:site), 
       data = cwm_leaf.thickness_summarize)

## look to see if a more simple model, get rid of 'singular' issue would be 
## resolved and make a better model, it doesn't look like it. 
# ## model (nested rdmn)
# lmer_leaf.thickness_site <- 
#   lmer(cwm_leaf.thickness ~ treatment * site + (1 | site), 
#        data = cwm_leaf.thickness_summarize)
# plot(lmer_leaf.thickness_site, which = 2)
# 
# 
# AIC(lmer_leaf.thickness, lmer_leaf.thickness_site)
# anova(lmer_leaf.thickness, lmer_leaf.thickness_site)


## plot
plot(lmer_leaf.thickness, which = 2) ## ---------------------------------------- HELP: cone-ish

## anova
Anova(lmer_leaf.thickness) # site ***

## emmeans
emmeans(lmer_leaf.thickness, ~site)
emmeans(lmer_leaf.thickness, ~treatment)
emmeans(lmer_leaf.thickness, ~treatment*site)

# looking deeper
pairs(emmeans(lmer_leaf.thickness, ~treatment, at = list(site = 'temple')))
pairs(emmeans(lmer_leaf.thickness, ~treatment, at = list(site = 'lubb')))
pairs(emmeans(lmer_leaf.thickness, ~treatment, at = list(site = "sevi")))

# visualize and interpret pairwise comparison (trt not sig. from each other)
cld(emmeans(lmer_leaf.thickness, ~treatment, at = list (site = 'temple')))

## visualizing
# raw data
ggplot(data_soft.spcomp, aes (site, leaf_thickness, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1)

# cwm data
ggplot(cwm_leaf.thickness_summarize, aes (site, cwm_leaf.thickness, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1)


## q-q plot
residuals_lmer_leaf.thickness <- residuals(lmer_leaf.thickness)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_leaf.thickness)
qqline(residuals_lmer_leaf.thickness)


### CWM: C:N ratio -------------------------------------------------------------

## histo & qq
hist(data_soft_full$cn_ratio)
hist(log(data_soft_full$cn_ratio))
qqnorm(data_soft_full$cn_ratio)
qqline(data_soft_full$cn_ratio)
qqnorm(log(data_soft_full$cn_ratio))
qqline(log(data_soft_full$cn_ratio))

## looking for NA values
any(is.na(data_soft.spcomp$cn_ratio))
unique(data_soft.spcomp$cn_ratio)

## cwm calculation
cwm_cn.ratio_summarize <- data_soft.spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_cn_ratio = weighted.mean(cn_ratio, average.cover,
                                                  na.rm = TRUE))

## model (nested rdmn)
lmer_cn.ratio <- 
  lmer(cwm_cn_ratio ~ treatment * site + (1 | block:site), ## no singular error
       data = cwm_cn.ratio_summarize)

## plot
plot(lmer_cn.ratio, which = 2)

## anova
Anova(lmer_cn.ratio) # site ***, treatment ***

## emmeans
emmeans(lmer_cn.ratio, ~site)
emmeans(lmer_cn.ratio, ~treatment)
emmeans(lmer_cn.ratio, ~treatment*site)

## looking just at the treatment
pairs(emmeans(lmer_cn.ratio, ~treatment))
cld(emmeans(lmer_cn.ratio, ~treatment))

# looking deeper, p-value
pairs(emmeans(lmer_cn.ratio, ~treatment, at = list(site = 'temple')))
pairs(emmeans(lmer_cn.ratio, ~treatment, at = list(site = 'lubb')))
pairs(emmeans(lmer_cn.ratio, ~treatment, at = list(site = 'sevi')))
pairs(emmeans(lmer_cn.ratio, ~treatment, at = list(site = 'arch')))

# visualize and interpret pairwise comparison
cld(emmeans(lmer_cn.ratio, ~treatment, at = list (site = 'temple'))) # 2 groups
cld(emmeans(lmer_cn.ratio, ~treatment, at = list (site = 'lubb')))
cld(emmeans(lmer_cn.ratio, ~treatment, at = list (site = 'sevi')))
cld(emmeans(lmer_cn.ratio, ~treatment, at = list (site = 'arch'))) # 2 groups



## just looking at the data... ugh
# Create a plot for the main effect of treatment
ggplot(cwm_cn.ratio_summarize, aes(x = treatment, y = cwm_cn_ratio, fill = treatment)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Effect of Treatment on C:N Ratio",
       x = "Treatment",
       y = "Community Weighted Mean C:N Ratio")

# Create a plot for the main effect of site
ggplot(cwm_cn.ratio_summarize, aes(x = site, y = cwm_cn_ratio, fill = site)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Effect of Site on C:N Ratio",
       x = "Site",
       y = "Community Weighted Mean C:N Ratio")


## q-q plot
residuals_lmer_cn.ratio <- residuals(lmer_cn.ratio)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_cn.ratio)
qqline(residuals_lmer_cn.ratio)


### CWM: nitrogen concentration (leaf) -----------------------------------------

## leaf nitrogen total
hist(data_soft_full$n_concentration)
hist(log(data_soft_full$n_concentration))
qqnorm(data_soft_full$n_concentration)
qqline(data_soft_full$n_concentration)
qqnorm(log(data_soft_full$n_concentration))
qqline(log(data_soft_full$n_concentration)) # ------------------------------------ HELP, kind of s shaped...

## looking for NA values
any(is.na(data_soft.spcomp$n_concentration))
unique(data_soft.spcomp$n_concentration)

## cwm calculation
cwm_n.concentration_summarize <- data_soft.spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_n.concentration = log(weighted.mean(n_concentration, 
                                                    average.cover, 
                                                    na.rm = TRUE)))

## model (nested rdmn)
lmer_n.concentration <- 
  lmer(cwm_n.concentration ~ treatment * site + 
         (1 | block:site), data = cwm_n.concentration_summarize)

## plot
plot(lmer_n.concentration, which = 2) # ---------------------------------------- HELP, kind of cone shaped

## anova
Anova(lmer_n.concentration) # site ***, treatment **

## emmeans
emmeans(lmer_n.concentration, ~site)
emmeans(lmer_n.concentration, ~treatment)
emmeans(lmer_n.concentration, ~treatment*site)

# looking deeper, p-value
pairs(emmeans(lmer_n.concentration, ~treatment, at = list(site = 'temple')))
pairs(emmeans(lmer_n.concentration, ~treatment, at = list(site = 'lubb')))
pairs(emmeans(lmer_n.concentration, ~treatment, at = list(site = 'sevi')))
pairs(emmeans(lmer_n.concentration, ~treatment, at = list(site = 'arch')))
pairs(emmeans(lmer_n.concentration, ~treatment))


# visualize and interpret pairwise comparison (trt not sig. from each other)
cld(emmeans(lmer_n.concentration, ~treatment, at = list (site = 'sevi')))
cld(emmeans(lmer_n.concentration, ~treatment))


## q-q plot
residuals_lmer_n.concentration <- residuals(lmer_n.concentration)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_n.concentration)
qqline(residuals_lmer_n.concentration)


### CWM: SSD  ------------------------------------------------------------------

## histo
hist(data_soft_full$ssd_dimensional)
hist(log(data_soft_full$ssd_dimensional)) # much better

qqnorm(data_soft_full$ssd_dimensional) ; qqline(data_soft_full$ssd_dimensional)
qqnorm(log(data_soft_full$ssd_dimensional))
qqline(log(data_soft_full$ssd_dimensional)) # much better

## looking for NA values
any(is.na(data_soft.spcomp$ssd_dimensional))
unique(data_soft.spcomp$ssd_dimensional)

## cwm calculation
cwm_ssd_summarize <- data_soft.spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_ssd = log(weighted.mean(ssd_dimensional, 
                                    average.cover, na.rm = TRUE)))

## model (nested rdmn)
lmer_ssd <- lmer(cwm_ssd ~ treatment * site + (1 | block:site), 
                 data = cwm_ssd_summarize)

## plot
plot(lmer_ssd, which = 2) #----------------------------------------------------- HELP, cone ish..

## anova
Anova(lmer_ssd) # site.

## emmeans
emmeans(lmer_ssd, ~site)
emmeans(lmer_ssd, ~treatment)
emmeans(lmer_ssd, ~treatment*site)

# looking deeper, p-value
pairs(emmeans(lmer_ssd, ~treatment, at = list(site = 'temple')))
pairs(emmeans(lmer_ssd, ~treatment, at = list(site = 'lubb')))
pairs(emmeans(lmer_ssd, ~treatment, at = list(site = 'sevi')))
pairs(emmeans(lmer_ssd, ~treatment, at = list(site = 'arch')))
pairs(emmeans(lmer_ssd, ~treatment))

# visualize and interpret pairwise comparison (trt not sig. from each other)
cld(emmeans(lmer_ssd, ~treatment, at = list (site = 'sevi')))
cld(emmeans(lmer_ssd, ~treatment))
cld(emmeans(lmer_ssd, ~treatment*site))

## q-q plot
residuals_lmer_ssd <- residuals(lmer_ssd)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_ssd); qqline(residuals_lmer_ssd) #------------------------ HELP, kind of wonky


### CWM: plant height ----------------------------------------------------------

## histo
hist(data_soft_full$plant_height)
hist(log(data_soft_full$plant_height))
qqnorm(data_soft.spcomp$plant_height) ; qqline(data_soft.spcomp$plant_height)
qqnorm(log(data_soft.spcomp$plant_height))
qqline(log(data_soft.spcomp$plant_height))


## looking for NA values
any(is.na(data_soft.spcomp$plant_height))
unique(data_soft.spcomp$plant_height)

## cwm calculation
cwm_plant.height_summarize <- data_soft.spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_plant.height = log(weighted.mean(plant_height, average.cover)))

## model (nested rdmn)
lmer_plant.height <- 
  lmer(cwm_plant.height ~ treatment * site + (1 | block:site), 
       data = cwm_plant.height_summarize)

## plot
plot(lmer_plant.height, which = 2)

## anova
Anova(lmer_plant.height) # site ***

## emmeans
emmeans(lmer_plant.height, ~site) # temple & lubb
emmeans(lmer_plant.height, ~treatment)
emmeans(lmer_plant.height, ~treatment*site)

# looking deeper, p-value
# temple close with a p of 0.0894
pairs(emmeans(lmer_plant.height, ~treatment, at = list(site = 'temple')))
pairs(emmeans(lmer_plant.height, ~treatment, at = list(site = 'lubb')))

# visualize and interpret pairwise comparison (trt not sig. from each other)
cld(emmeans(lmer_plant.height, ~treatment, at = list (site = 'temple')))
cld(emmeans(lmer_plant.height, ~treatment, at = list (site = 'lubb')))

## q-q plot
residuals_lmer_plant.height <- residuals(lmer_plant.height)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_plant.height) ; qqline(residuals_lmer_plant.height)


### CWM: LDMC ------------------------------------------------------------------

## leaf dry-matter content
hist(data_soft_full$ldmc)
hist(log(data_soft_full$ldmc))
qqnorm(data_soft.spcomp$ldmc) ; qqline(data_soft_full$ldmc)
qqnorm(log(data_soft.spcomp$ldmc)) ; qqline(log(data_soft.spcomp$ldmc))

## looking for NA values
any(is.na(data_soft.spcomp$ldmc))
unique(data_soft.spcomp$ldmc)

## log needed
hist(data_soft.spcomp$ldmc)
hist(log(data_soft.spcomp$ldmc))

## cwm calculation
cwm_ldmc_summarize <- data_soft.spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_ldmc = log(weighted.mean(ldmc, average.cover,
                                                 na.rm = TRUE)))

## model (nested rdmn)
lmer_ldmc <- 
  lmer(cwm_ldmc ~ treatment * site + (1 | block:site), 
       data = cwm_ldmc_summarize)

## plot
plot(lmer_ldmc, which = 2)

## anova
Anova(lmer_ldmc) # site ***

## emmeans
emmeans(lmer_ldmc, ~site)
emmeans(lmer_ldmc, ~treatment)
emmeans(lmer_ldmc, ~treatment*site)

# looking deeper, p-value
# temple close with a p of 0.0894
pairs(emmeans(lmer_ldmc, ~treatment, at = list(site = 'temple')))

# visualize and interpret pairwise comparison (trt not sig. from each other)
cld(emmeans(lmer_ldmc, ~treatment, at = list (site = 'temple')))


## q-q plot
residuals_lmer_ldmc <- residuals(lmer_ldmc)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_ldmc) ; qqline(residuals_lmer_ldmc)


## CWM: δ13C -------------------------------------------------------------------

### δ13C ---------- (data all together)
## data is bi-modal, split between c3 vs. c4 #----------------------------------HELP, bimodial 1/2
 
## looking at the data
hist(data_soft.spcomp$c_delta.13)
qqnorm(data_soft.spcomp$c_delta.13) ; qqline(data_soft.spcomp$c_delta.13) 

## looking for NA values
any(is.na(data_soft.spcomp$c_delta.13))
unique(data_soft.spcomp$c_delta.13) # negative values no log

## cwm calculation ----------
cwm_c.delta.13_summarize <- data_soft.spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_c_delta.13 = weighted.mean(c_delta.13, average.cover, 
                                           na.rm = TRUE))

## model (nested rdmn)
lmer_c_delta.13 <- 
  lmer(cwm_c_delta.13 ~ treatment * site + (1 | block:site), 
       data = cwm_c.delta.13_summarize)

## plot
plot(lmer_c_delta.13, which = 2)

## anova
Anova(lmer_c_delta.13) # site *

## emmeans
emmeans(lmer_c_delta.13, ~site)
emmeans(lmer_c_delta.13, ~treatment)
emmeans(lmer_c_delta.13, ~treatment*site)

# looking deeper, p-value
pairs(emmeans(lmer_c_delta.13, ~treatment, at = list(site = 'sevi')))
pairs(emmeans(lmer_c_delta.13, ~treatment, at = list(site = 'temple')))
pairs(emmeans(lmer_c_delta.13, ~treatment, at = list(site = 'arch')))
pairs(emmeans(lmer_c_delta.13, ~treatment, at = list(site = 'lubb')))

# visualize and interpret pairwise comparison (trt not sig. from each other)
cld(emmeans(lmer_c_delta.13, ~treatment, at = list (site = 'temple')))

## q-q plot
residuals_lmer_c_delta.13 <- residuals(lmer_c_delta.13)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_c_delta.13)
qqline(residuals_lmer_c_delta.13)


### δ13C (splitting c3 vs. c4 plants) ------------------------------------------HELP bimodal 2/2
## subsetting data by photo pathway, as the values can vary based on c3 vs. c4
data_c3.only <- subset(data_soft.spcomp, photo_pathway == "c3")
data_c4.only <- subset(data_soft.spcomp, photo_pathway == "c4")

hist(data_c3.only$c_delta.13)
hist(data_c4.only$c_delta.13)

qqnorm(data_c3.only$c_delta.13) ; qqline(data_c3.only$c_delta.13)
qqnorm(data_c4.only$c_delta.13) ; qqline(data_c4.only$c_delta.13)


## cwm calculation - photo. pathway included -----------
cwm_c.delta.13_summarize_photo <- data_soft.spcomp %>%
  group_by(site, plot, treatment, block, photo_pathway) %>% 
  summarise(cwm_c_delta.13 = weighted.mean(c_delta.13, average.cover, 
                                           na.rm = TRUE))

## model (nested rdmn) - photo. pathway included
lmer_c_delta.13_photo <- 
  lmer(cwm_c_delta.13 ~ treatment * site * photo_pathway + (1 | block:site), 
       data = cwm_c.delta.13_summarize_photo)

# plot
plot(lmer_c_delta.13_photo, which = 2)

# anova
Anova(lmer_c_delta.13_photo) # site, photo_pathway *** & site:photo_path

# emmeans
emmeans(lmer_c_delta.13_photo, ~site)




## CWM: δ15N -------------------------------------------------------------------

## histo
hist(data_soft.spcomp$n_delta.15)
qqnorm(data_soft.spcomp$n_delta.15) ; qqline(data_soft.spcomp$n_delta.15)

## looking for NA values
any(is.na(data_soft.spcomp$n_delta.15))
unique(data_soft.spcomp$n_delta.15) # negative values no log

## cwm calculation
cwm_n.delta.15_summarize <- data_soft.spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_n_delta.15 = weighted.mean(n_delta.15, average.cover, 
                                           na.rm = TRUE))

## model (nested rdmn)
lmer_n_delta.15 <- 
  lmer(cwm_n_delta.15 ~ treatment * site + (1 | block:site), 
       data = cwm_n.delta.15_summarize)

## plot
plot(lmer_n_delta.15, which = 2)

## anova
Anova(lmer_n_delta.15) # site ***, treatment:site ***

## emmeans
emmeans(lmer_n_delta.15, ~site)
emmeans(lmer_n_delta.15, ~treatment)
emmeans(lmer_n_delta.15, ~treatment*site)

# looking deeper, p-value
# temple close with a p of 0.0894
pairs(emmeans(lmer_n_delta.15, ~treatment, at = list(site = 'sevi')))
pairs(emmeans(lmer_n_delta.15, ~treatment, at = list(site = 'temple')))
pairs(emmeans(lmer_n_delta.15, ~treatment, at = list(site = 'arch')))
pairs(emmeans(lmer_n_delta.15, ~treatment, at = list(site = 'lubb'))) #p.0.0021

# visualize and interpret pairwise comparison (trt not sig. from each other)
cld(emmeans(lmer_n_delta.15, ~treatment, at = list (site = 'lubb'))) ## YAY!

## q-q plot
residuals_lmer_n_delta.15 <- residuals(lmer_n_delta.15)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_n_delta.15)
qqline(residuals_lmer_n_delta.15)



################################################################################
## soft traits - PCA
# https://www.youtube.com/watch?v=Tjxgd9FLeYc&t=79s
################################################################################ 
# ## not sure if this is right? -------------------------------------------------- PCA 1 video
# ## better option below (i think)
# 
# ## subset data, traits of interest
# pca_soft.traits <- data_soft.spcomp[,c("leaf_thickness", "sla", "n_concentration", 
#                               "cn_ratio", "ssd_dimensional", "plant_height", 
#                               "c_delta.13", "n_delta.15", "ldmc")]
# 
# ## qa/qc
# # write.csv(pca_soft.traits, "../data/03_rproducts/pca_soft.traits.csv")
# 
# # removing NA values
# any(is.na(pca_soft.traits))
# pca_soft_raw_clean <- na.omit(pca_soft.traits)
# any(is.na(pca_soft_raw_clean))
# 
# ## running the analysis
# pca1 <- rda(pca_soft_raw_clean)
# 
# ## look at total inertia (variance) and that explained by each PC
# pca1
# summary(pca1)
# 
# ## determine how many axes to retain for interpretation (BiodiversityR, package)
# # MRK: focus on % > bs%, only PC1 has "1.00000"
# PCAsignificance(pca1)
# 
# ## plot the data, really basic plot
# ordiplot(pca1)
# ordiplot(pca1, type = "t") # t = text
# 
# 
# ### plot the PCA using the package ----------
# ## using package 'ggvegan'
# autoplot(pca1, legend.position = "none") +
#   xlab("PC1 (92%)") +
#   ylab("PC2 (5%") +
#   geom_abline(intercept = 0, slope = 0, linetype = "dashed", size = 0.8) +
#   geom_vline(aes(xintercept=0), linetype = "dashed", size = 0.8) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"))


### PCA RAW with datacamp ---------------------------------------------------------- PCA 2 attempt (raw data)
# https://www.datacamp.com/tutorial/pca-analysis-r
# https://www.youtube.com/watch?v=5vgP05YpKdE 


## subset data, traits of interest
pca_soft_raw <- data_soft.spcomp[,c("leaf_thickness", "sla", "n_concentration", 
                               "cn_ratio", "ssd_dimensional", "plant_height", 
                               "c_delta.13", "n_delta.15", "ldmc")]
 

head(pca_soft_raw)
## saving/ looking over data
# write.csv(pca_soft.traits, "../data/03_rproducts/pca_soft_raw.csv")
 
# cleaning data for PCA, specifically removing NA values
any(is.na(pca_soft_raw))
colSums(is.na(pca_soft_raw))
colSums(is.na(pca_soft_raw))
pca_soft_raw_clean <- na.omit(pca_soft_raw) # pca to be used for raw data!
any(is.na(pca_soft_raw_clean))
colSums(is.na(pca_soft_raw_clean))

## Normalize the data (is this an issue for later? not sure)
pca_soft_raw_clean_normalized <- scale(pca_soft_raw_clean)
head(pca_soft_raw_clean_normalized)

## applying pca
pca_data.raw <- princomp(pca_soft_raw_clean_normalized)
summary(pca_data.raw )
# summary explained, 9 principal components have been generated (Comp.1 to 9)
# matches number of variables put in. Each component explains a percentage of
# the total variance in the data set. In the "Cumulative Proportion" area, the 
# first pca explains 28% of the total variance. 2 pca = 16% (total 44%) and so 
# on. This is great BUT what do they really mean?? Can be answered by the below.
# Which, is exploring how they related to each column using the loadings of 
# each principal component. 

## The Loading Matrix
pca_data.raw $loadings[,1:6] # Comp.1 - Comp.6 (cumulative 87%)

# The Loading Matrix shows that the first principal component has high positive
# values leaf_thickness, sla, n_concentration, n_delta.15 others are negative. 
# the example talks about is in response to like countires (row columns), here 
# not sure.. like just overall prefrence? 

### Visualization of the principal components ----------

## Scree Plot
# visualize the importance of each principal component and can be used to 
# determine the number of principal components to retain. 
fviz_eig(pca_data.raw , addlabels = TRUE)

# install.packages("factoextra")
# library(factoextra)


## Biplot of the Attributes
# with the biplot, it is possible to visualize the similarities and 
# dissimilarities between the samples, and further show the impact of each 
# attributes on each of the principal components
# Graph of the variables 
fviz_pca_var(pca_data.raw , col.var = "black")



fviz_pca_biplot(pca_data.raw, 
                geom.ind = "point", 
                col.ind = "blue",        # Customize color for individuals
                fill.ind = "red",        # Customize fill color for individuals
                col.var = "black",       # Customize color for variables
                repel = TRUE)           # Avoid overlapping text labels


# what you get from the plot: 
# 1) grouped variables = positively correlated with each other
# 2) higher the distance between the variable & the origin = the better 
#    represented that variable is. Example, leaf thickness better rep. than sla
# 3) variables negatively correlated are displayed opp. sides of the origin

## Contribution of each variable
# 3rd visual determine how much each variable is represented in a given 
# component (quality called Cos2) and corresponds to the square cosine. 
# low value = variable is not perfectly represented by that component. 
# high value = good representation of the variable on that component
fviz_cos2(pca_data.raw , choice = "var", axes = 1:2 )
# code above computed the square cosine value for each variable w/ respect to 
# the first two principal components. n_concentration, cn_ratio, plant_height
# n_delta15, leaf thickness contribute the most to PC1 and PC2..

## Biplot combined with cos2
# last 2 visualization approaches: biplot & attributes importance can be combined
# to create a single biplot, where attributes w/ similar cos2 scores will have
# similar colors. 
fviz_pca_var(pca_data.raw , col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)
# high co2 attributes = green (n_con, cn_ratio.. ), mid orange, low black



### PCA CW data with datacamp -------------------------------------------------- PCA 3 


## data, merging together cwm 
cwm_all_traits <- cwm_cn.ratio_summarize %>%
  full_join(cwm_ldmc_summarize) %>%
  full_join(cwm_leaf.thickness_summarize) %>%
  full_join(cwm_n.concentration_summarize) %>%
  full_join(cwm_n.delta.15_summarize) %>%
  full_join(cwm_plant.height_summarize) %>%
  full_join(cwm_sla_summarize) %>%
  full_join(cwm_ssd_summarize) %>%
  full_join(cwm_c.delta.13_summarize)


## subset data, traits of interest
pca_cwm_traits <- cwm_all_traits[,c("cwm_cn_ratio", "cwm_ldmc", 
                                  "cwm_leaf.thickness", "cwm_n.concentration",
                                  "cwm_n_delta.15", "cwm_plant.height", 
                                  "cwm_sla", "cwm_ssd", "cwm_c_delta.13")]

## updating column names for label
pca_cwm_traits <- pca_cwm_traits %>%
  rename(
    CN_ratio = cwm_cn_ratio,
    ldmc = cwm_ldmc,
    leaf_thickness = cwm_leaf.thickness,
    N_concentration = cwm_n.concentration, 
    N.delta_15 = cwm_n_delta.15, 
    plant_height = cwm_plant.height,
    sla = cwm_sla, 
    ssd = cwm_ssd,
    C.delta_13 = cwm_c_delta.13
  )



## checking for NA values
colSums(is.na(pca_cwm_traits))
any(is.na(pca_cwm_traits))

## normalizing the data
pca_cwm.traits_normalized <- scale(pca_cwm_traits)
head(pca_cwm.traits_normalized)

## applying PCA 
pca_cwm <- princomp(pca_cwm.traits_normalized)
summary(pca_cwm) # comp.1 (40%), comp.2 (27%), comp.3 (15%) = 83%

## loading matrix
pca_cwm$loadings[,1:4] # 1-4 = 90%

## Scree Plot (visual of the loading matrix)
fviz_eig(pca_cwm, addlabels = TRUE)

## Biplot of the attributes
biplot_pca.cwm <- fviz_pca_var(pca_cwm, col.var = "black")

## saving as an image
# png("../figures/pca_cwm.png",
#      width = 8, height = 8, units = 'in', res = 1500)
# biplot_pca.cwm
# dev.off()

## Contribute of each variable
fviz_cos2(pca_cwm, choice = "var", axes = 1:2)





## Biplot combined with cos2
fviz_pca_var(pca_cwm, col.var = "cos2",
             gradient.cols = c("lightblue", "blue", "black"), 
             repel = TRUE)





## attempting to plot more info

# Create the PCA biplot with both variables and data points
biplot <- fviz_pca_biplot(pca_cwm, 
                          col.var = "cos2", 
                          col.ind = "cos2",
                          gradient.cols = c("lightblue", "blue", "black"), 
                          repel = TRUE) +
  ggtitle("PCA Biplot: Variables and Individuals") +
  xlab("Principal Component 1") +
  ylab("Principal Component 2")


biplot


biplot <- fviz_pca_biplot(pca_cwm, 
                          geom.ind = "point",  # Use points for individuals
                          col.ind = "cos2",    # Color individuals by their squared cosine (cos2) values
                          col.var = "black",   # Color variables in black
                          gradient.cols = c("lightblue", "blue", "black"), 
                          repel = TRUE) +      # Avoid label overlap
  ggtitle("PCA Biplot: Community Weighted Means") +
  xlab("Principal Component 1 (40.5%)") +
  ylab("Principal Component 2 (27.9%")





### scoring species ------------------------------------------------------------
species_scores <- pca_cwm$



################################################################################
## species composition - calculating diversity indexes
################################################################################

### cleaning up dataset --------------------------------------------------------

## adding in plot meta data
spcomp_plots <- left_join(raw_data_species.comp, metadata_plots,
                          by = c("site_code","plot"))

## adding in acq. con ranking 
spcomp_plots <- left_join(spcomp_plots, metadata_rank, 
                          by = c("site_code", "taxon_code"))


## lowering column names and the treatment column Values
names(spcomp_plots) <- tolower(colnames(spcomp_plots))
spcomp_plots$treatment <- tolower(spcomp_plots$treatment)


## isolating plant only data
spcomp_plants <- subset(spcomp_plots, type_of_cover == "plant")

## just looking over data
# write.csv(spcomp_plants, "../data/03_rproducts/spcomp_plots_rank.csv")


## removing NA values, that don't consider acq. con.
spcomp_plants <- subset(spcomp_plants, !is.na(rank_sla))


## square percent covers for diversity estimate
spcomp_plants$percent.cover.squared <- 
  (spcomp_plants$percent_cover/100) * (spcomp_plants$percent_cover/100)


## grouping by site: plot | treatment
spcomp_grouped.site.plot <- group_by(spcomp_plants, site_code, plot)


# spcomp_grouped.site.treatment <- group_by(spcomp_plants, site_code, treatment)

spcomp_grouped.site.plot.rank <- 
  group_by(spcomp_plants, site_code, plot, rank_sla)

View(spcomp_grouped.site.plot.rank)

### richness -------------------------------------------------------------------

## site, plot
spcomp_richness_site.plot <- 
  summarise(spcomp_grouped.site.plot, 
            richness_plot = n_distinct(taxon_code))


## site, plot, rank
spcomp_grouped.site.plot.rank <- 
  summarise(spcomp_grouped.site.plot.rank, 
            richness_plot = n_distinct(taxon_code))


# ## site, treatment
# spcomp_richness_site.treatment <- 
#   summarise(spcomp_grouped.site.treatment, richness_treatment = 
#               n_distinct(taxon_code))



### diversity ------------------------------------------------------------------

## site, plot
spcomp_diversity_PLOT <- 
  summarise(spcomp_grouped.site.plot, 
            diversity_plot = 1/sum(percent.cover.squared))

## site, plot
spcomp_diversity_RANK  <- 
  summarise(spcomp_grouped.site.plot.rank, 
            diversity_plot = 1/sum(percent.cover.squared)) ##################### help is this even needed 


# ## site, treatment
# spcomp_diversity_site.treatment <- 
#   summarise(spcomp_grouped.site.treatment, 
#             diversity_treatment = 1/sum(percent.cover.squared))


### evenness -------------------------------------------------------------------

## merging diversity and richness for PLOTS ----------
spcomp_evenness_plots <- 
  left_join(spcomp_diversity_site.plot, spcomp_richness_site.plot, 
            by = c("site_code","plot"))


# merging diversity for plots with metadata plot treatment information
spcomp_evenness_plots <- 
  left_join(spcomp_evenness_plots, metadata_plot_treatment)


## merging diversity and richness for TREATMENT ----------
spcomp_evenness_treatment <- 
  left_join(spcomp_diversity_site.treatment, spcomp_richness_site.treatment)


## calculating evenness ----------

## plots
spcomp_evenness_plots$evenness <- 
  spcomp_evenness_plots$diversity_plot /
  spcomp_evenness_plots$richness_plot


## treatment
spcomp_evenness_treatment$evenness <- 
  spcomp_evenness_treatment$diversity_treatment /
  spcomp_evenness_treatment$richness_treatment


################################################################################
## species composition - models
################################################################################

## convert parameters; vector into a factor, (data categorical vs. continuous)
spcomp_evenness_plots$plotfac <- as.factor(spcomp_evenness_plots$plot)
spcomp_evenness_plots$sitefac <- as.factor(spcomp_evenness_plots$site_code)
spcomp_evenness_plots$trtfac <- as.factor(spcomp_evenness_plots$treatment)
spcomp_evenness_plots$block <- as.factor(spcomp_evenness_plots$block)

## updating data name for simplicity
spcomp_data_4lmer <- spcomp_evenness_plots

View(spcomp_data_4lmer)

### model 01: p.comp diversity across sites ------------------------------------ LEMON

# histos
hist(spcomp_data_4lmer$diversity_plot)
hist(spcomp_data_4lmer$richness_plot)
hist(spcomp_data_4lmer$evenness)

hist(log(spcomp_data_4lmer$diversity_plot))
hist(log(spcomp_data_4lmer$richness_plot))
hist(log(spcomp_data_4lmer$evenness))

# qq plots
qqnorm(log(spcomp_data_4lmer$diversity_plot))
qqline(log(spcomp_data_4lmer$diversity_plot))


## model
mod_div.site.trt <-lmer(log(diversity_plot) ~ sitefac * trtfac + 
                          (1 | plotfac) + (1 | block), 
                        data = (log*spcomp_data_4lmer))


## Q-Q plot for the residuals
plot(resid(mod_div.site.trt) ~fitted(mod_div.site.trt))
plot(mod_div.site.trt, which = 2) # plot residual too, but with a line
# MRK: looks like some heteroscedasticity? -------------------------------------- HELP


## anova
# MRK: site most influence
Anova(mod_div.site.trt) ## sitefac *** 

## post-hoc pairwise comparisons of the estimated marginal means (emmeans)
cld(emmeans(mod_div.site.trt, ~sitefac))


### model 02:  sp.comp evenness across sites -----------------------------------

## model 
mod_evenness.site.trt <-lmer(log(evenness) ~ sitefac * trtfac + 
                          (1 | plotfac) + (1 | block), 
                        data = (spcomp_data_4lmer))

## Q-Q plot for the residuals
plot(resid(mod_evenness.site.trt) ~fitted(mod_evenness.site.trt))
plot(mod_evenness.site.trt, which = 2)

## anova
Anova(mod_evenness.site.trt) ## sitefac ***

## emmeans 
cld(emmeans(mod_evenness.site.trt, ~sitefac)) 




################################################################################
## data visualization
################################################################################

### raw data -------------------------------------------------------------------

plot_raw.leaf_thickness <- 
  ggplot(data_soft_full, aes (site, leaf_thickness, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) +
  theme_minimal()


plot_raw.sla <- 
  ggplot(data_soft_full, aes (site, sla, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  theme_minimal() # temple, outliers might be okay..


## attempting with stats..SLA, looks amazing and has error bars!!
plot_raw.sla_errorbars <- ggplot() + 
    # Add error bars to the box plot for the specified subset of data
    stat_boxplot(
      data = data_soft.spcomp,
      aes(x = site, y = sla), 
      size = 0.75, 
      geom = "errorbar", 
      width = 0.2
    ) +
    # Add the box plot itself for the same subset of data
    geom_boxplot(
      data = data_soft.spcomp,
      aes(x = site, y = sla), 
      outlier.shape = NA
    ) + 
    theme_minimal()


## works okay facet by site, and treatment as the x.. works good
## likely looks better for other things/ sites
ggplot() + 
  # Add error bars to the box plot for the specified subset of data
  stat_boxplot(data = data_soft.spcomp,
    aes(x = treatment, y = sla), 
    size = 0.75, 
    geom = "errorbar", 
    width = 0.2 ) +
  # Add the box plot itself for the same subset of data
  geom_boxplot(data = data_soft.spcomp,
    aes(x = treatment, y = sla, fill = treatment), 
    outlier.shape = NA) + 
  # Facet by site (one row per site)
  facet_grid(rows = vars(site), scales = "free_y") +
  # Minimal theme
  theme_minimal()



plot_raw.n_concentraiton <- 
  ggplot(data_soft_full, aes (site, n_concentration, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  theme_minimal()


plot_raw.cn_ratio <- 
  ggplot(data_soft_full, aes (site, cn_ratio, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  theme_minimal()


plot_raw.ssd <- 
  ggplot(data_soft_full, aes (site, ssd_dimensional, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  theme_minimal()


plot_raw.plant_height <- 
  ggplot(data_soft_full, aes (site, plant_height, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  ylim (0, 125) + # lubb prgl2 = outliers (trees vs. saps)
  theme_minimal()


plot_raw.c_delta.13 <- ## plot created using the combined df_soft
  ggplot(data_soft.spcomp, aes (site, c_delta.13, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  facet_grid(rows = vars (photo_pathway), scales = "free_y") +
  theme_minimal()


plot_raw.n_delta.15 <- ### LUBB HAS A SIGNIFICANT VALUE BETWEEN TREATMENT
  ggplot(data_soft_full, aes (site, n_delta.15, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  theme_minimal()


plot_raw.ldmc <-
  ggplot(data_soft_full, aes (site, ldmc, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  theme_minimal()
# temple has some weirdly high ldmc, likely due to the scale at temple not being
# able to register below 0.01 grams. it is what it is though.


### saving the plots ----------

# ## raw 
# png("../figures/plot_raw.leaf_thickness.png",
#     width = 8, height = 8, units = 'in', res = 1500)
# plot_raw.leaf_thickness
# dev.off()
# 
# png("../figures/plot_raw.sla.png",
#     width = 8, height = 8, units = 'in', res = 1500)
# plot_raw.sla
# dev.off()
# 
# png("../figures/plot_raw.n_concentraiton.png",
#     width = 8, height = 8, units = 'in', res = 1500)
# plot_raw.n_concentraiton
# dev.off()
# 
# png("../figures/plot_raw.cn_ratio.png",
#     width = 8, height = 8, units = 'in', res = 1500)
# plot_raw.cn_ratio
# dev.off()
# 
# png("../figures/plot_raw.ssd.png",
#     width = 8, height = 8, units = 'in', res = 1500)
# plot_raw.ssd
# dev.off()
# 
# png("../figures/plot_raw.plant_height.png",
#     width = 8, height = 8, units = 'in', res = 1500)
# plot_raw.plant_height
# dev.off()
# 
# png("../figures/plot_raw.c_delta.13.png",
#     width = 8, height = 8, units = 'in', res = 1500)
# plot_raw.c_delta.13
# dev.off()
# 
# png("../figures/plot_raw.n_delta.15.png",
#     width = 8, height = 8, units = 'in', res = 1500)
# plot_raw.n_delta.15
# dev.off()
# 
# png("../figures/plot_raw.ldmc.png",
#     width = 8, height = 8, units = 'in', res = 1500)
# plot_raw.ldmc
# dev.off()
#
# png("../figures/plot_raw.sla_errorbars.png",
#     width = 8, height = 8, units = 'in', res = 1500)
# plot_raw.sla_errorbars
# dev.off()


