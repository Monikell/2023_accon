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
library(BiodiversityR)
library(ggvegan)
# install.packages("remotes")
# remotes::install_github("gavinsimpson/ggvegan")




### load data ------------------------------------------------------------------
## data
# soft data = field and lab data collected
raw_data_soft <- read.csv("../data/02_cleaned/field/accon_field.lab_2023.csv")
raw_data_species.comp <- read.csv("../data/02_cleaned/field/accon_species.comp_2023.csv")

## metadata
metadata_plots <- read.csv("../data/00_meta/plot-descriptions-02-August-2019.csv")
metadata_species.list <- read.csv("../data/00_meta/species_list.csv")


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
df_soft_full <- data_soft_calc

## looking over data for qc/qa
# write.csv(df_soft_full, "../data/03_rproducts/df_soft_full.csv")
# MRK: looks good, except sample "sevi_39_boer_1_491", the CN isotope values
# seemed to have an issue from the lab, removing data error. 

# removing CN isotope issue sample
df_soft_full <- df_soft_full[df_soft_full$id_full != "sevi_39_boer_1_491", ]


################################################################################
## soft traits - models
################################################################################




################################################################################
## community weighted means (CWM) - models
################################################################################

### merging dfs soft and species comp together for CWM -------------------------

## verifying photosynthesis pathways based on delta c 13 values
# Bender 1971; Smith & Epstein 1971: C3 = 13C (-24‰ to -34‰), c4 = (-6‰ to -19‰)
df_soft_full <- df_soft_full %>%
  mutate(photo_pathway_delta.c.13 = case_when(
    c_delta.13 >= -35 & c_delta.13 <= -23 ~ "c3",
    c_delta.13 >= -19 & c_delta.13 <= -6 ~ "c4"
  ))
# MRK: looks good, only issue is one "boer", that has a weird C13 value, lab
# noted that the sample was very low, and that likely explains the value.


## mergining soft and species comp. average cover
df_soft_spcomp <- left_join(df_soft_full, df_grouped_spcomp_avg.cover, 
                   by = c("site_code", "plot", "taxon_code"))

## mergining species list information into dataset 
df_soft_spcomp <- left_join(df_soft_spcomp, metadata_species.list, 
                         by = c("site_code", "taxon_code"))

## removing site.y, and changing site.x to just site, for simplicity
# remove column site.y
df_soft_spcomp$site.y <- NULL

# update name 
names(df_soft_spcomp)[names(df_soft_spcomp) == "site.x"] <- "site"


## sp.comp % cover update NA values from "NA" to "0", top 5 species data
# collected from all plots even if not observed in 1x1s
df_soft_spcomp$average.cover[is.na(df_soft_spcomp$average.cover)] <- 0


## qa/qc real quick
# write.csv(df_soft_spcomp, "../data/03_rproducts/data_soft.spcomp.csv")
# MRK: looks good/ okay! 


### CWM: SLA -------------------------------------------------------------------

## looking for NA values
any(is.na(df_soft_spcomp$sla))

## looking over raw data
hist(df_soft_spcomp$sla)
hist(log(df_soft_spcomp$sla))

qqnorm(df_soft_spcomp$sla) ; qqline(df_soft_spcomp$sla, col = "red")
qqnorm(log(df_soft_spcomp$sla)) ; qqline(log(df_soft_spcomp$sla))


## cwm calculation
summarize_cwm_sla <- df_soft_spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_sla = log(weighted.mean(sla, average.cover)))


## model (nested rdmn) _________________________________________________________ HELP, 'is singluar'
# overcomplicated random effects, potentailyl remove block?
lmer_sla <- 
  lmer(cwm_sla ~ treatment * site + (1 | block:site), 
       data = summarize_cwm_sla)

# # Simplified model with random effect for site only
# lmer_sla_simplified <- lmer(cwm_sla ~ treatment * site + (1 | site), 
#                             data = summarize_cwm_sla)


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
ggplot(df_soft_spcomp, aes (site, sla, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1)

# cwm data
# MRK: arch and temple (wet sites) the SLA is significant
ggplot(summarize_cwm_sla, aes (site, cwm_sla, fill = treatment)) +
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
hist(df_soft_full$leaf_thickness)
hist(log(df_soft_full$leaf_thickness))
qqnorm(df_soft_full$leaf_thickness)
qqline(df_soft_full$leaf_thickness)
qqnorm(log(df_soft_full$leaf_thickness))
qqline(log(df_soft_full$leaf_thickness))

## looking for NA values
any(is.na(df_soft_spcomp$leaf_thickness))
unique(df_soft_spcomp$leaf_thickness)

## cwm calculation
summarize_leaf.thickness <- df_soft_spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_leaf.thickness = 
              log(weighted.mean(leaf_thickness, average.cover)))

## model (nested rdmn)
lmer_leaf.thickness <- 
  lmer(cwm_leaf.thickness ~ treatment * site + (1 | block:site), 
       data = summarize_leaf.thickness)

## look to see if a more simple model, get rid of 'singular' issue would be 
## resolved and make a better model, it doesn't look like it. 
# ## model (nested rdmn)
# lmer_leaf.thickness_site <- 
#   lmer(cwm_leaf.thickness ~ treatment * site + (1 | site), 
#        data = summarize_leaf.thickness)
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
ggplot(df_soft_spcomp, aes (site, leaf_thickness, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1)

# cwm data
ggplot(summarize_leaf.thickness, aes (site, cwm_leaf.thickness, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1)


## q-q plot
residuals_lmer_leaf.thickness <- residuals(lmer_leaf.thickness)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_leaf.thickness)
qqline(residuals_lmer_leaf.thickness)


### CWM: C:N ratio -------------------------------------------------------------

## histo & qq
hist(df_soft_full$cn_ratio)
hist(log(df_soft_full$cn_ratio))
qqnorm(df_soft_full$cn_ratio)
qqline(df_soft_full$cn_ratio)
qqnorm(log(df_soft_full$cn_ratio))
qqline(log(df_soft_full$cn_ratio))

## looking for NA values
any(is.na(df_soft_spcomp$cn_ratio))
unique(df_soft_spcomp$cn_ratio)

## cwm calculation
summarize_cn_ratio <- df_soft_spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_cn_ratio = weighted.mean(cn_ratio, average.cover,
                                                  na.rm = TRUE))

## model (nested rdmn)
lmer_cn.ratio <- 
  lmer(cwm_cn_ratio ~ treatment * site + (1 | block:site), ## no singular error
       data = summarize_cn_ratio)

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
ggplot(summarize_cn_ratio, aes(x = treatment, y = cwm_cn_ratio, fill = treatment)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Effect of Treatment on C:N Ratio",
       x = "Treatment",
       y = "Community Weighted Mean C:N Ratio")

# Create a plot for the main effect of site
ggplot(summarize_cn_ratio, aes(x = site, y = cwm_cn_ratio, fill = site)) +
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
hist(df_soft_full$n_concentration)
hist(log(df_soft_full$n_concentration))
qqnorm(df_soft_full$n_concentration)
qqline(df_soft_full$n_concentration)
qqnorm(log(df_soft_full$n_concentration))
qqline(log(df_soft_full$n_concentration)) # ------------------------------------ HELP, kind of s shaped...

## looking for NA values
any(is.na(df_soft_spcomp$n_concentration))
unique(df_soft_spcomp$n_concentration)

## cwm calculation
summarize_n.concentration <- df_soft_spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_n.concentration = log(weighted.mean(n_concentration, 
                                                    average.cover, 
                                                    na.rm = TRUE)))

## model (nested rdmn)
lmer_n.concentration <- 
  lmer(cwm_n.concentration ~ treatment * site + 
         (1 | block:site), data = summarize_n.concentration)

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
hist(df_soft_full$ssd_dimensional)
hist(log(df_soft_full$ssd_dimensional)) # much better

qqnorm(df_soft_full$ssd_dimensional) ; qqline(df_soft_full$ssd_dimensional)
qqnorm(log(df_soft_full$ssd_dimensional))
qqline(log(df_soft_full$ssd_dimensional)) # much better

## looking for NA values
any(is.na(df_soft_spcomp$ssd_dimensional))
unique(df_soft_spcomp$ssd_dimensional)

## cwm calculation
summarize_ssd <- df_soft_spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_ssd = log(weighted.mean(ssd_dimensional, 
                                    average.cover, na.rm = TRUE)))

## model (nested rdmn)
lmer_ssd <- lmer(cwm_ssd ~ treatment * site + (1 | block:site), 
                 data = summarize_ssd)

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
hist(df_soft_full$plant_height)
hist(log(df_soft_full$plant_height))
qqnorm(df_soft_spcomp$plant_height) ; qqline(df_soft_spcomp$plant_height)
qqnorm(log(df_soft_spcomp$plant_height))
qqline(log(df_soft_spcomp$plant_height))


## looking for NA values
any(is.na(df_soft_spcomp$plant_height))
unique(df_soft_spcomp$plant_height)

## cwm calculation
summarize_plant.height <- df_soft_spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_plant.height = log(weighted.mean(plant_height, average.cover)))

## model (nested rdmn)
lmer_plant.height <- 
  lmer(cwm_plant.height ~ treatment * site + (1 | block:site), 
       data = summarize_plant.height)

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
hist(df_soft_full$ldmc)
hist(log(df_soft_full$ldmc))
qqnorm(df_soft_spcomp$ldmc) ; qqline(df_soft_full$ldmc)
qqnorm(log(df_soft_spcomp$ldmc)) ; qqline(log(df_soft_spcomp$ldmc))

## looking for NA values
any(is.na(df_soft_spcomp$ldmc))
unique(df_soft_spcomp$ldmc)

## log needed
hist(df_soft_spcomp$ldmc)
hist(log(df_soft_spcomp$ldmc))

## cwm calculation
summarize_ldmc <- df_soft_spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_ldmc = log(weighted.mean(ldmc, average.cover,
                                                 na.rm = TRUE)))

## model (nested rdmn)
lmer_ldmc <- 
  lmer(cwm_ldmc ~ treatment * site + (1 | block:site), 
       data = summarize_ldmc)

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

## histo & qq, = bimodal
hist(df_soft_spcomp$c_delta.13)
qqnorm(df_soft_spcomp$c_delta.13) ; qqline(df_soft_spcomp$c_delta.13) #--------- HELP, bimodial 
# MRK: data is bimodal, split between c3 and c4 plants.. what to do about it?? 


## subsetting data by photo pathway, as the values can vary based on c3 vs. c4
df_data_c3 <- subset(df_soft_spcomp, photo_pathway == "c3")
df_data_c4 <- subset(df_soft_spcomp, photo_pathway == "c4")

hist(df_data_c3$c_delta.13)
hist(df_data_c4$c_delta.13)

qqnorm(df_data_c3$c_delta.13) ; qqline(df_data_c3$c_delta.13)
qqnorm(df_data_c4$c_delta.13) ; qqline(df_data_c4$c_delta.13)


## looking for NA values
any(is.na(df_soft_spcomp$c_delta.13))
unique(df_soft_spcomp$c_delta.13) # negative values no log


## cwm calculation ----------
summarize_c_delta.13 <- df_soft_spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_c_delta.13 = weighted.mean(c_delta.13, average.cover, 
                                           na.rm = TRUE))

## model (nested rdmn)
lmer_c_delta.13 <- 
  lmer(cwm_c_delta.13 ~ treatment * site + (1 | block:site), 
       data = summarize_c_delta.13)


## cwm calculation - photo. pathway included -----------
summarize_c_delta.13_photo <- df_soft_spcomp %>%
  group_by(site, plot, treatment, block, photo_pathway) %>% 
  summarise(cwm_c_delta.13 = weighted.mean(c_delta.13, average.cover, 
                                           na.rm = TRUE))

## model (nested rdmn) - photo. pathway included
lmer_c_delta.13_photo <- 
  lmer(cwm_c_delta.13 ~ treatment * site * photo_pathway + (1 | block:site), 
       data = summarize_c_delta.13_photo)


## plot
plot(lmer_c_delta.13, which = 2)
plot(lmer_c_delta.13_photo, which = 2)

## anova
Anova(lmer_c_delta.13) # site *
Anova(lmer_c_delta.13_photo) # site, photo_pathway *** & site:photo_path

## emmeans
emmeans(lmer_c_delta.13, ~site)
emmeans(lmer_c_delta.13, ~treatment)
emmeans(lmer_c_delta.13, ~treatment*site)

emmeans(lmer_c_delta.13_photo, ~site)

# looking deeper, p-value
# temple close with a p of 0.0894
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


## CWM: δ15N -------------------------------------------------------------------

## histo
hist(df_soft_spcomp$n_delta.15)
qqnorm(df_soft_spcomp$n_delta.15) ; qqline(df_soft_spcomp$n_delta.15)

## looking for NA values
any(is.na(df_soft_spcomp$n_delta.15))
unique(df_soft_spcomp$n_delta.15) # negative values no log

## cwm calculation
summarize_n_delta.15 <- df_soft_spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_n_delta.15 = weighted.mean(n_delta.15, average.cover, 
                                           na.rm = TRUE))

## model (nested rdmn)
lmer_n_delta.15 <- 
  lmer(cwm_n_delta.15 ~ treatment * site + (1 | block:site), 
       data = summarize_n_delta.15)

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

## not sure if this is right? -------------------------------------------------- PCA 1 video
## need to add in metadata? Site, treatment, species?

## subset data, traits of interest
pca_soft.traits <- df_soft_spcomp[,c("leaf_thickness", "sla", "n_concentration", 
                              "cn_ratio", "ssd_dimensional", "plant_height", 
                              "c_delta.13", "n_delta.15", "ldmc")]

## qa/qc
# write.csv(pca_soft.traits, "../data/03_rproducts/pca_soft.traits.csv")

# removing NA values
any(is.na(pca_soft.traits))
pca_soft.traits_clean <- na.omit(pca_soft.traits)
any(is.na(pca_soft.traits_clean))

## running the analysis
pca1 <- rda(pca_soft.traits_clean)

## look at total inertia (variance) and that explained by each PC
pca1
summary(pca1)

## determine how many axes to retain for interpretation (BiodiversityR, package)
# MRK: focus on % > bs%, only PC1 has "1.00000"
PCAsignificance(pca1)

## plot the data, really basic plot
ordiplot(pca1)
ordiplot(pca1, type = "t") # t = text


### plot the PCA using the package ----------
## using package 'ggvegan'
autoplot(pca1, legend.position = "none") +
  xlab("PC1 (92%)") +
  ylab("PC2 (5%") +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", size = 0.8) +
  geom_vline(aes(xintercept=0), linetype = "dashed", size = 0.8) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))


### PCA with datacamp ---------------------------------------------------------- PCA 2 attempt (raw data)
# https://www.datacamp.com/tutorial/pca-analysis-r

install.packages("corrr")
library('corrr')

install.packages("ggcorrplot")
library(ggcorrplot)

install.packages("FactoMineR")
library("FactoMineR")

## data cleaning/ checks
# check for na values in the columns (also only numerical varaibles)
colSums(is.na(pca_soft.traits_clean)) 

## Normalize the data (is this an issue for later? not sure)
pca_soft.traits_clean_normalized <- scale(pca_soft.traits_clean)
head(pca_soft.traits_clean_normalized)

## applying pca
pca2 <- princomp(pca_soft.traits_clean_normalized)
summary(pca2)
# summary explained, 9 principal components have been generated (Comp.1 to 9)
# matches number of variables put in. Each component explains a percentage of
# the total variance in the data set. In the "Cumulative Proportion" area, the 
# first pca explains 28% of the total variance. 2 pca = 16% (total 44%) and so 
# on. This is great BUT what do they really mean?? Can be answered by the below.
# Which, is exploring how they related to each column using the loadings of 
# each principal component. 

## The Loading Matrix
pca2$loadings[,1:6] # Comp.1 - Comp.6 (cumulative 87%)

# The Loading Matrix shows that the first principal component has high positive
# values leaf_thickness, sla, n_concentration, n_delta.15 others are negative. 
# the example talks about is in response to like countires (row columns), here 
# not sure.. like just overall prefrence? 

### Visualization of the principal components ----------

## scree plot
# visualize the importance of each principal component and can be used to 
# determine the number of principal components to retain. 
fviz_eig(data.pca, addlabels = TRUE) ### LEMON








### plot the PCA ggplot -------------------------------------------------------- PCA 3 ChatGPT (kinda bad)


pca_fort <- fortify(pca1, axes = 1:2)


ggplot() +
  geom_point(data = subset(pca_fort, score == 'sites'),
             mapping = aes(x = PC1, y = PC2), 
             colour = "darkgray", 
             alpha = 0.5) +
  geom_segment(data = subset(pca_fort, score == 'species'),
               mapping = aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.03, "npc"), type = "closed"),
               colour = "darkgray", 
               size = 0.8) +
  geom_text(data = subset(pca_fort, score == 'species'), 
            mapping = aes(label = label, x = PC1 * 1.1, y = PC2 * 1.1)) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", 
              size = 0.8, colour = "gray") + 
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.8, colour = "gray") +
  xlab("PC1 (92%)") + 
  ylab("PC2 (5%)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ylim(c(-10, 20)) +
  xlim(c(-25, 25))





### PCA attempt for site, treatment, taxon code -------------------------------- HELP, PCA, with shaddy AI help..
## trying to plot the traits, sites, species (taxon_code)


## subset data, traits of interest
pca_soft.meta <- df_soft_spcomp[,c("site", "treatment", "taxon_code", 
                                "leaf_thickness", "sla", "c_total", "n_total",
                                  "ssd_dimensional", "plant_height", 
                                "c_delta.13", "n_delta.15", "ldmc")]

## removing NA values
any(is.na(pca_soft.meta))
pca_soft.meta_clean <- na.omit(pca_soft.meta)
any(is.na(pca_soft.meta_clean))


## PCA on trait data ## -------------------------------------------------------- LEMON
pca_result <- prcomp(pca_soft.meta_clean %>%
                       select(leaf_thickness, sla, c_total, n_total, 
                              ssd_dimensional, plant_height, 
                              c_delta.13, n_delta.15, ldmc),
                     center = TRUE,
                     scale. = TRUE)


## looking at stuff
pca_result
summary(pca_result)

## merging together
pca_df <- as.data.frame(pca_result$x) %>%
  bind_cols(pca_soft.meta_clean %>% select(site, treatment, taxon_code))


## PLOTTING: dots for species, across sites and their PCA, (okay)
ggplot(pca_df, aes (x = PC1, y = PC2)) +
  geom_point(aes(color = site), size = 3, alpha = 0.7) +
  labs(title = "PCA of Traits",
       x = paste("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)", sep = "")) +
  theme_minimal()


## PLOTTING: arrows for species... not super informative? (not great?)
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = site), size = 3, alpha = 0.7) +
  geom_segment(data = pca_df %>%
                 group_by(taxon_code) %>%
                 summarize(PC1 = mean(PC1), PC2 = mean(PC2)), 
               aes(x = 0, y = 0, xend = PC1, yend = PC2, color = taxon_code), 
               arrow = arrow(length = unit(0.3, "inches")), size = 0.7) +
  labs(title = "PCA of Traits with Species Arrows",
       x = paste("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)", sep = "")) +
  theme_minimal()


## PLOTTING: dots for species, color for site, and with trait arrows
pca_loadings <- as.data.frame(pca_result$rotation)  # PCA loadings (traits)

# Assuming you have trait names in rownames(pca_loadings) and loading columns are named PC1, PC2
pca_loadings$trait <- rownames(pca_loadings)
pca_loadings$PC1 <- pca_loadings$PC1
pca_loadings$PC2 <- pca_loadings$PC2


# 
# # Plot with arrows for traits and dots for taxon_code
# ggplot(pca_df, aes(x = PC1, y = PC2)) +
#   # Points for taxon_code, colored by site
#   geom_point(aes(color = site), size = 3, alpha = 0.7) +
#   
#   # Arrows for traits
#   geom_segment(data = pca_loadings, 
#                aes(x = 0, y = 0, xend = PC1, yend = PC2, color = trait), 
#                arrow = arrow(length = unit(0.3, "inches"), type = "closed"), 
#                size = 0.7) +
#   
#   # Add labels for traits
#   geom_text(data = pca_loadings, 
#             aes(label = trait, x = PC1 * 1.1, y = PC2 * 1.1), 
#             size = 3, hjust = 0.5, vjust = 0.5) +
#   
#   # Labels and theme
#   labs(title = "PCA of Traits with Trait Arrows",
#        x = paste("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)", sep = ""),
#        y = paste("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)", sep = "")) +
#   theme_minimal()



## trying with bigger arrows
scaling_factor <- 10  # Adjust this factor as needed

pca_loadings_scaled <- pca_loadings %>%
  mutate(PC1 = PC1 * scaling_factor,
         PC2 = PC2 * scaling_factor)

# Plot with scaled trait arrows
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  # Points for taxon_code, colored by site
  geom_point(aes(color = site), size = 3, alpha = 0.7) +
  
  # Arrows for traits with scaling
  geom_segment(data = pca_loadings_scaled, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "inches"), type = "closed"), 
               size = 0.5) +
  
  # Add labels for traits
  geom_text(data = pca_loadings_scaled, 
            aes(label = trait, x = PC1 * 1.1, y = PC2 * 1.1), 
            size = 3, hjust = 0.5, vjust = 0.5) +
  
  # Add zero lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  
  # Labels and theme
  labs(title = "PCA of Traits with Scaled Trait Arrows",
       x = paste("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)", sep = ""),
       y = paste("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)", sep = "")) +
  theme_minimal()


################################################################################
## species composition - calculating diversity indexes
################################################################################

### cleaning up dataset --------------------------------------------------------

## adding in plot meta data
spcomp_plots <- left_join(raw_data_species.comp, metadata_plots,
                          by = c("site_code","plot"))

## lowering column names and the treatment column Values
names(spcomp_plots) <- tolower(colnames(spcomp_plots))
spcomp_plots$treatment <- tolower(spcomp_plots$treatment)


## isolating plant only data
spcomp_plants <- subset(spcomp_plots, type_of_cover == "plant")


## square percent covers for diversity estimate
spcomp_plants$percent.cover.squared <- 
  (spcomp_plants$percent_cover/100) * (spcomp_plants$percent_cover/100)


## grouping by site: plot | treatment
spcomp_grouped.site.plot <- group_by(spcomp_plants, site_code, plot)
spcomp_grouped.site.treatment <- group_by(spcomp_plants, site_code, treatment)


### richness -------------------------------------------------------------------

## site, plot
spcomp_richness_site.plot <- 
  summarise(spcomp_grouped.site.plot, richness_plot = n_distinct(taxon_code))


## site, treatment
spcomp_richness_site.treatment <- 
  summarise(spcomp_grouped.site.treatment, richness_treatment = 
              n_distinct(taxon_code))


### diversity ------------------------------------------------------------------

## site, plot
spcomp_diversity_site.plot <- 
  summarise(spcomp_grouped.site.plot, 
            diversity_plot = 1/sum(percent.cover.squared))

## site, treatment
spcomp_diversity_site.treatment <- 
  summarise(spcomp_grouped.site.treatment, 
            diversity_treatment = 1/sum(percent.cover.squared))


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

## updating data name for simplicity
spcomp_data_4lmer <- spcomp_evenness_plots


### exploring data -------------------------------------------------------------

# raw
hist(spcomp_data_4lmer$diversity_plot)
hist(spcomp_data_4lmer$richness_plot)
hist(spcomp_data_4lmer$evenness)

# log
hist(log(spcomp_data_4lmer$diversity_plot))
hist(log(spcomp_data_4lmer$richness_plot))
hist(log(spcomp_data_4lmer$evenness))



### model 01: p.comp diversity across sites ------------------------------------

## model
mod_div.site.trt <-lmer(log(diversity_plot) ~ sitefac * trtfac + 
                          (1 | plotfac) + (1 | block), 
                        data = (spcomp_data_4lmer))


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
## species composition - diversity index
################################################################################








################################################################################
## data visualization
################################################################################

### raw data -------------------------------------------------------------------

plot_raw.leaf_thickness <- 
  ggplot(df_soft_full, aes (site, leaf_thickness, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) +
  theme_minimal()


plot_raw.sla <- 
  ggplot(df_soft_full, aes (site, sla, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  theme_minimal() # temple, outliers might be okay..


## attempting with stats..SLA, looks amazing and has error bars!!
plot_raw.sla_errorbars <- ggplot() + 
    # Add error bars to the box plot for the specified subset of data
    stat_boxplot(
      data = df_soft_spcomp,
      aes(x = site, y = sla), 
      size = 0.75, 
      geom = "errorbar", 
      width = 0.2
    ) +
    # Add the box plot itself for the same subset of data
    geom_boxplot(
      data = df_soft_spcomp,
      aes(x = site, y = sla), 
      outlier.shape = NA
    ) + 
    theme_minimal()


## works okay facet by site, and treatment as the x.. works good
## likely looks better for other things/ sites
ggplot() + 
  # Add error bars to the box plot for the specified subset of data
  stat_boxplot(data = df_soft_spcomp,
    aes(x = treatment, y = sla), 
    size = 0.75, 
    geom = "errorbar", 
    width = 0.2 ) +
  # Add the box plot itself for the same subset of data
  geom_boxplot(data = df_soft_spcomp,
    aes(x = treatment, y = sla, fill = treatment), 
    outlier.shape = NA) + 
  # Facet by site (one row per site)
  facet_grid(rows = vars(site), scales = "free_y") +
  # Minimal theme
  theme_minimal()



plot_raw.n_concentraiton <- 
  ggplot(df_soft_full, aes (site, n_concentration, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  theme_minimal()


plot_raw.cn_ratio <- 
  ggplot(df_soft_full, aes (site, cn_ratio, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  theme_minimal()


plot_raw.ssd <- 
  ggplot(df_soft_full, aes (site, ssd_dimensional, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  theme_minimal()


plot_raw.plant_height <- 
  ggplot(df_soft_full, aes (site, plant_height, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  ylim (0, 125) + # lubb prgl2 = outliers (trees vs. saps)
  theme_minimal()


plot_raw.c_delta.13 <- ## plot created using the combined df_soft
  ggplot(df_soft_spcomp, aes (site, c_delta.13, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  facet_grid(rows = vars (photo_pathway), scales = "free_y") +
  theme_minimal()


plot_raw.n_delta.15 <- ### LUBB HAS A SIGNIFICANT VALUE BETWEEN TREATMENT
  ggplot(df_soft_full, aes (site, n_delta.15, fill = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, size = 2, color = "blue", alpha = 0.1) + 
  theme_minimal()


plot_raw.ldmc <-
  ggplot(df_soft_full, aes (site, ldmc, fill = treatment)) +
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


