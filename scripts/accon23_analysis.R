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
metadata_block.sevi.update <- read.csv("../data/00_meta/plot-descriptions-02-August-2019_sevi.block_updates.csv")


################################################################################
## soft traits - calculating soft plant trait information

# traits: leaf thickness, SLA, leaf carbon content, leaf nitrogen content, 
# SSD, plant height, δ13C, δ15N, LDMC (phosphorous and potassium to come!)
################################################################################

### cleaning data --------------------------------------------------------------

## renaming data set to be modified
df_soft_cleaned <- raw_data_soft

## lowering column names and the treatment column Values
names(df_soft_cleaned) <- tolower(colnames(df_soft_cleaned))

## renaming data set to be modified
df_soft_calc. <- df_soft_cleaned


### adding percent cover to the soft dataset ----------

## creating shorter version of plot metadata
metadata_plot_treatment <- 
  metadata_plots[,c("site_code","Treatment","plot","block")]

# easier to merge col. names
names(metadata_plot_treatment) <- tolower(colnames(metadata_plot_treatment))

## merging plot information with species comp information
df_spcomp_avg.cover <- left_join(raw_data_species.comp, metadata_plot_treatment)


## getting sp.comp cover averages
# MRK: sp.comp average between plots, good?? should be?? #---------------------- HELP sp.cover averages?
 df_grouped_spcomp_avg.cover <- df_spcomp_avg.cover %>%
   group_by(site_code, plot, taxon_code) %>%
   summarise(average.cover = mean(percent_cover, na.rm = TRUE))


## adding in updated sevi. block info ------------------------------------------ HELP (sevi block 1/3)
# MRK: check with Nick, not sure if a good solution? 
df_soft_calc. <- left_join(df_soft_calc., metadata_block.sevi.update)


### specific leaf area (SLA) ---------------------------------------------------
## calculating SLA: area/oven-dry mass
df_soft_calc.$sla <- 
  df_soft_calc.$leaf_area / df_soft_calc.$leaf_dry_g



### leaf dry-matter content (LDMC) ---------------------------------------------
## calculating LDMC: oven-dry mass / water-saturated fresh mass
df_soft_calc.$ldmc <- 
  df_soft_calc.$leaf_dry_g / df_soft_calc.$leaf_wet_g



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
df_soft_calc. <- df_soft_calc. %>%
  rename(
  stem_hollow.or.solid = stem_hollow, # clarifying if stem solid or hollow
  stem_imagej_area.whole = stem_area_whole,
  stem_imagej_area.hollow = stem_area_hollow,
  ddh_average_mm = ddh_average
  )


## step 1.b: calculating average % of hollow stem for each species by site
stem_hollow.average <- df_soft_calc. %>%
  group_by(site, taxon_code) %>%
  summarise(stem_hollow.percent.average = mean(stem_area_percent_hollow, 
                                               na.rm = TRUE))

## step 1.c: merge the calculated averages back into the original data set
df_soft_calc. <- df_soft_calc. %>%
  left_join(stem_hollow.average, by = c("site", "taxon_code"))


## step 1.d: estimating hollow diameter
# calculate the solid and hollow area of the stem using the ddh average
df_soft_calc.$stem_total.area <- 
  pi * (df_soft_calc.$ddh_average_mm / 2)^2

# calculate the hollow area based on the calculated percentage
df_soft_calc.$stem_hollow.area <- 
  (df_soft_calc.$stem_hollow.percent.average / 100) * 
  df_soft_calc.$stem_total.area

# calculate the estimated diameter of the hollow area
# purpose: used to check math against field ddh values (reasonable numbers)
df_soft_calc.$stem_hollow.diameter_estimate <- 
  2 * sqrt(df_soft_calc.$stem_hollow.area / pi)


## step 2: finalizing stem area ----------

# step 2.a: replacing NA values with 0 for the stem_hollow.area
# purpose: highlights solid vs. hollow stems + good for math in step 2.b
df_soft_calc. <- df_soft_calc. %>% 
  mutate(stem_hollow.area = replace_na(stem_hollow.area, 0))


## step 2.b: calculating solid cross-sectional area: whole.area - hollow.area
## last step for dealing with hollow stems! can move on to volume calculations 
df_soft_calc.$stem_area.solid <- 
  df_soft_calc.$stem_total.area - df_soft_calc.$stem_hollow.area



## step 3: calculating volume  ----------
# formula V = (stem_area.solid) X length of stem

# converting the area from mm^2 to cm^2
df_soft_calc.$stem_area.solid_cm2 <- df_soft_calc.$stem_area.solid / 100

# calculating volume: V = (0.5D)^2 x pi x L
df_soft_calc.$ stem_volume.dimensional <- 
  df_soft_calc.$stem_area.solid_cm2 * df_soft_calc.$stem_length_cm


## step 4: calculating SSD

df_soft_calc.$ssd_dimensional <- 
  df_soft_calc.$stem_dry_g / df_soft_calc.$stem_volume.dimensional


## 02) water-displacement method calculations ----------

## step 1: updating column names for clarity
df_soft_calc.$stem_volume.displacement <- df_soft_calc.$stem_water.displaced.g


## step 2: calculating ssd 
# volume physically calculated in lab, no math required other than ssd calc. 
df_soft_calc.$ssd_displacement <- 
  df_soft_calc.$stem_dry_g / df_soft_calc.$stem_volume.displacement



### soft plant data calculated creating new data frame  ------------------------
soft_full <- df_soft_calc.


################################################################################
## soft traits - community weighted means (CWM)
################################################################################

### merging dfs soft and species comp together for CWM -------------------------

## verifying photosynthesis pathways based on delta c 13 values
# Bender 1971; Smith & Epstein 1971: C3 = 13C (-24‰ to -34‰), c4 = (-6‰ to -19‰)
soft_full <- soft_full %>%
  mutate(photo_pathway_delta.c.13 = case_when(
    c_delta.13 >= -35 & c_delta.13 <= -23 ~ "c3",
    c_delta.13 >= -19 & c_delta.13 <= -6 ~ "c4"
  ))
# MRK: looks good, only issue is one "boer", that has a weird C13 value, lab
# noted that the sample was very low, and that likely explains the value.


## mergining soft and species comp. average cover
soft_spcomp <- left_join(soft_full, df_grouped_spcomp_avg.cover, 
                   by = c("site_code", "plot", "taxon_code"))

## mergining species list information into dataset 
soft_spcomp <- left_join(soft_spcomp, metadata_species.list, 
                         by = c("site_code", "taxon_code"))

## removing site.y, and changing site.x to just site, for simplicity
# remove column site.y
soft_spcomp$site.y <- NULL

# update name 
names(soft_spcomp)[names(soft_spcomp) == "site.x"] <- "site"


## sp.comp % cover update NA values from "NA" to "0", top 5 species data
# collected from all plots even if not observed in 1x1s
soft_spcomp$average.cover[is.na(soft_spcomp$average.cover)] <- 0


## looking over data
write.csv(soft_spcomp, "../data/03_rproducts/data_soft.spcomp.csv")



### CWM: SLA ------------------------------------------------------------------- HELP (sevi block 2/3)
# MRK: this is using the old block for sevi, all "block 1"

## looking for NA values
any(is.na(soft_spcomp$sla))

## histo
hist(soft_full$sla)
hist(log(soft_full$sla))

## cwm calculation
summarize_cwm_sla <- soft_spcomp %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cwm_sla = log(weighted.mean(sla, average.cover)))

## model (nested rdmn)
# MRK: "isSingular", model too complex/ not enough data #----------------------- HELP, "isSingular"
lmer_sla <- lmer(cwm_sla ~ treatment * site + (1 | block:site), 
                 data = summarize_cwm_sla)

## plot
plot(lmer_sla, which = 2)

## anova
Anova(lmer_sla) ## site ***, treatment:site *

## emmeans
emmeans(lmer_sla, ~site)


### CWM: SLA (with updated sevi block) ----------------------------------------- HELP (sevi block 3/3)
# MRK: this is using the new block for sevi
# MRK: used the "block_sevi.update" going forward with other traits.. ask Nick!

## looking for NA values
any(is.na(soft_spcomp$sla))

## cwm calculation
summarize_cwm_sla_block.update <- soft_spcomp %>%
  group_by(site, plot, treatment, block_sevi.update) %>% 
  summarise(cwm_sla = log(weighted.mean(sla, average.cover)))

## model (nested rdmn)
# MRK: "isSingular", model too complex?? 
lmer_sla_block.update <- 
  lmer(cwm_sla ~ treatment * site + (1 | block_sevi.update:site), 
       data = summarize_cwm_sla_block.update)

## plot
# MRK: cone structure... okay? some large outliers
plot(lmer_sla_block.update, which = 2)

## anova
Anova(lmer_sla_block.update) ## site ***, treatment:site ***

## emmeans
emmeans(lmer_sla_block.update, ~site) #temple highest emmean
emmeans(lmer_sla_block.update, ~treatment)
emmeans(lmer_sla_block.update, ~treatment*site)

# looking deeper, p-value (so close, p = 0.0659)
pairs(emmeans(lmer_sla, ~treatment, at = list(site = 'temple'))) #p 0.1061
pairs(emmeans(lmer_sla, ~treatment, at = list(site = 'arch'))) #p 0.1003
pairs(emmeans(lmer_sla, ~treatment, at = list(site = 'lubb')))
pairs(emmeans(lmer_sla, ~treatment, at = list(site = 'sevi')))


# visualize and interpret pairwise comparison (trt not sig. from each other)
cld(emmeans(lmer_sla, ~treatment, at = list (site = 'temple')))

## q-q plot
residuals_lmer_sla <- residuals(lmer_sla)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_sla)
qqline(residuals_lmer_sla)


### CWM: Leaf thickness --------------------------------------------------------

## histo
hist(soft_full$leaf_thickness)
hist(log(soft_full$leaf_thickness))
qqnorm(log(soft_full$leaf_thickness))
qqline(log(soft_full$leaf_thickness))

## looking for NA values
any(is.na(soft_spcomp$leaf_thickness))
unique(soft_spcomp$leaf_thickness)

## cwm calculation
summarize_leaf.thickness <- soft_spcomp %>%
  group_by(site, plot, treatment, block_sevi.update) %>% 
  summarise(cwm_leaf.thickness = 
              log(weighted.mean(leaf_thickness, average.cover)))

## model (nested rdmn)
lmer_leaf.thickness <- 
  lmer(cwm_leaf.thickness ~ treatment * site + (1 | block_sevi.update:site), 
       data = summarize_leaf.thickness)

## plot
plot(lmer_leaf.thickness, which = 2)

## anova
Anova(lmer_leaf.thickness) # site **

## emmeans
emmeans(lmer_leaf.thickness, ~site)
emmeans(lmer_leaf.thickness, ~treatment)
emmeans(lmer_leaf.thickness, ~treatment*site)

# looking deeper
pairs(emmeans(lmer_leaf.thickness, ~treatment, at = list(site = 'temple')))
pairs(emmeans(lmer_leaf.thickness, ~treatment, at = list(site = 'lubb')))

# visualize and interpret pairwise comparison (trt not sig. from each other)
cld(emmeans(lmer_leaf.thickness, ~treatment, at = list (site = 'temple')))

## q-q plot
residuals_lmer_leaf.thickness <- residuals(lmer_leaf.thickness)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_leaf.thickness)
qqline(residuals_lmer_leaf.thickness)


### CWM: leaf carbon total -----------------------------------------------------

## histo & qq
hist(soft_full$c_total) # looks alright
qqnorm(soft_full$c_total)
qqline(soft_full$c_total)

## looking for NA values
any(is.na(soft_spcomp$c_total))
unique(soft_spcomp$c_total)

## cwm calculation
summarize_leaf.carbon.total <- soft_spcomp %>%
  group_by(site, plot, treatment, block_sevi.update) %>% 
  summarise(cwm_leaf.carbon.total = weighted.mean(c_total, average.cover,
                                                  na.rm = TRUE))

## model (nested rdmn)
lmer_leaf.carbon.total <- 
  lmer(cwm_leaf.carbon.total ~ treatment * site + (1 | block_sevi.update:site), 
       data = summarize_leaf.carbon.total)

## plot
plot(lmer_leaf.carbon.total, which = 2) 

## anova
Anova(lmer_leaf.carbon.total) # site **

## emmeans
emmeans(lmer_leaf.carbon.total, ~site)
emmeans(lmer_leaf.carbon.total, ~treatment)
emmeans(lmer_leaf.carbon.total, ~treatment*site)

# looking deeper, p-value
pairs(emmeans(lmer_leaf.carbon.total, ~treatment, at = list(site = 'temple')))
pairs(emmeans(lmer_leaf.carbon.total, ~treatment, at = list(site = 'lubb')))
pairs(emmeans(lmer_leaf.carbon.total, ~treatment, at = list(site = 'sevi')))
pairs(emmeans(lmer_leaf.carbon.total, ~treatment, at = list(site = 'arch')))

# visualize and interpret pairwise comparison (trt not sig. from each other)
cld(emmeans(lmer_leaf.carbon.total, ~treatment, at = list (site = 'temple')))

## q-q plot
residuals_lmer_leaf.carbon.total <- residuals(lmer_leaf.carbon.total)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_leaf.carbon.total)
qqline(residuals_lmer_leaf.carbon.total)


### CWM: leaf nitrogen total ---------------------------------------------------

## leaf nitrogen total
hist(soft_full$n_total)
hist(log(soft_full$n_total))
qqnorm(soft_full$n_total)
qqline(soft_full$n_total)
qqnorm(log(soft_full$n_total))
qqline(log(soft_full$n_total))

## looking for NA values
any(is.na(soft_spcomp$n_total))
unique(soft_spcomp$n_total)

## cwm calculation
summarize_leaf.nitrogen.total <- soft_spcomp %>%
  group_by(site, plot, treatment, block_sevi.update) %>% 
  summarise(cwm_leaf.nitrogen.total = log(weighted.mean(n_total, average.cover, 
                                                    na.rm = TRUE)))

## model (nested rdmn)
lmer_leaf.nitrogen.total <- 
  lmer(cwm_leaf.nitrogen.total ~ treatment * site + 
         (1 | block_sevi.update:site), data = summarize_leaf.nitrogen.total)

## plot
plot(lmer_leaf.nitrogen.total, which = 2)

## anova
Anova(lmer_leaf.nitrogen.total) # site ***, treatment *

## emmeans
emmeans(lmer_leaf.nitrogen.total, ~site) # lubb with highest emmean
emmeans(lmer_leaf.nitrogen.total, ~treatment)
emmeans(lmer_leaf.nitrogen.total, ~treatment*site)

# looking deeper, p-value
pairs(emmeans(lmer_leaf.nitrogen.total, ~treatment, at = list(site = 'temple')))
pairs(emmeans(lmer_leaf.nitrogen.total, ~treatment, at = list(site = 'lubb')))
# MRK, sevi has 0.0557
pairs(emmeans(lmer_leaf.nitrogen.total, ~treatment, at = list(site = 'sevi')))
pairs(emmeans(lmer_leaf.nitrogen.total, ~treatment, at = list(site = 'arch')))

# visualize and interpret pairwise comparison (trt not sig. from each other)
cld(emmeans(lmer_leaf.nitrogen.total, ~treatment, at = list (site = 'sevi')))

## q-q plot
residuals_lmer_leaf.nitrogen.total <- residuals(lmer_leaf.nitrogen.total)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_leaf.nitrogen.total)
qqline(residuals_lmer_leaf.nitrogen.total)


### CWM: SSD  ------------------------------------------------------------------

## histo
hist(soft_full$ssd_dimensional)
hist(log(soft_full$ssd_dimensional))

## looking for NA values
any(is.na(soft_spcomp$ssd_dimensional))
unique(soft_spcomp$ssd_dimensional)

## cwm calculation
summarize_ssd <- soft_spcomp %>%
  group_by(site, plot, treatment, block_sevi.update) %>% 
  summarise(cwm_ssd = log(weighted.mean(ssd_dimensional, 
                                    average.cover, na.rm = TRUE)))

## model (nested rdmn)
lmer_ssd <- lmer(cwm_ssd ~ treatment * site + (1 | block_sevi.update:site), 
                 data = summarize_ssd)

## plot
plot(lmer_ssd, which = 2)

## anova
Anova(lmer_ssd) # none

## q-q plot
residuals_lmer_ssd <- residuals(lmer_ssd)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_ssd)
qqline(residuals_lmer_ssd)


### CWM: plant height ----------------------------------------------------------

## histo
hist(soft_full$plant_height)
hist(log(soft_full$plant_height))
qqnorm(soft_spcomp$plant_height)
qqnorm(log(soft_spcomp$plant_height))


## looking for NA values
any(is.na(soft_spcomp$plant_height))

## cwm calculation
summarize_plant.height <- soft_spcomp %>%
  group_by(site, plot, treatment, block_sevi.update) %>% 
  summarise(cwm_plant.height = log(weighted.mean(plant_height, average.cover)))

## model (nested rdmn)
lmer_plant.height <- 
  lmer(cwm_plant.height ~ treatment * site + (1 | block_sevi.update:site), 
       data = summarize_plant.height)

## plot
plot(lmer_plant.height, which = 2)

## anova
Anova(lmer_plant.height) # site ***, treatment .

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
qqnorm(residuals_lmer_plant.height)
qqline(residuals_lmer_plant.height)


### CWM: LDMC ------------------------------------------------------------------

## leaf dry-matter content
hist(soft_full$ldmc)
hist(log(soft_full$ldmc))
qqnorm(soft_spcomp$ldmc)
qqnorm(log(soft_spcomp$ldmc))
qqline(log(soft_spcomp$ldmc))

## looking for NA values
any(is.na(soft_spcomp$ldmc))
unique(soft_spcomp$ldmc)

## log needed
hist(soft_spcomp$ldmc)
hist(log(soft_spcomp$ldmc))

## cwm calculation
summarize_ldmc <- soft_spcomp %>%
  group_by(site, plot, treatment, block_sevi.update) %>% 
  summarise(cwm_ldmc = log(weighted.mean(ldmc, average.cover,
                                                 na.rm = TRUE)))

## model (nested rdmn)
lmer_ldmc <- 
  lmer(cwm_ldmc ~ treatment * site + (1 | block_sevi.update:site), 
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
qqnorm(residuals_lmer_ldmc)
qqline(residuals_lmer_ldmc)


## CWM: δ13C -------------------------------------------------------------------

## histo & qq, = bimodal
hist(soft_spcomp$c_delta.13) 
qqnorm(soft_spcomp$c_delta.13) 

## looking for NA values
any(is.na(soft_spcomp$c_delta.13))
unique(soft_spcomp$c_delta.13) # negative values no log


## cwm calculation
summarize_c_delta.13 <- soft_spcomp %>%
  group_by(site, plot, treatment, block_sevi.update) %>% 
  summarise(cwm_c_delta.13 = weighted.mean(c_delta.13, average.cover, 
                                           na.rm = TRUE))

## model (nested rdmn)
lmer_c_delta.13 <- 
  lmer(cwm_c_delta.13 ~ treatment * site + (1 | block_sevi.update:site), 
       data = summarize_c_delta.13)

## plot
plot(lmer_c_delta.13, which = 2)

## anova
Anova(lmer_c_delta.13) # site *

## emmeans
emmeans(lmer_c_delta.13, ~site)
emmeans(lmer_c_delta.13, ~treatment)
emmeans(lmer_c_delta.13, ~treatment*site)

# looking deeper, p-value
# temple close with a p of 0.0894
pairs(emmeans(lmer_c_delta.13, ~treatment, at = list(site = 'sevi')))
pairs(emmeans(lmer_c_delta.13, ~treatment, at = list(site = 'temple')))

# visualize and interpret pairwise comparison (trt not sig. from each other)
cld(emmeans(lmer_c_delta.13, ~treatment, at = list (site = 'temple')))

## q-q plot
residuals_lmer_c_delta.13 <- residuals(lmer_c_delta.13)
# Create Q-Q plot of residuals
qqnorm(residuals_lmer_c_delta.13)
qqline(residuals_lmer_c_delta.13)


## CWM: δ15N -------------------------------------------------------------------

## histo
hist(soft_spcomp$n_delta.15)
qqnorm(soft_spcomp$n_delta.15)
qqline(soft_spcomp$n_delta.15)

## looking for NA values
any(is.na(soft_spcomp$n_delta.15))
unique(soft_spcomp$n_delta.15) # negative values no log

## cwm calculation
summarize_n_delta.15 <- soft_spcomp %>%
  group_by(site, plot, treatment, block_sevi.update) %>% 
  summarise(cwm_n_delta.15 = weighted.mean(n_delta.15, average.cover))

## model (nested rdmn)
lmer_n_delta.15 <- 
  lmer(cwm_n_delta.15 ~ treatment * site + (1 | block_sevi.update:site), 
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

## not sure if this is right? -------------------------------------------------- HELP PCA
## need to add in metadata? Site, treatment, species?

## subset data, traits of interest
pca_soft.traits <- soft_spcomp[,c("leaf_thickness", "sla", "c_total", "n_total",
                                  "ssd_dimensional", "plant_height", 
                                  "c_delta.13", "n_delta.15", "ldmc")]

## looking for what values have NA
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


### plot the PCA ggplot ----------


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
        axis.line = element_line(colour = "black"))





### PCA attempt for site, treatment, taxon code -------------------------------- HELP, PCA, with shaddy AI help..
## trying to plot the traits, sites, species (taxon_code)


## subset data, traits of interest
pca_soft.meta <- soft_spcomp[,c("site", "treatment", "taxon_code", 
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

