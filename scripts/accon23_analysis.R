################################################################################
## Author: Monika Kelley
## Date: 2024/05/16 to _______
## Purpose: Summer 2023, Acquisitive vs. Conservative Grasslands Communities
## Project code/ name: ACCON (AC-quisitive CON-servative)
################################################################################

### Notes 
# LEMON = last stopped
# help = issue/ not sure

### load libraries -------------------------------------------------------------

library(tidyverse)
library(dplyr)
library(lme4)
library(ggplot2)
library(car) # "Anova"
library(emmeans)
library(multcompView) # "cld" for emmeans
library(multcomp) # "cld" for emmeans 


### load data ------------------------------------------------------------------
## data
# soft data = field and lab data collected
raw_data_soft <- read.csv("../data/02_cleaned/field/accon_field.lab_2023.csv")
raw_data_species.comp <- read.csv("../data/02_cleaned/field/accon_species.comp_2023.csv")

## metadata
metadata_plots <- read.csv("../data/00_meta/plot-descriptions-02-August-2019.csv")


################################################################################
## soft traits - calculating soft plant trait information
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


## getting sp.comp cover average
# MRK: HELP, make sure getting the average between plot is okay
 df_grouped_spcomp_avg.cover <- df_spcomp_avg.cover %>%
   group_by(site_code, plot, taxon_code) %>%
   summarise(average.cover = mean(percent_cover, na.rm = TRUE))


# mergining soft and species composition percent cover together
df_soft_calc. <- left_join(df_soft_calc., df_grouped_spcomp_avg.cover, 
                  by = c("site_code", "plot", "taxon_code"))


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


### exploring data -------------------------------------------------------------

## sla
hist(soft_full$sla)
hist(log(soft_full$sla))

## leaf thickness
hist(soft_full$leaf_thickness)
hist(log(soft_full$leaf_thickness))

## leaf carbon total
hist(soft_full$c_total) # looks normal
hist(log(soft_full$c_total)) # not needed?

## leaf carbon 13
hist(soft_full$c_delta.13)
# hist(log(soft_full$c_delta.13)) # "NaNs produced"

## leaf nitrogen total
hist(soft_full$n_total)
hist(log(soft_full$n_total))

## leaf nitrogen 15
hist(soft_full$n_delta.15)
hist(log(soft_full$n_delta.15))

## leaf dry-matter content
hist(soft_full$ldmc)
hist(log(soft_full$ldmc))

## stem-specific density
hist(soft_full$ssd_dimensional)
hist(log(soft_full$ssd_dimensional))

## plant height
hist(soft_full$plant_height)
hist(log(soft_full$plant_height))

## leaf angle
hist(soft_full$leaf.angle)
hist(log(soft_full$leaf.angle))


################################################################################
## soft traits - community weighted means (cwm)
################################################################################

### LEMON

## calculating cwm
summarize_cwm_sla <- soft_full %>%
  group_by(site_code, treatment, block) %>% 
  summarise(cwm_sla = weighted.mean(sla, average.cover))


################################################################################
## soft traits - PCA
################################################################################


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

# # updating metadata plot info. to merge with species comp. data
# metadata_plot_treatment <- 
#   metadata_plots[,c("site_code","Treatment","plot","block")]
# names(metadata_plot_treatment) <- tolower(colnames(metadata_plot_treatment))

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
# MRK: HELP, looks like some heteroscedasticity?


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
