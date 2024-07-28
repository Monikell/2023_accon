################################################################################
## Author: Monika Kelley
## Date: 2024/05/16 to _______
## Purpose: Summer 2023, Acquisitive vs. Conservative Grasslands Communities
## Project code/ name: ACCON (AC-quisitive CON-servative)
################################################################################


### load libraries -------------------------------------------------------------

library(tidyverse)
library(dplyr)
library(lme4)
library(ggplot2)
library(car) ## "Anova"
library(emmeans)



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
df_soft_full <- df_soft_calc.



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

# creating df for metadata plot info. to merge with plot data
metadata_plot_treatment <- 
  metadata_plots[,c("site_code","Treatment","plot","block","N","P","K")]
names(metadata_plot_treatment) <- tolower(colnames(metadata_plot_treatment))

# creating binary for treatment type
metadata_plot_treatment <- metadata_plot_treatment %>%
  mutate(treatment_binary = 
           ifelse(treatment == "NPK" | treatment == "NPK+Fence", 1, 0))

# merging diversity for plots with metadata
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

## add parameters as factors
## mk: probably just need to have either NPK or control.
spcomp_evenness_plots$nfac <- as.factor(spcomp_evenness_plots$n)
spcomp_evenness_plots$pfac <- as.factor(spcomp_evenness_plots$p)
spcomp_evenness_plots$kfac <- as.factor(spcomp_evenness_plots$k)
spcomp_evenness_plots$plotfac <- as.factor(spcomp_evenness_plots$plot)
spcomp_evenness_plots$sitefac <- as.factor(spcomp_evenness_plots$site_code)

## updating data name for simplicity
spcomp_data_4lmer <- spcomp_evenness_plots


######################################################################### LEMON

## got rid of the flugg extra treatment types


## add in treatments as binary factors, only focusing on NPK or control 
## trying to focus just on treatment with this plot..
spcomp_evenness_trtbinary <- spcomp_evenness_plots

spcomp_evenness_trtbinary$trtfac <- as.factor(spcomp_evenness_trtbinary$treatment_binary)
spcomp_evenness_trtbinary$plotfac <- as.factor(spcomp_evenness_trtbinary$plot)
spcomp_evenness_trtbinary$sitefac <- as.factor(spcomp_evenness_trtbinary$site_code)
spcomp_evenness_trtbinary$treatmentfac <- as.factor(spcomp_evenness_trtbinary$treatment)

## updating data name for simplicity/ to keep track! 
spcomp_data_4lmer_LEMON <- spcomp_evenness_trtbinary

######################################################################### LEMON


## statistical models of diversity across sites --------------------------------

## trying with first 
mod_div.site.trt <-lmer(log(diversity_plot) ~ sitefac * nfac * pfac * kfac + 
                          (1 | plotfac) + (1 | block), 
                        data = (spcomp_data_4lmer))


## residual plot
plot(resid(mod_div.site.trt) ~fitted(mod_div.site.trt))


## anova info
Anova(mod_div.site.trt)


## eemeans - LEMON Stopped here, issue with eemeans

######################################################################### LEMON

## got rid of the unnecessary extra treatment types

mod_div.site.trt_LEMON <-lmer(log(diversity_plot) ~ sitefac * treatmentfac + 
                          (1 | plotfac) + (1 | block), 
                        data = (spcomp_data_4lmer_LEMON))

### try smaller model for "isSingular" issue
# https://stackoverflow.com/questions/60028673/lme4-error-boundary-singular-fit-see-issingular


## residual plot
plot(resid(mod_div.site.trt_LEMON) ~fitted(mod_div.site.trt_LEMON))


## anova info
Anova(mod_div.site.trt_LEMON)


######################################################################### LEMON

