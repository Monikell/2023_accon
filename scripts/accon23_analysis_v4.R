################################################################################
## Author: Monika Kelley
## Date: 2024/05/16 to _______
## Purpose: Summer 2023, Acquisitive vs. Conservative Grasslands Communities
## Project code/ name: ACCON (AC-quisitive CON-servative)
################################################################################

### ctrl f notes 
# LEMON = last stopped
# HELP = issue/ not sure

setwd("C:/Users/monik/Documents/git/2023_accon/scripts")

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

### load data ------------------------------------------------------------------
## data
# soft data = field and lab data collected
soft_raw <- read.csv("../data/02_cleaned/field/accon_field.lab_2023.csv")
spcomp_raw <- read.csv("../data/02_cleaned/field/accon_species.comp_2023.csv")

## metadata
metadata_plots <- read.csv("../data/00_meta/plot-descriptions-02-August-2019.csv")
metadata_species.list <- read.csv("../data/00_meta/species_list.csv")

## r code for figures and things
# color pallet for figures
colors <- c("Control" = "#E0FFFF", "NPK" = "#7A8B8B")



################################################################################
## cleaning data
################################################################################

### soft data ------------------------------------------------------------------

## renaming data set to be modified
soft_calculations <- soft_raw

## lowering col names to make merging later easier
names(soft_calculations) <- tolower(colnames(soft_calculations))

# creating shorter version of plot metadata
metadata_plot_treatment <- 
  metadata_plots[,c("site_code","site","Treatment","plot","block")]

# easier to merge col. names
names(metadata_plot_treatment) <- tolower(colnames(metadata_plot_treatment))

## adding in plot information
soft_calculations <- left_join(soft_calculations, metadata_plot_treatment)


### species comp ---------------------------------------------------------------
## merging plot information with species comp information
spcomp_avg.cover <- 
  left_join(spcomp_raw, metadata_plot_treatment)

## removing the "non-plant" cover percentages
unique(spcomp_avg.cover$type_of_cover) 
spcomp_avg.cover <- subset(spcomp_avg.cover, type_of_cover == "plant")

## species comp. averages by site, plot, and per species (taxon_code)
# also condenses Temple plots into one subplot average.
spcomp_avg.cover <- spcomp_avg.cover %>%
  group_by(site_code, plot, taxon_code) %>%
  summarise(average.cover = mean(percent_cover, na.rm = TRUE))

## adding sp. comp data to the soft data. might be better to add 
## to average values per taxon code.. later? 
# soft_calculations <- left_join(soft_calculations, spcomp_avg.cover)



################################################################################
## soft traits - calculating soft plant traits

# traits: leaf thickness, SLA, carbon concentration, nitrogen concentration, 
# C:N ratio, SSD, plant height, δ13C, δ15N, LDMC 
# Coming: phosphorous and potassium
################################################################################


### specific leaf area (SLA) ---------------------------------------------------
## calculating SLA: area/oven-dry mass
soft_calculations$sla_cm2.g <- 
  soft_calculations$leaf_area / soft_calculations$leaf_dry_g

## convert sla cm^2/g to m^2/kg 
soft_calculations$sla_m2.kg <- soft_calculations$sla_cm2.g / 10


### leaf dry-matter content (LDMC) ---------------------------------------------
# expressed in units of mg g^-1, measured in grams conversion needed first. 

## dry weight g to mg 
soft_calculations$leaf_dry_mg <- soft_calculations$leaf_dry_g * 1000

## calculating LDMC: oven-dry mass / water-saturated fresh mass: units mg.g^-1 
soft_calculations$ldmc <- 
  soft_calculations$leaf_dry_mg / soft_calculations$leaf_wet_g


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
soft_calculations <- soft_calculations %>%
  rename(
    stem_hollow.or.solid = stem_hollow, # clarifying if stem solid or hollow
    stem_imagej_area.whole = stem_area_whole,
    stem_imagej_area.hollow = stem_area_hollow,
    ddh_average_mm = ddh_average
  )


## step 1.b: calculating average % of hollow stem for each species by site
df_stem_hollow.average <- soft_calculations %>%
  group_by(site, taxon_code) %>%
  summarise(stem_hollow.percent.average = mean(stem_area_percent_hollow, 
                                               na.rm = TRUE))

## step 1.c: merge the calculated averages back into the original data set
soft_calculations <- soft_calculations %>%
  left_join(df_stem_hollow.average, by = c("site", "taxon_code"))

## step 1.d: estimating hollow diameter
# calculate the solid and hollow area of the stem using the ddh average
soft_calculations$stem_total.area <- 
  pi * (soft_calculations$ddh_average_mm / 2)^2

# calculate the hollow area based on the calculated percentage
soft_calculations$stem_hollow.area <- 
  (soft_calculations$stem_hollow.percent.average / 100) * 
  soft_calculations$stem_total.area

# calculate the estimated diameter of the hollow area
# purpose: used to check math against field ddh values (reasonable numbers)
soft_calculations$stem_hollow.diameter_estimate <- 
  2 * sqrt(soft_calculations$stem_hollow.area / pi)


## step 2: finalizing stem area ----------

# step 2.a: replacing NA values with 0 for the stem_hollow.area
# purpose: highlights solid vs. hollow stems + good for math in step 2.b
soft_calculations <- soft_calculations %>% 
  mutate(stem_hollow.area = replace_na(stem_hollow.area, 0))


## step 2.b: calculating solid cross-sectional area: whole.area - hollow.area
## last step for dealing with hollow stems! can move on to volume calculations 
soft_calculations$stem_area.solid <- 
  soft_calculations$stem_total.area - soft_calculations$stem_hollow.area


## step 3: calculating volume  ----------
# formula V = (stem_area.solid) X length of stem, units typically in 
# mgmm–3

# converting length from cm to mm
soft_calculations$stem_length_mm <- soft_calculations$stem_length_cm * 10

## converting the area from mm^2 to cm^2 MRK: want to be in mm^2
# soft_calculations$stem_area.solid_cm2 <- soft_calculations$stem_area.solid/100

# calculating volume: V = (0.5D)^2 x pi x L
soft_calculations$stem_volume.dimensional_mm3 <- 
  soft_calculations$stem_area.solid * soft_calculations$stem_length_mm

# converting stem dry weight g to mg
soft_calculations$stem_dry_mg <- 
  soft_calculations$stem_dry_g * 1000

## step 4: calculating SSD
soft_calculations$ssd_dimensional <- 
  soft_calculations$stem_dry_mg / soft_calculations$stem_volume.dimensional_mm3


## 02) water-displacement method calculations ---------- (might be off/ wrong)
# small samples, = needed more reps? Not really sure what doing wrong here. 
# use dimensional values, has full data anyway (mostly)


## step 2: calculating ssd 
# volume physically calculated in lab, no math required other than ssd calc. 
 soft_calculations$ssd_displacement <- 
   soft_calculations$stem_dry_mg / soft_calculations$stem_water.displaced.g


### carbon and nitrogen concentration ------------------------------------------
## total of sample (µg)/ sample weight (mg)

# carbon
soft_calculations$c_concentration <- 
  soft_calculations$c_total / soft_calculations$cn_sample.weight

# nitrogen
soft_calculations$n_concentration <- 
  soft_calculations$n_total / soft_calculations$cn_sample.weight


### carbon nitrogen ration (c:n) -----------------------------------------------
## carbon concentration / nitrogen concentration

soft_calculations$cn_ratio <- 
  soft_calculations$c_concentration / soft_calculations$n_concentration


### soft plant data calculated creating new data frame  ------------------------
soft_full <- soft_calculations



### merging values for community weighted means (cwm) --------------------------

spcomp_avg.cover

soft_spcomp_full <- left_join(soft_full, metadata_species.list, 
                              by = c("site", "taxon_code"))

# creating a new column for merging that is just "site".
spcomp_avg.cover$site <- gsub("\\.us","",spcomp_avg.cover$site_code)

# merging the soft full, with all the sp. info with the average cover info.
soft_spcomp_full <- left_join(soft_spcomp_full, spcomp_avg.cover, 
                              by = c("site", "plot", "taxon_code"))

################################################################################
## soft traits - models (raw and cwm)
################################################################################

## stopped here!! Create modesl for both the raw and cwm, and create figures
## and plots, follow tempalte like that was used at ESA. 

### SLA ------------------------------------------------------------------------

## SLA raw ----------

# visualize data
hist(soft_spcomp_full$sla_m2.kg)
hist(log(soft_spcomp_full$sla_m2.kg))
qqnorm(log(soft_spcomp_full$sla_m2.kg))
qqline(log(soft_spcomp_full$sla_m2.kg))
# log needed

# model
sla_lmer_raw <- lmer(log(sla_m2.kg) ~ treatment * site + (1|plot) + (1|block),
                     data = (soft_spcomp_full))

# residual plots
plot(resid(sla_lmer_raw) ~fitted(sla_lmer_raw))
plot(sla_lmer_raw, which = 2)

# anova
Anova(sla_lmer_raw) # site ***

# emmeans 
emmeans(sla_lmer_raw, ~site)
cld(emmeans(sla_lmer_raw, ~treatment))


## plotting
fig_sla_raw <- ggplot(soft_spcomp_full, aes(x = treatment, y = sla_m2.kg, 
                                            fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("sla raw") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(1, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  facet_wrap( ~ site, scale = "free_y")

## save the image
# png('../figures/fig_sla_raw.png',
#     width = 12, height = 8, units = 'in', res = 1500)
# fig_sla_raw
# dev.off()


fig_sla_raw.sites <- ggplot(soft_spcomp_full, aes(x = treatment, y = sla_m2.kg, 
                                            fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("sla raw sites combined") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(1, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors)


png('../figures/fig_sla_raw.sites.png',
    width = 12, height = 8, units = 'in', res = 1500)
fig_sla_raw.sites
dev.off()


## SLA cwm ----------


## percentage of NA to zero, top 5 individuals 0% cover in some plots
# will use for all community weighted means
cwm_soft_spcomp <- soft_spcomp_full
cwm_soft_spcomp$average.cover[is.na(cwm_soft_spcomp$average.cover)] <- 0

cwm_spcomp_averages <- subset(cwm_soft_spcomp, 
                              select = c("site", "plot", "taxon_code", 
                                         "average.cover"))

## lemon want the average cover for each species, including zeros 
cwm_spcomp_averages <- cwm_spcomp_averages %>%
  group_by(site, plot, taxon_code) %>%
  summarise(average.cover = mean(average.cover, na.rm = TRUE))
                          
## average sla
sla_cwm <- cwm_soft_spcomp %>%
  group_by(site_code, plot, taxon_code) %>%
  summarise(sla_avg = mean(sla_m2.kg, na.rm = TRUE))

## merging back in plot information and the cover percentages  
sla_cwm <- left_join(sla_cwm, metadata_plot_treatment)
sla_cwm <- left_join(sla_cwm, cwm_spcomp_averages, 
                     by = c("site","plot","taxon_code"))

## cwm calculation ############################################################### LEMON NOT SURE ABOUT THE COMMUNITY WEIGHTED MEANS, the log??
sla_cwm <- sla_cwm %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(sla_avg = log(weighted.mean(sla_m2.kg, average.cover)))

hist(sla_cwm$average.cover)
hist(sla_cwm$sla_avg)

# visualize data
hist(soft_spcomp_full$sla_m2.kg)
hist(log(soft_spcomp_full$sla_m2.kg))
qqnorm(log(soft_spcomp_full$sla_m2.kg))
qqline(log(soft_spcomp_full$sla_m2.kg))
# log needed

# model
sla_lmer_raw <- lmer(log(sla_m2.kg) ~ treatment * site + (1|plot) + (1|block),
                     data = (soft_spcomp_full))

# residual plots
plot(resid(sla_lmer_raw) ~fitted(sla_lmer_raw))
plot(sla_lmer_raw, which = 2)

# anova
Anova(sla_lmer_raw) # site ***

# emmeans 
emmeans(sla_lmer_raw, ~site)
cld(emmeans(sla_lmer_raw, ~treatment))



# plotting
fig_sla_raw <- ggplot(soft_spcomp_full, aes(x = treatment, y = sla_m2.kg, 
                                            fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("sla raw") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(1, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  facet_wrap( ~ site, scale = "free_y")


png('../figures/fig_sla_raw.png',
    width = 12, height = 8, units = 'in', res = 1500)
fig_sla_raw
dev.off()


fig_sla_raw.sites <- ggplot(soft_spcomp_full, aes(x = treatment, y = sla_m2.kg, 
                                                  fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("sla raw sites combined") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(1, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors)


png('../figures/fig_sla_raw.sites.png',
    width = 12, height = 8, units = 'in', res = 1500)
fig_sla_raw.sites
dev.off()








## save this for the CN analysis
# removing CN isotope issue sample, extremely high (and likely bad/wrong) value
soft_full <- 
  soft_full[soft_full$id_full != "sevi_39_boer_1_491", ]

