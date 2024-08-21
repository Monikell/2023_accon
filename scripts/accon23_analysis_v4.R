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

## ggplot
library(ggpattern) # hatch marks
library(patchwork) # easily combine plots

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





################################################################################
## cleaning data
################################################################################

## r code for figures ("z")
# color pallet for figures
colors <- c("Control" = "#E0FFFF", "NPK" = "#7A8B8B")

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
# C:N ratio, SSD, plant height, δ13C (c_delta.13), δ15N, LDMC 
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
sla_raw_lmer <- lmer(log(sla_m2.kg) ~ treatment * site + (1|plot) + (1|block),
                     data = (soft_spcomp_full))

# residual plots
plot(resid(sla_raw_lmer) ~fitted(sla_raw_lmer))
plot(sla_raw_lmer, which = 2)

# anova
Anova(sla_raw_lmer) # site ***

# emmeans 
emmeans(sla_raw_lmer, ~site)
cld(emmeans(sla_raw_lmer, ~treatment))
pairs(emmeans(sla_raw_lmer, ~site))

## plotting
z_sla_raw_sites <- ggplot(soft_spcomp_full, aes(x = treatment, y = sla_m2.kg, 
                                            fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("sla individuals") +
  ylab(expression("m"^2*"kg")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  facet_grid( ~ site) +
  coord_cartesian(ylim = c(0, 30))


# # save image
# png('../figures/z_sla_raw_sites.png',
#     width = 9, height = 8, units = 'in', res = 1000)
#    z_sla_raw_sites
#    dev.off()


 
z_sla_raw_all <- ggplot(soft_spcomp_full, aes(x = treatment, y = sla_m2.kg, 
                                            fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("sla individuals") +
  ylab(expression("m"^2*"kg")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(1, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(0, 31))

# save image
 # png('../figures/z_sla_raw_all.png',
 #     width = 9, height = 8, units = 'in', res = 1000)
 # z_sla_raw_all
 # dev.off()


## SLA cwm ---------- 

## percentage of NA to zero, top 5 individuals 0% cover in some plots
# will use for all community weighted means
cwm_soft_spcomp <- soft_spcomp_full
cwm_soft_spcomp$average.cover[is.na(cwm_soft_spcomp$average.cover)] <- 0


cwm_spcomp_averages <- subset(cwm_soft_spcomp, 
                              select = c("site", "plot", "taxon_code", 
                                         "average.cover"))

## average cover for each species, including zeros 
cwm_spcomp_averages <- cwm_spcomp_averages %>%
  group_by(site, plot, taxon_code) %>%
  summarise(average.cover = mean(average.cover, na.rm = TRUE))



## average sla
sla_cwm <- cwm_soft_spcomp %>%
  group_by(site, plot, taxon_code) %>%
  summarise(sla_avg = mean(sla_m2.kg, na.rm = TRUE))

## merging back in plot information and the cover percentages  
sla_cwm <- left_join(sla_cwm, metadata_plot_treatment)
sla_cwm <- left_join(sla_cwm, cwm_spcomp_averages, 
                     by = c("site","plot","taxon_code"))

## cwm calculation 
sla_cwm_calculated <- sla_cwm %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(sla_cwm = weighted.mean(sla_avg, average.cover))


## visualizing
hist(log(sla_cwm_calculated$sla_cwm))
qqnorm(log(sla_cwm_calculated$sla_cwm))
qqline(log(sla_cwm_calculated$sla_cwm))
qqnorm(sla_cwm_calculated$sla_cwm)
qqline(sla_cwm_calculated$sla_cwm)
# log looks a little better

# model
sla_cwm_lmer <- lmer(log(sla_cwm) ~ treatment * site + (1|plot) + (1|block),
                     data = (sla_cwm_calculated))

# residual plots
plot(resid(sla_cwm_lmer) ~fitted(sla_cwm_lmer))
plot(sla_cwm_lmer, which = 2)

# anova
Anova(sla_cwm_lmer) # site ***

# emmeans 
emmeans(sla_cwm_lmer, ~site)
cld(emmeans(sla_cwm_lmer, ~treatment))

# plotting
z_sla_cwm_sites <- ggplot(sla_cwm_calculated, 
                          aes(x = treatment, y = sla_cwm, fill = treatment)) + 
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  ggtitle("sla cwm") +
  ylab(expression("m"^2*"kg")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) +
  facet_grid(~ site) +
  coord_cartesian(ylim = c(0, 30))



## all sites
z_sla_cwm_all <- ggplot(sla_cwm_calculated, aes(x = treatment, y = sla_cwm, 
                                                  fill = treatment)) +
  geom_boxplot_pattern(color = "black", outlier.color = "black", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  ggtitle("sla cwm") +
  ylab(expression("m"^2*"kg")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(1, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors)+
  coord_cartesian(ylim = c(0, 31))


# merging plots together
fig1_sla_sites <- z_sla_raw_sites / z_sla_cwm_sites
fig2_sla_all <- z_sla_raw_all + z_sla_cwm_all

# # saving the image - sites
png('../figures/fig1_sla_sites.png',
    width = 18, height = 8, units = 'in', res = 1000)
fig1_sla_sites
dev.off()

# saving the image
 png('../figures/fig2_sla_all.png',
     width = 15, height = 7, units = 'in', res = 1000)
 fig2_sla_all
 dev.off()




### C:N ratio ------------------------------------------------------------------

## CN raw ----------

# visualize data
hist(soft_spcomp_full$cn_ratio)
hist(log(soft_spcomp_full$cn_ratio))
qqnorm(soft_spcomp_full$cn_ratio)
qqline(soft_spcomp_full$cn_ratio)
qqnorm(log(soft_spcomp_full$cn_ratio))
qqline(log(soft_spcomp_full$cn_ratio))
# log needed

# model
cn_raw_lmer <- lmer(log(cn_ratio) ~ treatment * site + (1|plot) + (1|block),
                     data = (soft_spcomp_full))

# residual plots
plot(resid(cn_raw_lmer) ~fitted(cn_raw_lmer))
plot(cn_raw_lmer, which = 2)

# anova
Anova(cn_raw_lmer) # site ***, treatment *

# emmeans 
emmeans(cn_raw_lmer, ~site)
emmeans(cn_raw_lmer, ~treatment)
emmeans(cn_raw_lmer, ~site*treatment)
cld(emmeans(cn_raw_lmer, ~treatment))

pairs(emmeans(cn_raw_lmer, ~treatment)) # 0.03 !!
pairs(emmeans(cn_raw_lmer, ~treatment, at = list(site = 'temple'))) #0.21
pairs(emmeans(cn_raw_lmer, ~treatment, at = list(site = 'arch'))) #0.10
pairs(emmeans(cn_raw_lmer, ~treatment, at = list(site = 'lubb'))) #0.92
pairs(emmeans(cn_raw_lmer, ~treatment, at = list(site = 'sevi'))) #0.06


# figure sites
z_cn_ratio_raw_sites <- ggplot(soft_spcomp_full, 
                               aes(x = treatment, y = cn_ratio, 
                                   fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("C:N individual") +
  ylab("C:N ratio") +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  facet_grid( ~ site) + 
  coord_cartesian(ylim = c(0,80))


# figure all 
z_cn_ratio_raw_all <- ggplot(soft_spcomp_full, aes(x = treatment, y = cn_ratio, 
                                                fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("C:N individual") +
  ylab("C:N ratio") +
  xlab ("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(1, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(0, 70))

 
## CN cwm ---------- 

## percentage of NA to zero, top 5 individuals 0% cover in some plots
# will use for all community weighted means
cwm_soft_spcomp <- soft_spcomp_full
cwm_soft_spcomp$average.cover[is.na(cwm_soft_spcomp$average.cover)] <- 0

# new cover that includes all species, even those with percent cover "0", which
# should be the "top 5" species, as those were sampled for soft traits
# regardless of presence for each plot
cwm_spcomp_averages <- subset(cwm_soft_spcomp, 
                              select = c("site", "plot", "taxon_code", 
                                         "average.cover"))

# merging together "duplicate" values, those with replicates. 
cwm_spcomp_averages <- cwm_spcomp_averages %>%
  group_by(site, plot, taxon_code) %>%
  summarise(average.cover = mean(average.cover, na.rm = TRUE))


## average cn_ratio, per site, plot, & species
cn_ratio_cwm <- cwm_soft_spcomp %>%
  group_by(site, plot, taxon_code) %>%
  summarise(cn_ratio_avg = mean(cn_ratio, na.rm = TRUE))

## merging back in plot information and the cover percentages  
cn_ratio_cwm <- left_join(cn_ratio_cwm, metadata_plot_treatment)
cn_ratio_cwm <- left_join(cn_ratio_cwm, cwm_spcomp_averages, 
                     by = c("site","plot","taxon_code"))

## cwm calculation,
cn_ratio_cwm_calculated <- cn_ratio_cwm %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cn_ratio_cwm = 
              weighted.mean(cn_ratio_avg, average.cover, na.rm=TRUE))


## visualizing
hist(cn_ratio_cwm_calculated$cn_ratio_cwm)
hist(log(cn_ratio_cwm_calculated$cn_ratio_cwm))
qqnorm(log(cn_ratio_cwm_calculated$cn_ratio_cwm))
qqline(log(cn_ratio_cwm_calculated$cn_ratio_cwm))
qqnorm(cn_ratio_cwm_calculated$cn_ratio_cwm)
qqline(cn_ratio_cwm_calculated$cn_ratio_cwm)
# log looks a little better

# model
cn_ratio_cwm_lmer <- lmer(log(cn_ratio_cwm) ~ treatment *
                            site + (1|plot) + (1|block),
                     data = (cn_ratio_cwm_calculated))

# residual plots
plot(resid(cn_ratio_cwm_lmer) ~fitted(cn_ratio_cwm_lmer))
plot(cn_ratio_cwm_lmer, which = 2)

# anova
Anova(cn_ratio_cwm_lmer)

# emmeans 
cld(emmeans(cn_ratio_cwm_lmer, ~site))
cld(emmeans(cn_ratio_cwm_lmer, ~treatment))
cld(emmeans(cn_ratio_cwm_lmer, ~treatment*site))

pairs(emmeans(cn_ratio_cwm_lmer, ~treatment, at = list(site='temple'))) #0.0129
pairs(emmeans(cn_ratio_cwm_lmer, ~treatment, at = list(site='sevi'))) #0.0104
pairs(emmeans(cn_ratio_cwm_lmer, ~treatment, at = list(site='arch'))) #0.3624
pairs(emmeans(cn_ratio_cwm_lmer, ~treatment, at = list(site='lubb'))) #0.7572
pairs((emmeans(cn_ratio_cwm_lmer, ~treatment))) #0.0119


# plotting - sites
z_cn_ratio_cwm_sites <- ggplot(cn_ratio_cwm_calculated, 
                                 aes(x = treatment, y = cn_ratio_cwm, 
                                     fill = treatment)) +
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  ggtitle("C:N cwm") +
  ylab("C:N ratio") +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  facet_grid( ~ site) + 
  coord_cartesian(ylim = c(0, 80))


# plotting - all
z_cn_ratio_cwm_all <- ggplot(cn_ratio_cwm_calculated, 
                          aes(x = treatment, y = cn_ratio_cwm, 
                                                  fill = treatment)) +
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  ggtitle("C:N cwm") +
  ylab("C:N ratio") +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(1, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(0, 70))


# combining plots
fig3_cn_ratio_sites <- z_cn_ratio_raw_sites / z_cn_ratio_cwm_sites 
fig4_cn_ratio_all <- z_cn_ratio_raw_all + z_cn_ratio_cwm_all


# saving the plots
sites
# png('../figures/fig3_cn_ratio_sites.png',
#     width = 18, height = 8, units = 'in', res = 1000)
# fig3_cn_ratio_sites
# dev.off()
# 
# # all
# png('../figures/fig4_cn_ratio_all.png',
#     width = 15, height = 7, units = 'in', res = 1000)
# fig4_cn_ratio_all
# dev.off()


### δ13C, (c_delta.13) ---------------------------------------------------------

## δ13C, (c_delta.13) raw ----------
# bi-modal data, c3 and c4 plants have diff. delta 13 values, separating.
cd13_c3 <- subset(soft_spcomp_full, photo_pathway == "c3")
cd13_c4 <- subset(soft_spcomp_full, photo_pathway == "c4")

# removing issue sample, "contained less C than smallest reference"
cd13_c4 <- cd13_c4[cd13_c4$id_full != "sevi_39_boer_1_491", ]

# removing NA rows, data not sent/ issues with sample
cd13_c3 <- subset(cd13_c3, !is.na(c_delta.13))


# visualize data
hist(cd13_c3$c_delta.13) #c3
qqnorm(cd13_c3$c_delta.13) #c3
qqline(cd13_c3$c_delta.13) #c3

hist(cd13_c4$c_delta.13) #c4
qqnorm(cd13_c4$c_delta.13) #c4
qqline(cd13_c4$c_delta.13) #c4 

# model C3
cd13_c3_raw_lmer <- lmer((c_delta.13) ~ treatment * site + 
                           (1|plot) + (1|block), data = (cd13_c3)) #c3

# model c4
cd13_c4_raw_lmer <- lmer((c_delta.13) ~ treatment * site + 
                           (1|plot) + (1|block), data = (cd13_c4)) #c4


# residual plots
plot(cd13_c3_raw_lmer, which = 2) #c3
plot(cd13_c4_raw_lmer, which = 2) #c4


# anova
Anova(cd13_c3_raw_lmer) # site ***, treatment:site . #c3
Anova(cd13_c4_raw_lmer) # site ***,  site #c4

plot(cd13_c3_raw_lmer, which = 2)
plot(cd13_c4_raw_lmer, which = 2)


# emmeans C3
cld(emmeans(cd13_c3_raw_lmer, ~site))
cld(emmeans(cd13_c3_raw_lmer, ~treatment))
cld(emmeans(cd13_c3_raw_lmer, ~site*treatment))

pairs(emmeans(cd13_c3_raw_lmer, ~treatment)) # 0.9419
pairs(emmeans(cd13_c3_raw_lmer, ~site))
cld(pairs(emmeans(cd13_c3_raw_lmer, ~site*treatment)))


pairs(emmeans(cd13_c3_raw_lmer, ~treatment, at = list(site = 'sevi'))) #0.31
pairs(emmeans(cd13_c3_raw_lmer, ~treatment, at = list(site = 'temple'))) #0.07 
pairs(emmeans(cd13_c3_raw_lmer, ~treatment, at = list(site = 'arch'))) #0.26
pairs(emmeans(cd13_c3_raw_lmer, ~treatment, at = list(site = 'lubb'))) #0.51



## figure sites raw - c3
z_cd13_c3_raw_sites <- ggplot(cd13_c3, aes(x = treatment, y = c_delta.13, 
                                                 fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("δ13C individual - C3") +
  ylab(expression("δ"^13*"C‰")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  facet_grid( ~ site) +
  coord_cartesian(ylim = c(-33, -25))

 
 ## figure sites raw - c4
z_cd13_c4_raw_sites <- ggplot(cd13_c4, aes(x = treatment, y = c_delta.13, 
                                                 fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("δ13C individual - C4") +
  ylab(expression("δ"^13*"C‰")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  facet_grid( ~ site) +
  coord_cartesian(ylim = c(-19, -11))



## figure all raw - c3
z_cd13_c3_raw_all <- ggplot(cd13_c3, aes(x = treatment, y = c_delta.13, 
                                               fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("δ13C individual - C3") +
  ylab(expression("δ"^13*"C‰")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(1, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(-33, -25))


## figure all raw - c4  
z_cd13_c4_raw_all <- ggplot(cd13_c4, aes(x = treatment, y = c_delta.13, 
                                             fill = treatment)) +
    geom_boxplot(color = "black", outlier.color = "red", outlier.shape = NA, 
                 outlier.fill = "white") +
    ggtitle("δ13C individual - C4") +
  ylab(expression("δ"^13*"C‰")) +
  xlab("") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.0, size = 20),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.text = element_text(size = 13),
      strip.text = element_text(size = 24),  # Increase size of facet titles
      legend.position = "none",  # Remove the legend
      panel.spacing = unit(1, "lines"),  # Increase space between the plots
      panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
    ) +
    scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(-17, -11))




## δ13C, (c_delta.13) cwm ------------------------------------------------------ 

## average per site, plot, & species - c3
cd13_c3_cwm <- cd13_c3 %>%
  group_by(site, plot, taxon_code) %>%
  summarise(cd13_avg = mean(c_delta.13, na.rm = TRUE))

## average per site, plot, & species - c4
cd13_c4_cwm <- cd13_c4 %>%
  group_by(site, plot, taxon_code) %>%
  summarise(cd13_avg = mean(c_delta.13, na.rm = TRUE))

## merging back in plot information and the cover percentages - c3
cd13_c3_cwm <- left_join(cd13_c3_cwm, metadata_plot_treatment)
cd13_c3_cwm <- left_join(cd13_c3_cwm, cwm_spcomp_averages, 
                          by = c("site","plot","taxon_code"))

## merging back in plot information and the cover percentages - c4
cd13_c4_cwm <- left_join(cd13_c4_cwm, metadata_plot_treatment)
cd13_c4_cwm <- left_join(cd13_c4_cwm, cwm_spcomp_averages, 
                         by = c("site","plot","taxon_code"))

## cwm calculation ----------- 

# c3
cd13_c3_cwm_calculated <- cd13_c3_cwm %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cd13_c3_cwm = weighted.mean(cd13_avg, average.cover))
# removing NA values, which is sevi 8(C) & 15(NPK) - c3
cd13_c3_cwm_calculated <- subset(cd13_c3_cwm_calculated, !is.na(cd13_c3_cwm))


## c4
cd13_c4_cwm_calculated <- cd13_c4_cwm %>%
  group_by(site, plot, treatment, block) %>%
  summarise(cd13_c4_cwm = weighted.mean(cd13_avg, average.cover))
## removing NA value, lubb 28(C)
cd13_c4_cwm_calculated <- subset(cd13_c4_cwm_calculated, !is.na(cd13_c4_cwm))

## visualizing c3
hist(cd13_c3_cwm_calculated$cd13_c3_cwm)
qqnorm(cd13_c3_cwm_calculated$cd13_c3_cwm)
qqline(cd13_c3_cwm_calculated$cd13_c3_cwm)

## visualizing c4
hist(cd13_c4_cwm_calculated$cd13_c4_cwm)
qqnorm(cd13_c4_cwm_calculated$cd13_c4_cwm)
qqline(cd13_c4_cwm_calculated$cd13_c4_cwm)


unique(cd13_c3_cwm_calculated$cd13_c3_cwm)

## model c3
cd13_c3_cwm_lmer <- lmer((cd13_c3_cwm) ~ treatment * site +
                           (1|plot) + (1|block), 
                         data = (cd13_c3_cwm_calculated))

## model c4
cd13_c4_cwm_lmer <- lmer((cd13_c4_cwm) ~ treatment * site +
                           (1|plot) + (1|block), 
                         data = (cd13_c4_cwm_calculated))
                          
# residual plots
plot(cd13_c3_cwm_lmer, which = 2) #c3
plot(cd13_c4_cwm_lmer, which = 2) #c4


# anova
Anova(cd13_c3_cwm_lmer) # site ***, treatment:site * #c3
Anova(cd13_c4_cwm_lmer) # site ***, c4

plot(cd13_c3_cwm_lmer, which = 2)
plot(cd13_c4_cwm_lmer, which = 2)


# emmeans 
cld(emmeans(cd13_c3_cwm_lmer, ~site))
cld(emmeans(cd13_c3_cwm_lmer, ~treatment))
cld(emmeans(cd13_c3_cwm_lmer, ~treatment*site)) 

pairs(emmeans(cd13_c3_cwm_lmer, ~treatment, at = list(site='temple'))) #0.14
pairs(emmeans(cd13_c3_cwm_lmer, ~treatment, at = list(site='sevi'))) #0.79
pairs(emmeans(cd13_c3_cwm_lmer, ~treatment, at = list(site='arch'))) #0.46
pairs(emmeans(cd13_c3_cwm_lmer, ~treatment, at = list(site='lubb'))) #0.95
pairs((emmeans(cd13_c3_cwm_lmer, ~treatment))) #0.77
cld(pairs(emmeans(cd13_c3_cwm_lmer, ~treatment*site)))


## showing site differences
# want plots to be listed in order from "wet" to dry.
cd13_c3_site_diffs <- cd13_c3_cwm_calculated
cd13_c3_site_diffs$site <- factor(cd13_c3_site_diffs$site, 
                                  levels = c("arch", "temple", "lubb", "sevi"))

ggplot(cd13_c3_site_diffs, aes(x = site, y = cd13_c3_cwm, 
                                  fill = treatment)) +
  geom_boxplot() +
  facet_wrap(~ treatment)
  


## figure site cwm - c3 
# sites - c3
z_cd13_c3_cwm_sites <- ggplot(cd13_c3_cwm_calculated, 
                                 aes(x = treatment, y = cd13_c3_cwm, 
                                     fill = treatment)) +
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  ggtitle("δ13C cwm - C3") +
  ylab(expression("δ"^13*"C‰")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  facet_grid( ~ site) +
  coord_cartesian(ylim = c(-33, -25))

# all - c3
z_cd13_c3_cwm_all <- ggplot(cd13_c3_cwm_calculated, 
                              aes(x = treatment, y = cd13_c3_cwm, 
                                  fill = treatment)) +
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  ggtitle("δ13C cwm - C3") +
  ylab(expression("δ"^13*"C‰")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  coord_cartesian(ylim = c(-33, -25))

# sites - c4
z_cd13_c4_cwm_sites <- ggplot(cd13_c4_cwm_calculated, 
                              aes(x = treatment, y = cd13_c4_cwm, 
                                  fill = treatment)) +
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  ggtitle("δ13C cwm - C4") +
  ylab(expression("δ"^13*"C‰")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  facet_grid( ~ site) +
  coord_cartesian(ylim = c(-19, -11))


# all - c4
z_cd13_c4_cwm_all <- ggplot(cd13_c4_cwm_calculated, 
                              aes(x = treatment, y = cd13_c4_cwm, 
                                  fill = treatment)) +
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  ggtitle("δ13C cwm - C4") +
  ylab(expression("δ"^13*"C‰")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  coord_cartesian(ylim = c(-17, -11))



## merging plots together

## c3
fig5_cd13_c3_sites <- z_cd13_c3_raw_sites / z_cd13_c3_cwm_sites
fig7_cd13_c3_all <- z_cd13_c3_raw_all + z_cd13_c3_cwm_all

## c3
fig6_cd13_c4_sites <- z_cd13_c4_raw_sites / z_cd13_c4_cwm_sites
fig8_cd13_c4_all <- z_cd13_c4_raw_all + z_cd13_c4_cwm_all


## saving the images
# sites - c3
# png('../figures/fig5_cd13_c3_sites.png',
#     width = 18, height = 8, units = 'in', res = 1000)
# fig5_cd13_c3_sites
# dev.off()
# 
# # sites - c4
# png('../figures/fig6_cd13_c4_sites.png',
#     width = 18, height = 8, units = 'in', res = 1000)
# fig6_cd13_c4_sites
# dev.off()
# 
# 
# # all - c3
# png('../figures/fig7_cd13_c3_all.png',
#     width = 15, height = 7, units = 'in', res = 1000)
# fig7_cd13_c3_all
# dev.off()
# 
# # all - c4
# png('../figures/fig8_cd13_c4_all.png',
#     width = 15, height = 7, units = 'in', res = 1000)
# fig8_cd13_c4_all
# dev.off()

## save this for the CN analysis -----------------------------------
# removing CN isotope issue sample, extremely high (and likely bad/wrong) value
# also remove NA values
soft_full <- 
  soft_full[soft_full$id_full != "sevi_39_boer_1_491", ]



ggplot(soft_spcomp_full, aes(x = treatment, y = cn_ratio, 
                                               fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", 
               outlier.fill = "red") +
  ggtitle("cn raw sites combined") +
  theme_minimal()

##mk probably remove the NA values, and the really large value??
## actually data looks okay for the ratio, might be okay
write.csv(soft_spcomp_full, "../data/03_rproducts/soft_spcomp_full.csv")



### site counts ----------------------------------------------------------------
unique_taxon_counts <- cwm_spcomp_averages %>%
  group_by(site) %>%
  summarize(unique_taxon_count = n_distinct(taxon_code))

unique_samples <- soft_spcomp_full %>%
  group_by(site) %>%
  summarise(samples = n_distinct(id_unique))
  

