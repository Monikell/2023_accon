################################################################################
## Author: Monika Kelley
## Date: 2024/05/16 to _______
## Purpose: Summer 2023, Acquisitive vs. Conservative Grasslands Communities
## Project code/ name: ACCON (AC-quisitive CON-servative)
################################################################################

### ctrl f notes 
# LEMON = last stopped
# HELP = issue/ not sure

setwd("~/git/2023_accon/scripts")

# libraries --------------------------------------------------------------------
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


## plots
library(purrr)
library(broom.mixed)



## PCA
#install.packages("corrr")
library('corrr')
#install.packages("ggcorrplot")
library(ggcorrplot)
#install.packages("FactoMineR")
library(FactoMineR)
# install.packages("factoextra")
library(factoextra)

# data_v0 ----------------------------------------------------------------------
## data
# soft data = field and lab data collected
soft_raw <- read.csv("../data/02_cleaned/field/accon_field.lab_2023.csv")
spcomp_raw <- read.csv("../data/02_cleaned/field/accon_species.comp_2023.csv")

## metadata
metadata_plots <- read.csv("../data/00_meta/plot-descriptions-02-August-2019.csv")
metadata_species.list <- read.csv("../data/00_meta/species_list.csv")





## data: clean -----------------------------------------------------------------


## r code for figures ("z")
# color pallet for figures
colors <- c("Control" = "#E0FFFF", "NPK" = "#7A8B8B")

## data: soft ------------------------------------------------------------------


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


## species comp ----------------------------------------------------------------
### merging plot information with species comp information
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



# A.) calculations -------------------------------------------------------------
## soft traits - calculating soft plant traits

# traits: leaf thickness, SLA, carbon concentration, nitrogen concentration, 
# C:N ratio, SSD, plant height, δ13C (c_delta.13), δ15N, LDMC 
# Coming: phosphorous and potassium



## A1.) SLA --------------------------------------------------------------------
# specific leaf area
# calculating SLA: area/oven-dry mass
soft_calculations$sla_cm2.g1 <- 
  soft_calculations$leaf_area / soft_calculations$leaf_dry_g

# convert sla cm^2/g to m^2/kg 
# soft_calculations$sla_m2.kg <- soft_calculations$sla_cm2.g / 10


## A2.) LDMC -------------------------------------------------------------------
# leaf dry-matter content
# expressed in units of mg g^-1, measured in grams conversion needed first. 


# calculating LDMC: oven-dry mass / water-saturated fresh mass: units mg.g^-1 
soft_calculations$ldmc <- 
  (soft_calculations$leaf_dry_g / soft_calculations$leaf_wet_g) * 1000


## A3.) SSD --------------------------------------------------------------------
# stem-specific density
## calculating SSD: oven-dry mass / volume of fresh stem

## notes
# 2 methods for calculating stem volume: 1) dimensional, 2) water-displacement
# 01) dimensional method: V = (0.5D)^2 X pi X L
# 01) notes) D = diameter (ddh average), L = length of section of stem
# 02) water-displacement: physically calculated in the lab


### SSD.a dimensional method calculations --------------------------------------

################################################ step 1: addressing hollow stems 

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


################################################### step 2: finalizing stem area

# step 2.a: replacing NA values with 0 for the stem_hollow.area
# purpose: highlights solid vs. hollow stems + good for math in step 2.b
soft_calculations <- soft_calculations %>% 
  mutate(stem_hollow.area = replace_na(stem_hollow.area, 0))


## step 2.b: calculating solid cross-sectional area: whole.area - hollow.area
## last step for dealing with hollow stems! can move on to volume calculations 
soft_calculations$stem_area.solid <- 
  soft_calculations$stem_total.area - soft_calculations$stem_hollow.area


##################################################### step 3: calculating volume
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


### SSB.b water-displacement method calculations ------------------------------- 
## (might be off/ wrong)
# small samples, = needed more reps? Not really sure what doing wrong here. 
# use dimensional values, has full data anyway (mostly)


## step 2: calculating ssd 
# volume physically calculated in lab, no math required other than ssd calc. 
 soft_calculations$ssd_displacement <- 
   soft_calculations$stem_dry_mg / soft_calculations$stem_water.displaced.g


## A4.) CN concentration -------------------------------------------------------
## total of sample (µg)/ sample weight (mg)

# carbon
soft_calculations$c_concentration <- 
  soft_calculations$c_total / soft_calculations$cn_sample.weight

# nitrogen
soft_calculations$n_concentration <- 
  soft_calculations$n_total / soft_calculations$cn_sample.weight


## A5.) C:N --------------------------------------------------------------------
# Carbon & Nitrogen ratio ration (c:n) 
## carbon concentration / nitrogen concentration

soft_calculations$cn_ratio <- 
  soft_calculations$c_concentration / soft_calculations$n_concentration


# data_v1 ----------------------------------------------------------------------
#soft plant data calculated creating new data frame  
soft_full <- soft_calculations


# merging values for community weighted means (cwm) 

spcomp_avg.cover

soft_spcomp_full <- left_join(soft_full, metadata_species.list, 
                              by = c("site", "taxon_code"))

# creating a new column for merging that is just "site".
spcomp_avg.cover$site <- gsub("\\.us","",spcomp_avg.cover$site_code)

# merging the soft full, with all the sp. info with the average cover info.
soft_spcomp_full <- left_join(soft_spcomp_full, spcomp_avg.cover, 
                              by = c("site", "plot", "taxon_code"))






# B.) models -------------------------------------------------------------------

## B1.) SLA --------------------------------------------------------------------

### SLA raw --------------------------------------------------------------------

# visualize data
hist(soft_spcomp_full$sla_cm2.g1)
hist(log(soft_spcomp_full$sla_cm2.g1))
qqnorm(log(soft_spcomp_full$sla_cm2.g1))
qqline(log(soft_spcomp_full$sla_cm2.g1))
# log needed

# model v0
sla_raw_lmer <- lmer(log(sla_cm2.g1) ~ treatment * site + (1|block : site),
                     data = (soft_spcomp_full))

# model v1
sla_raw_lmerv1 <- lmer(log(sla_cm2.g1) ~ treatment * site + (1|plot) + (1|block),
                     data = (soft_spcomp_full))


# model v2
sla_raw_lmerv2 <- lmer(log(sla_cm2.g1) ~ treatment * site + 
                         (1|site/block) + (1|plot), data = soft_spcomp_full)


Anova(sla_raw_lmer)
Anova(sla_raw_lmerv1)
Anova(sla_raw_lmerv2)



##NGS: need to include species, are likely to have similar sites. taxon code.
## phylogeny constrained models, 
## those measurements were take on the same species.


# residual plots
plot(resid(sla_raw_lmer) ~fitted(sla_raw_lmer))
plot(sla_raw_lmer, which = 2)

# anova
Anova(sla_raw_lmer) # site ***

unique(soft_spcomp_full$sla_cm2.g1)

# emmeans 
emmeans(sla_raw_lmer, ~site)
cld(emmeans(sla_raw_lmer, ~treatment))
pairs(emmeans(sla_raw_lmer, ~site))

## plotting
z_sla_raw_sites <- ggplot(soft_spcomp_full, aes(x = treatment, y = sla_cm2.g1, 
                                            fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "red", outlier.shape = NA, 
               outlier.fill = "red") +
  ggtitle("sla individuals") +
  ylab(expression("sla "*"(cm"^2*"g"^-1*")")) +
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
  coord_cartesian(ylim = c(0, 300))


# # save image
# png('../figures/z_sla_raw_sites.png',
#     width = 9, height = 8, units = 'in', res = 1000)
#    z_sla_raw_sites
#    dev.off()


 
z_sla_raw_all <- ggplot(soft_spcomp_full, aes(x = treatment, y = sla_cm2.g1, 
                                            fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("sla individuals") +
  ylab(expression("sla "*"(cm"^2*"g"^-1*")")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 30),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 25),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  coord_cartesian(ylim = c(0, 300))

# save image
 # png('../figures/z_sla_raw_all.png',
 #     width = 9, height = 8, units = 'in', res = 1000)
 # z_sla_raw_all
 # dev.off()


### SLA cwm --------------------------------------------------------------------

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
  summarise(sla_avg = mean(sla_cm2.g1, na.rm = TRUE))

## merging back in plot information and the cover percentages  
sla_cwm <- left_join(sla_cwm, metadata_plot_treatment)
sla_cwm <- left_join(sla_cwm, cwm_spcomp_averages, 
                     by = c("site","plot","taxon_code"))

## cwm calculation 
sla_cwm_calculated <- sla_cwm %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(sla_cwm = weighted.mean(sla_avg, average.cover))


## 150% cover and add them all up 75% in that plot. Weight the weights. by the 
# total amout of plant cover

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
  ylab(expression("sla "*"(cm"^2*"g"^-1*")")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 30),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 25),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  coord_cartesian(ylim = c(0, 300))




# merging plots together
fig1_sla_sites <- z_sla_raw_sites / z_sla_cwm_sites
fig2_sla_all <- z_sla_raw_all + z_sla_cwm_all

# # saving the image - sites
# png('../figures/fig1_sla_sites.png',
#     width = 18, height = 8, units = 'in', res = 1000)
# fig1_sla_sites
# dev.off()

# saving the image - all
 # png('../figures/fig2_sla_all.png',
 #     width = 15, height = 7, units = 'in', res = 1000)
 # fig2_sla_all
 # dev.off()




## B2.) C:N --------------------------------------------------------------------

### C:N raw --------------------------------------------------------------------

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

# testing to see if the dropping block is fine? 
cn_raw_lmerv2 <- lmer(log(cn_ratio) ~ treatment * site + (1|plot), data = soft_spcomp_full)

Anova(cn_raw_lmerv2)
Anova(cn_raw_lmer)

# model w/ species
cn_raw_species_lmer <- lmer(log(cn_ratio) ~ taxon_code * treatment * site + 
                              (1|plot) + (1|block), data = (soft_spcomp_full))

# residual plots
plot(resid(cn_raw_lmer) ~fitted(cn_raw_lmer))
plot(cn_raw_lmer, which = 2)

# residual plots w/ species
plot(resid(cn_raw_species_lmer) ~fitted(cn_raw_species_lmer))
plot(cn_raw_species_lmer, which = 2)

# anova
Anova(cn_raw_lmer) # site ***, treatment *

# anova with species
Anova(cn_raw_species_lmer) # taxon ***, treatment ***


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
    plot.title = element_text(hjust = 0.0, size = 30),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 25),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(1, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(0, 70))



 
### C:N cwm --------------------------------------------------------------------

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
cn_ratio_cwm_lmer <- lmer(log(cn_ratio_cwm) ~ treatment * site +
                            (1|plot) + (1|block),
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
    plot.title = element_text(hjust = 0.0, size = 30),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 25),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  facet_grid( ~ site) + 
  coord_cartesian(ylim = c(0, 80))


## just looking at lubb
z_lubb_cn_ratio_cwm_sites <- ggplot(data = subset(cn_ratio_cwm_calculated, 
                                                  site == "lubb"), 
                               aes(x = treatment, y = cn_ratio_cwm, 
                                   fill = treatment)) +
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  ggtitle("C:N cwm - lubb site") +
  ylab("C:N ratio") +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 30),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 25),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors)


z_lubb_cn_ratio_raw_site <- ggplot(data = subset(soft_spcomp_full, 
                                                 site == "lubb"), 
                                    aes(x = treatment, y = cn_ratio, 
                                        fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("C:N individual - lubb site") +
  ylab("C:N ratio") +
  xlab ("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 30),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 25),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(1, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors)




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
    plot.title = element_text(hjust = 0.0, size = 30),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 25),
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
#sites
# png('../figures/fig3_cn_ratio_sites.png',
#     width = 18, height = 8, units = 'in', res = 1000)
# fig3_cn_ratio_sites
# dev.off()
# 
# all
# png('../figures/fig4_cn_ratio_all.png',
#     width = 15, height = 7, units = 'in', res = 1000)
# fig4_cn_ratio_all
# dev.off()

Anova(cn_raw_lmer)
Anova(cn_ratio_cwm_lmer)
## B3.) δ13C -------------------------------------------------------------------

### δ13C raw -------------------------------------------------------------------
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
    plot.title = element_text(hjust = 0.0, size = 30),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 25),
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
      plot.title = element_text(hjust = 0.0, size = 30),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.text = element_text(size = 25),
      strip.text = element_text(size = 24),  # Increase size of facet titles
      legend.position = "none",  # Remove the legend
      panel.spacing = unit(1, "lines"),  # Increase space between the plots
      panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
    ) +
    scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(-17, -11))


## combining raw (individual) delta 13 c plots

fig11_cd13_c3.c4_raw_all <- z_cd13_c3_raw_all + z_cd13_c4_raw_all


# saving the plots
# ##sites
# png('../figures/fig11_cd13_c3.c4_raw_all.png',
#     width = 18, height = 8, units = 'in', res = 1000)
# fig11_cd13_c3.c4_raw_all
# dev.off()
# 
# ##all
# png('../figures/fig4_cn_ratio_all.png',
#     width = 15, height = 7, units = 'in', res = 1000)
# fig4_cn_ratio_all
# dev.off()


### δ13C cwm -------------------------------------------------------------------

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

## cwm calculation 

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

## c4
fig6_cd13_c4_sites <- z_cd13_c4_raw_sites / z_cd13_c4_cwm_sites
fig8_cd13_c4_all <- z_cd13_c4_raw_all + z_cd13_c4_cwm_all



# saving the plots
#sites
# png('../figures/fig3_cn_ratio_sites.png',
#     width = 18, height = 8, units = 'in', res = 1000)
# fig3_cn_ratio_sites
# dev.off()
# 
# all
# png('../figures/fig4_cn_ratio_all.png',
#     width = 15, height = 7, units = 'in', res = 1000)
# fig4_cn_ratio_all
# dev.off()


## cwm for both c3 and c4 combined 
## photopathway 

cd13_pathways_cwm <- cwm_soft_spcomp %>%
  group_by(site, plot, taxon_code, photo_pathway) %>%
  summarise(cd13_pathway_averages = mean(c_delta.13, na.rm = TRUE))

## merging back in plot information and the cover percentages  
cd13_pathways_cwm <- left_join(cd13_pathways_cwm, metadata_plot_treatment)
cd13_pathways_cwm <- left_join(cd13_pathways_cwm, cwm_spcomp_averages, 
                          by = c("site","plot","taxon_code"))

## cwm calculation,
cd13_pathways_cwm_calculated <- cd13_pathways_cwm %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(cd13_pathways_cwm = 
              weighted.mean(cd13_pathway_averages, average.cover, na.rm=TRUE))



## visualizing
hist(cd13_pathways_cwm_calculated$cd13_pathways_cwm)
qqnorm(cd13_pathways_cwm_calculated$cd13_pathways_cwm)
qqline(cd13_pathways_cwm_calculated$cd13_pathways_cwm)
# log looks a little better


# model
cd13_pathways_cwm_lmer <- lmer((cd13_pathways_cwm) ~ treatment * 
                            site + (1|plot) + (1|block),
                          data = (cd13_pathways_cwm_calculated))


# residual plots
plot(cd13_pathways_cwm_lmer, which = 2)

# anova
Anova(cd13_pathways_cwm_lmer) 


# order sites wet to dry
cd13_pathways_cwm_calculated$site <- factor(cd13_pathways_cwm_calculated$site, 
                                  levels = c("arch", "temple", "lubb", "sevi"))


# figure 
fig10_cd13_pathways_cwm_all <- ggplot(cd13_pathways_cwm_calculated, 
                              aes(x = treatment, y = cd13_pathways_cwm, 
                                  fill = treatment)) +
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  ggtitle("δ13C cwm - all photopathways") +
  ylab(expression("δ"^13*"C‰")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 30),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 25),
    strip.text = element_text(size = 30),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors)


### saving the plot
# sites
# png('../figures/fig9_cd13_pathways_cwm_sites.png',
#     width = 18, height = 8, units = 'in', res = 1000)
# fig9_cd13_pathways_cwm_sites
# dev.off()

# # all
# png('../figures/fig10_cd13_pathways_cwm_all.png',
#     width = 15, height = 7, units = 'in', res = 1000)
# fig10_cd13_pathways_cwm_all
# dev.off()




## B4.) X (chi) ----------------------------------------------------------------

# want both individuals and cwm

x_calc <- soft_spcomp_full

## calculate big delta C13 from small delta 
x_calc$delta = ((-0.008 - x_calc$c_delta.13 * 0.001) / 
                  (1 + x_calc$c_delta.13 * 0.001)) * 1000

hist(x_calc$delta)

x_calc$photo_pathway

## calculate chi for C3 plants 
x_calc$chi[x_calc$photo_pathway == 'c3'] = 
  (x_calc$delta[x_calc$photo_pathway == 'c3'] * 
     0.001 - 0.0044) / (0.027 - 0.0044)
hist(x_calc$chi)


## calculate chi for C4 plants 
x_calc$chi[x_calc$photo_pathway == 'c4'] = 
  (x_calc$delta[x_calc$photo_pathway == 'c4'] *
     0.001 - 0.0044) / ((-0.0057 + 0.03*0.4) - 0.0044)




### X raw ----------------------------------------------------------------------

# removing NA values, those that did not get a delta 13 value back
# and problem sample #491
x_calc <- subset(x_calc, !is.na(chi))
x_calc <- x_calc[x_calc$id_full != "sevi_39_boer_1_491", ]

hist(x_calc$chi)
qqnorm(x_calc$chi) ; qqline(x_calc$chi) 
# model
x_raw_lmer <- lmer((chi) ~ treatment * site + (1|plot) + (1|block),
                    data = (x_calc))

# residual plots
plot(resid(x_raw_lmer) ~fitted(x_raw_lmer))
plot(x_raw_lmer, which = 2)

# anova
Anova(x_raw_lmer) # site ***

# emmeans 
emmeans(x_raw_lmer, ~site)
emmeans(x_raw_lmer, ~treatment)
emmeans(x_raw_lmer, ~site*treatment)
cld(emmeans(x_raw_lmer, ~treatment))
cld(emmeans(x_raw_lmer, ~site))

pairs(emmeans(x_raw_lmer, ~treatment))
pairs(emmeans(x_raw_lmer, ~treatment, at = list(site = 'temple')))
pairs(emmeans(x_raw_lmer, ~treatment, at = list(site = 'arch')))
pairs(emmeans(x_raw_lmer, ~treatment, at = list(site = 'lubb')))
pairs(emmeans(x_raw_lmer, ~treatment, at = list(site = 'sevi')))

# order sites wet to dry
x_calc$site <- factor(x_calc$site, 
                                levels = c("arch", "temple", "lubb", "sevi"))


# figure sites
z_x_c3_raw_all <- ggplot(data = subset(x_calc, photo_pathway == "c3"), 
                               aes(x = treatment, y = chi, 
                                   fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("χ individuals - C3") +
  ylab("χ") +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 30),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 25),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors)


# z_x_c4_raw_all <- ggplot(data = subset(x_c4_raw_correct_values, 
#                                          photo_pathway == "c4" &
#                                            chi >= 0.05 & 
#                                            chi <= 0.95), 
#                            aes(x = treatment, y = chi, 
#                                fill = treatment)) +
#   geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
#                outlier.fill = "white") +
#   ggtitle("χ individuals - C4") +
#   ylab("χ") +
#   xlab("") +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.0, size = 30),
#     axis.title.x = element_text(size = 15),
#     axis.title.y = element_text(size = 15),
#     axis.text = element_text(size = 25),
#     strip.text = element_text(size = 24),  # Increase size of facet titles
#     legend.position = "none",  # Remove the legend
#     panel.spacing = unit(4, "lines"),  # Increase space between the plots
#     panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
#   ) +
#   scale_fill_manual(values = colors)


## merging plots together

# fig12_x_all <- z_x_c3_raw_all + z_x_c4_raw_all

## saving plot figure
# png('../figures/fig12_x_all.png',
#     width = 18, height = 8, units = 'in', res = 1000)
# fig12_x_all
# dev.off()

### X cwm ---------------------------------------------------------------------- 


## average cn_ratio, per site, plot, & species
x_cwm <- x_calc %>%
  group_by(site, plot, taxon_code) %>%
  summarise(chi_avg = mean(chi, na.rm = TRUE))

## merging back in plot information and the cover percentages  
x_cwm <- left_join(x_cwm, metadata_plot_treatment)
x_cwm <- left_join(x_cwm, cwm_spcomp_averages, 
                          by = c("site","plot","taxon_code"))

## cwm calculation,
x_cwm_calculated <- x_cwm %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(x_cwm = 
              weighted.mean(chi_avg, average.cover, na.rm=TRUE))

## visualizing
hist(x_cwm_calculated$x_cwm) 
hist(log(x_cwm_calculated$x_cwm))
qqnorm(x_cwm_calculated$x_cwm) ; qqline(x_cwm_calculated$x_cwm)
qqnorm(log(x_cwm_calculated$x_cwm))
qqline(log(x_cwm_calculated$x_cwm))
# log looks a little better??

# model, nested lemer 
x_cwm_lmer <- lmer(log(x_cwm) ~ treatment *
                            site + (1|block : site),
                          data = (x_cwm_calculated))


# residual plots
plot(x_cwm_lmer, which = 2)

# anova
Anova(x_cwm_lmer) # site ***

# emmeans 
cld(emmeans(x_cwm_lmer, ~site))
cld(emmeans(x_cwm_lmer, ~treatment))
cld(emmeans(x_cwm_lmer, ~treatment*site))

# order sites wet to dry
x_cwm_calculated$site <- factor(x_cwm_calculated$site, 
                                  levels = c("arch", "temple", "lubb", "sevi"))

# plotting - sites
z_x_cwm_sites <- ggplot(x_cwm_calculated, 
                               aes(x = treatment, y = x_cwm, 
                                   fill = treatment)) +
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  ggtitle("X (chi) cwm") +
  ylab("X(chi)") +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 30),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 25),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  facet_grid( ~ site) +
  coord_cartesian(ylim = c(-0.5, 3))




x_cwm_calculated$x_cwm



# plotting - all
fig13_x_cwm_all <- ggplot(data = subset(x_cwm_calculated, 
                                    x_cwm >= 0.05 &
                                      x_cwm <= 0.95),
                      aes(x = treatment, y = x_cwm, 
                          fill = treatment)) +
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  ggtitle("χ cwm - all photopathways") +
  ylab("χ") +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 30),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 25),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(1, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors)


# combining plots
# fig10_x_sites <- z_x_raw_sites / z_x_cwm_sites 
# fig11_x_all <- z_x_raw_all + z_x_cwm_all


## saving the plots
# sites
# png('../figures/fig10_x_sites.png',
#     width = 18, height = 8, units = 'in', res = 1000)
# fig10_x_sites
# dev.off()

# # all
# png('../figures/fig13_x_cwm_all.png',
#     width = 15, height = 7, units = 'in', res = 1000)
# fig13_x_cwm_all
# dev.off()


# site counts ------------------------------------------------------------------
unique_taxon_counts <- cwm_spcomp_averages %>%
  group_by(site) %>%
  summarize(unique_taxon_count = n_distinct(taxon_code))


# Your data
unique_taxon_counts <- tibble::tibble(
  site = c("arch", "lubb", "sevi", "temple"),
  unique_taxon_count = c(11, 6, 7, 23)
)

# Assign factor levels from driest to wettest
unique_taxon_counts <- unique_taxon_counts %>%
  mutate(site = factor(site, levels = c("sevi", "lubb", "temple", "arch")))

# Plot
ggplot(unique_taxon_counts, aes(x = site, y = unique_taxon_count, fill = site)) +
  geom_col() +
  scale_fill_manual(values = c("sevi" = "#E69F00",   # driest
                               "lubb" = "#F0E442",
                               "temple" = "#56B4E9",
                               "arch" = "#0072B2")) +  # wettest
  labs(x = "Site", y = "Unique Taxa Count",
       title = "Unique Taxa per Site",
       fill = "Site (Dry to Wet)") +
  theme_minimal()
11 + 6 + 7 + 23

unique_samples <- soft_spcomp_full %>%
  group_by(site) %>%
  summarise(samples = n_distinct(id_unique))
  
soft_spcomp_full$photo_pathway

c3_c4_counts <- soft_spcomp_full %>%
  group_by(photo_pathway) %>%
  summarize(unique_taxon = n_distinct(taxon_code))

c3_c4_counts_site <- soft_spcomp_full %>%
  group_by(photo_pathway, site) %>%
  summarize(unique_taxon = n_distinct(taxon_code))




# 07/11/2025 updates -----------------------------------------------------------

## B5.) plant height -----------------------------------------------------------

### plant height raw -----------------------------------------------------------

soft_spcomp_full$plant_height

# visualize data
hist(soft_spcomp_full$plant_height)
hist(log(soft_spcomp_full$plant_height)) # logged
qqline(log(soft_spcomp_full$plant_height))
# log needed

# model
pheight_raw_lmer <- lmer(log(plant_height) ~ treatment * site + 
                           (1|plot) + (1|block),
                     data = (soft_spcomp_full))



# residual plots
plot(resid(pheight_raw_lmer) ~fitted(pheight_raw_lmer))
plot(pheight_raw_lmer, which = 2)

# anova
Anova(pheight_raw_lmer) # site ***, treatment * (0.04553)

unique(soft_spcomp_full$plant_height)

# emmeans 
emmeans(pheight_raw_lmer, ~site)
cld(emmeans(pheight_raw_lmer, ~treatment))
pairs(emmeans(pheight_raw_lmer, ~site))

## plotting
z_pheight_raw_sites <- ggplot(soft_spcomp_full, aes(x = treatment, y = plant_height, 
                                                fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "red", outlier.shape = NA, 
               outlier.fill = "red") +
  ggtitle("plant height individuals") +
  ylab(expression("plant height (cm)")) +
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
  coord_cartesian(ylim = c(0, 300))


# # save image
# png('../figures/z_sla_raw_sites.png',
#     width = 9, height = 8, units = 'in', res = 1000)
#    z_sla_raw_sites
#    dev.off()



z_pheight_raw_all <- ggplot(soft_spcomp_full, aes(x = treatment, y = plant_height, 
                                              fill = treatment)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = NA, 
               outlier.fill = "white") +
  ggtitle("plant height individuals") +
  ylab(expression("plant height (cm)")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 30),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 25),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  coord_cartesian(ylim = c(0, 300))

# save image
# png('../figures/z_sla_raw_all.png',
#     width = 9, height = 8, units = 'in', res = 1000)
# z_sla_raw_all
# dev.off()


### plant height cwm -----------------------------------------------------------

## average plant height
pheight_cwm <- cwm_soft_spcomp %>%
  group_by(site, plot, taxon_code) %>%
  summarise(pheight_avg = mean(plant_height, na.rm = TRUE))

## merging back in plot information and the cover percentages  
pheight_cwm <- left_join(pheight_cwm, metadata_plot_treatment)
pheight_cwm <- left_join(pheight_cwm, cwm_spcomp_averages, 
                     by = c("site","plot","taxon_code"))

## cwm calculation 
pheight_cwm_calculated <- pheight_cwm %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(pheight_cwm = weighted.mean(pheight_avg, average.cover))

unique(pheight_cwm_calculated$pheight_cwm)

## visualizing
hist(log(pheight_cwm_calculated$pheight_cwm))
qqnorm(log(pheight_cwm_calculated$pheight_cwm))
qqline(log(pheight_cwm_calculated$pheight_cwm))
qqnorm(pheight_cwm_calculated$pheight_cwm)
qqline(pheight_cwm_calculated$pheight_cwm)
# log looks a little better

# model
pheight_cwm_lmer <- lmer(log(pheight_cwm) ~ treatment * site +
                           (1|plot) + (1|block),
                     data = (pheight_cwm_calculated))



# residual plots
plot(resid(pheight_cwm_lmer) ~fitted(pheight_cwm_lmer))
plot(pheight_cwm_lmer, which = 2)

# anova
Anova(pheight_cwm_lmer) # site ***, treatment .
Anova(pheight_raw_lmer) # site ***, treatment *

emmeans(pheight_cwm_lmer, ~site*treatment)



# emmeans 
emmeans(pheight_cwm_lmer, ~site)
cld(emmeans(pheight_cwm_lmer, ~treatment))

# plotting
z_pheight_cwm_sites <- ggplot(pheight_cwm_calculated, 
                          aes(x = treatment, y = pheight_cwm, fill = treatment)) + 
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  ggtitle("plant height cwm") +
  ylab(expression("plant height (cm)")) +
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
  facet_grid(~ site)



## all sites
z_pheight_cwm_all <- ggplot(pheight_cwm_calculated, aes(x = treatment, y = pheight_cwm, 
                                                fill = treatment)) +
  geom_boxplot_pattern(color = "black", outlier.color = "black", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  ggtitle("plant height cwm") +
  ylab(expression("plant height (cm)")) +
  xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 30),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 25),
    strip.text = element_text(size = 24),  # Increase size of facet titles
    legend.position = "none",  # Remove the legend
    panel.spacing = unit(4, "lines"),  # Increase space between the plots
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  scale_fill_manual(values = colors) + 
  coord_cartesian(ylim = c(0, 300))


## B6.) leaf thickness ---------------------------------------------------------

### leaf thickness raw ---------------------------------------------------------

unique(soft_spcomp_full$leaf_thickness)

# visualize data
hist(soft_spcomp_full$leaf_thickness)
hist(log(soft_spcomp_full$leaf_thickness)) # logged
qqnorm(log(soft_spcomp_full$leaf_thickness))
qqline(log(soft_spcomp_full$leaf_thickness))
# log needed

# model
leaf_thickness_raw_lmer <- lmer(log(leaf_thickness) ~ treatment * site + 
                           (1|plot) + (1|block),
                         data = (soft_spcomp_full))



# residual plots
plot(resid(leaf_thickness_raw_lmer) ~fitted(leaf_thickness_raw_lmer))
plot(leaf_thickness_raw_lmer, which = 2)

# anova
Anova(leaf_thickness_raw_lmer) # site ***

unique(soft_spcomp_full$leaf_thickness)




## plotting sites
ggplot(soft_spcomp_full, aes(x = treatment, y = leaf_thickness, 
                             fill = treatment)) + 
  geom_boxplot() +
  theme_minimal() +   
  facet_grid( ~ site) +
  coord_cartesian(ylim = c(0, 2))


## plotting all
ggplot(soft_spcomp_full, aes(x = treatment, y = leaf_thickness, 
                             fill = treatment)) + 
  geom_boxplot() +
  theme_minimal() + 
  coord_cartesian(ylim = c(0, 2))


### leaf thickness cwm ---------------------------------------------------------

## average leaf thickness
leaf_thickness_cwm <- cwm_soft_spcomp %>%
  group_by(site, plot, taxon_code) %>%
  summarise(leaf_thickness_avg = mean(leaf_thickness, na.rm = TRUE))

## merging back in plot information and the cover percentages  
leaf_thickness_cwm <- left_join(leaf_thickness_cwm, metadata_plot_treatment)
leaf_thickness_cwm <- left_join(leaf_thickness_cwm, cwm_spcomp_averages, 
                         by = c("site","plot","taxon_code"))

## cwm calculation 
leaf_thickness_cwm_calculated <- leaf_thickness_cwm %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(leaf_thickness_cwm = weighted.mean(leaf_thickness_avg, average.cover))

unique(leaf_thickness_cwm_calculated$leaf_thickness_cwm)

## visualizing
hist(leaf_thickness_cwm_calculated$leaf_thickness_cwm)
hist(log(leaf_thickness_cwm_calculated$leaf_thickness_cwm))
qqnorm(log(leaf_thickness_cwm_calculated$leaf_thickness_cwm))
qqline(log(leaf_thickness_cwm_calculated$leaf_thickness_cwm))
qqnorm(leaf_thickness_cwm_calculated$leaf_thickness_cwm)
qqline(leaf_thickness_cwm_calculated$leaf_thickness_cwm)
# log looks a little better

# model
leaf_thickness_cwm_lmer <- lmer(log(leaf_thickness_cwm) ~ treatment * site +
                           (1|plot) + (1|block),
                         data = (leaf_thickness_cwm_calculated))



# residual plots
plot(resid(leaf_thickness_cwm_lmer) ~fitted(leaf_thickness_cwm_lmer))
plot(leaf_thickness_cwm_lmer, which = 2)

# anova
Anova(leaf_thickness_cwm_lmer) # site ***, trt:site **
Anova(leaf_thickness_raw_lmer) # site ***

emmeans(leaf_thickness_cwm_lmer, ~site*treatment)


# emmeans 
emmeans(leaf_thickness_cwm_lmer, ~site)
cld(emmeans(leaf_thickness_cwm_lmer, ~treatment))


# plotting
ggplot(leaf_thickness_cwm_calculated, 
       aes(x = treatment, 
           y = leaf_thickness_cwm, 
           fill = treatment)) + 
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  theme_minimal() +
  scale_fill_manual(values = colors) +
  facet_grid(~ site)



## all sites
ggplot(leaf_thickness_cwm_calculated, 
       aes(x = treatment, 
           y = leaf_thickness_cwm, 
           fill = treatment)) + 
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  theme_minimal() +
  scale_fill_manual(values = colors)


## B7.) LDMC -------------------------------------------------------------------

### LDMC raw -------------------------------------------------------------------

unique(soft_spcomp_full$ldmc)
mean(soft_spcomp_full$ldmc)
mean(soft_spcomp_full$ldmc)

# removing na values for "ldmc"
soft_spcomp_ldmc_cleaned <- soft_spcomp_full[!is.na(soft_spcomp_full$ldmc), ]
unique(soft_spcomp_ldmc_cleaned$ldmc)
mean(soft_spcomp_ldmc_cleaned$ldmc)
min(soft_spcomp_ldmc_cleaned$ldmc)
max(soft_spcomp_ldmc_cleaned$ldmc)



# visualize data
hist(soft_spcomp_ldmc_cleaned$ldmc)
hist(log(soft_spcomp_ldmc_cleaned$ldmc)) # logged
qqnorm(log(soft_spcomp_ldmc_cleaned$ldmc))
qqline(log(soft_spcomp_ldmc_cleaned$ldmc))
# log needed

# model
ldmc_raw_lmer <- lmer(log(ldmc) ~ treatment * site + 
                                  (1|plot) + (1|block),
                                data = (soft_spcomp_ldmc_cleaned))



# residual plots
plot(resid(ldmc_raw_lmer) ~fitted(ldmc_raw_lmer))
plot(ldmc_raw_lmer, which = 2)

# anova
Anova(ldmc_raw_lmer) # site ***



## plotting sites
ggplot(soft_spcomp_ldmc_cleaned, aes(x = treatment, y = ldmc, 
                             fill = treatment)) + 
  geom_boxplot() +
  theme_minimal() +   
  facet_grid( ~ site) +
  coord_cartesian(ylim = c(0, 1000))


## plotting all
ggplot(soft_spcomp_ldmc_cleaned, aes(x = treatment, y = ldmc, 
                             fill = treatment)) + 
  geom_boxplot() +
  theme_minimal() + 
  coord_cartesian(ylim = c(0, 1000))


### LDMC cwm ---------------------------------------------------------

## average leaf thickness
ldmc_cwm <- cwm_soft_spcomp %>%
  group_by(site, plot, taxon_code) %>%
  summarise(ldmc_avg = mean(ldmc, na.rm = TRUE))

## merging back in plot information and the cover percentages  
ldmc_cwm <- left_join(ldmc_cwm, metadata_plot_treatment)
ldmc_cwm <- left_join(ldmc_cwm, cwm_spcomp_averages, 
                                by = c("site","plot","taxon_code"))


## cwm calculation 
ldmc_cwm_calculated <- ldmc_cwm %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(ldmc_cwm = weighted.mean(ldmc_avg, average.cover))

unique(ldmc_cwm_calculated$ldmc_cwm)

## visualizing
hist(ldmc_cwm_calculated$ldmc_cwm)
hist(log(ldmc_cwm_calculated$ldmc_cwm))
qqnorm(log(ldmc_cwm_calculated$ldmc_cwm))
qqline(log(ldmc_cwm_calculated$ldmc_cwm))

# log looks a little better

# model
ldmc_cwm_lmer <- lmer(log(ldmc_cwm) ~ treatment * site +
                                  (1|plot) + (1|block),
                                data = (ldmc_cwm_calculated))



# residual plots
plot(resid(ldmc_cwm_lmer) ~fitted(ldmc_cwm_lmer))
plot(ldmc_cwm_lmer, which = 2)

# anova
Anova(ldmc_cwm_lmer) # site ***
Anova(ldmc_cwm_lmer) # site ***

emmeans(ldmc_cwm_lmer, ~site)


# emmeans 
emmeans(ldmc_cwm_lmer, ~site)
cld(emmeans(ldmc_cwm_lmer, ~site))


# plotting
ggplot(ldmc_cwm_calculated, 
       aes(x = treatment, 
           y = ldmc_cwm, 
           fill = treatment)) + 
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  theme_minimal() +
  scale_fill_manual(values = colors) +
  facet_grid(~ site)



## all sites
ggplot(ldmc_cwm_calculated, 
       aes(x = treatment, 
           y = ldmc_cwm, 
           fill = treatment)) + 
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  theme_minimal() +
  scale_fill_manual(values = colors)



## B8.) SSD --------------------------------------------------------------------

### SSD raw --------------------------------------------------------------------

unique(soft_spcomp_full$ssd_dimensional)


# removing na values for "ldmc"
soft_spcomp_ssd_cleaned <- soft_spcomp_full[!is.na(soft_spcomp_full$ssd_dimensional), ]
unique(soft_spcomp_ssd_cleaned$ssd_dimensional)
mean(soft_spcomp_ssd_cleaned$ssd_dimensional)
min(soft_spcomp_ssd_cleaned$ssd_dimensional)
max(soft_spcomp_ssd_cleaned$ssd_dimensional)



# visualize data
hist(soft_spcomp_ssd_cleaned$ssd_dimensional)
hist(log(soft_spcomp_ssd_cleaned$ssd_dimensional))
qqnorm(log(soft_spcomp_ssd_cleaned$ssd_dimensional))
qqline(log(soft_spcomp_ssd_cleaned$ssd_dimensional))
# log needed

# model
ssd_raw_lmer <- lmer(log(ssd_dimensional) ~ treatment * site + 
                        (1|plot) + (1|block),
                      data = (soft_spcomp_ssd_cleaned))



# residual plots
plot(resid(ssd_raw_lmer) ~fitted(ssd_raw_lmer))
plot(ssd_raw_lmer, which = 2)

# anova
Anova(ssd_raw_lmer) # site **, treatment 0.09



## plotting sites
ggplot(soft_spcomp_ssd_cleaned, aes(x = treatment, y = ssd_dimensional, 
                                     fill = treatment)) + 
  geom_boxplot() +
  theme_minimal() +   
  facet_grid( ~ site) +
  coord_cartesian(ylim = c(0, 1.5))


## plotting all
ggplot(soft_spcomp_ssd_cleaned, aes(x = treatment, y = ssd_dimensional, 
                                     fill = treatment)) + 
  geom_boxplot() +
  theme_minimal() + 
  coord_cartesian(ylim = c(0, 1.5))


### SSD cwm ---------------------------------------------------------

## averages
ssd_cwm <- cwm_soft_spcomp %>%
  group_by(site, plot, taxon_code) %>%
  summarise(ssd_avg = mean(ssd_dimensional, na.rm = TRUE))

## merging back in plot inf## merging back in plot inf## merging back in plot information and the cover percentages  
ssd_cwm <- left_join(ssd_cwm, metadata_plot_treatment)
ssd_cwm <- left_join(ssd_cwm, cwm_spcomp_averages, 
                      by = c("site","plot","taxon_code"))


## cwm calculation 
ssd_cwm_calculated <- ssd_cwm %>%
  group_by(site, plot, treatment, block) %>% 
  summarise(ssd_cwm = weighted.mean(ssd_avg, average.cover))

unique(ssd_cwm_calculated$ssd_cwm)

## visualizing
hist(ssd_cwm_calculated$ssd_cwm)
hist(log(ssd_cwm_calculated$ssd_cwm))
qqnorm(log(ssd_cwm_calculated$ssd_cwm))
qqline(log(ssd_cwm_calculated$ssd_cwm))

# log looks a little better

# model
ssd_cwm_lmer <- lmer(log(ssd_cwm) ~ treatment * site +
                        (1|plot) + (1|block),
                      data = (ssd_cwm_calculated))



# residual plots
plot(resid(ssd_cwm_lmer) ~fitted(ssd_cwm_lmer))
plot(ssd_cwm_lmer, which = 2)

# anova
Anova(ssd_cwm_lmer) # site *


emmeans(ssd_cwm_lmer, ~site)


# emmeans 
emmeans(ssd_cwm_lmer, ~site)
cld(emmeans(ssd_cwm_lmer, ~site))


# plotting
ggplot(ssd_cwm_calculated, 
       aes(x = treatment, 
           y = ssd_cwm, 
           fill = treatment)) + 
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  theme_minimal() +
  scale_fill_manual(values = colors) +
  facet_grid(~ site) + 
  coord_cartesian(ylim = c(0.20, 1))



## all sites
ggplot(ssd_cwm_calculated, 
       aes(x = treatment, 
           y = ssd_cwm, 
           fill = treatment)) + 
  geom_boxplot_pattern(color = "black", outlier.color = "red", 
                       outlier.shape = NA,
                       pattern = "stripe",  # same pattern for each 
                       pattern_fill = "black", pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_angle = 45) +
  theme_minimal() +
  scale_fill_manual(values = colors)





## Forest plots? ---------------------------------------------------------------


# Traits to model
traits <- c("plant_height", "leaf_thickness", "ldmc", "ssd", "sla_cm2.g1", 
            "ssd_dimensional", "cn_ratio")

# Filter out traits with all NA values or insufficient data
traits <- traits[traits %in% names(soft_spcomp_full)]  # Keep only traits that exist
traits <- traits[sapply(soft_spcomp_full[traits], function(x) sum(!is.na(x)) > 10)]  # More than 10 non-missing

# Fit models and extract treatment effects
model_results <- traits %>%
  set_names() %>%
  map_df(~ {
    trait_data <- soft_spcomp_full %>%
      filter(!is.na(.data[[.x]])) %>%
      filter(if (.x == "ldmc") .data[[.x]] >= 100 & .data[[.x]] < 600 else TRUE)  # Exclude extreme ldmc values
    
    mod <- tryCatch(
      lmer(as.formula(paste0(.x, " ~ treatment * site + (1|plot) + (1|block)")),
           data = trait_data),
      error = function(e) NULL
    )
    
    if (!is.null(mod)) {
      tidy(mod, effects = "fixed", conf.int = TRUE) %>%
        filter(term == "treatmentNPK") %>%
        mutate(trait = .x)
    }
  }, .id = "trait") %>%
  filter(!is.na(estimate))

# Identify significant results based on confidence intervals
model_results <- model_results %>%
  mutate(significant = !(conf.low <= 0 & conf.high >= 0))  # CI does not overlap zero

# Forest plot
ggplot(model_results, aes(x = estimate, y = reorder(trait, estimate), color = significant)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("gray30", "red3")) +
  labs(
    title = "Effect of NPK Treatment on Plant Traits",
    x = "Estimate (NPK vs Control)",
    y = "Trait",
    color = "Significant\n(p < 0.05)"
  ) +
  theme_minimal(base_size = 13)



### now the community weighted means
cwm_merged <- leaf_thickness_cwm_calculated %>%
  full_join(ldmc_cwm_calculated) %>%
  full_join(sla_cwm_calculated) %>%
  full_join(cn_ratio_cwm_calculated) %>%
  full_join(ssd_cwm_calculated) %>%
  full_join(pheight_cwm_calculated)
cwm_merged  

colnames(cwm_merged)

traits <- c("leaf_thickness_cwm", "ldmc_cwm", "sla_cwm", "cn_ratio_cwm",
            "ssd_cwm", "pheight_cwm")

model_results_cwm <- traits %>%
  set_names() %>%
  map_df(~ {
    trait_data <- cwm_merged %>%
      filter(!is.na(.data[[.x]])) %>%
      filter(if (.x == "ldmc") .data[[.x]] >= 150 & .data[[.x]] < 600 else TRUE)
    
    mod <- tryCatch(
      lmer(as.formula(paste0(.x, " ~ treatment * site + (1|plot) + (1|block)")),
           data = trait_data),
      error = function(e) NULL
    )
    
    if (!is.null(mod)) {
      tidy(mod, effects = "fixed", conf.int = TRUE) %>%
        filter(term == "treatmentNPK") %>%
        mutate(trait = .x)
    }
  }, .id = "trait") %>%
  filter(!is.na(estimate)) %>%
  mutate(significant = !(conf.low <= 0 & conf.high >= 0))




## forest
ggplot(model_results_cwm, aes(x = estimate, y = reorder(trait, estimate), color = significant)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("gray30", "red3")) +
  labs(
    title = "Effect of NPK Treatment on Community-Weighted Means",
    x = "Estimate (NPK vs Control)",
    y = "Trait",
    color = "Significant\n(p < 0.05)"
  ) +
  theme_minimal(base_size = 13)


### new attempt
trait_labels <- c(
  plant_height = "Plant Height",
  leaf_thickness = "Leaf Thickness",
  ldmc = "LDMC",
  ssd = "Stem Specific Density",
  sla_cm2.g1 = "SLA",
  ssd_dimensional = "SSD",
  cn_ratio = "Leaf C:N"
)

model_results <- model_results %>%
  mutate(trait = trait_labels[trait])

# Updated forest plot
ggplot(model_results, aes(x = estimate, y = reorder(trait, estimate), color = significant)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("gray30", "red3")) +
  labs(
    title = "Effect of NPK Treatment on Plant Traits",
    x = "Estimate (NPK vs Control)",
    y = "Trait",
    color = "Significant\n(CI excludes 0)"
  ) +
  theme_minimal(base_size = 13)







z_forest <- ggplot(model_results, aes(x = estimate, y = reorder(trait, estimate), color = significant)) +
  geom_point(size = 5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("gray30", "red3")) +
  labs(
    title = "Effect of NPK Treatment on Plant Traits",
    x = "Estimate (NPK vs Control)",
    y = NULL,
    color = "Significant\n(CI excludes 0)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.0, size = 20),
    axis.title.x = element_text(size = 30, color = "#535353"),
    axis.title.y = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )



# saving the plot

png('../figures/fig_forest.png',
    width = 17, height = 8, units = 'in', res = 1000)
z_forest
dev.off()


### try to add the delta 13C stuff...



# model C3
cd13_c3_raw_lmer <- lmer((c_delta.13) ~ treatment * site + 
                           (1|plot) + (1|block), data = (cd13_c3)) #c3

# model c4
cd13_c4_raw_lmer <- lmer((c_delta.13) ~ treatment * site + 
                           (1|plot) + (1|block), data = (cd13_c4)) #c4


c3_result <- tidy(cd13_c3_raw_lmer, effects = "fixed", conf.int = TRUE) %>%
  filter(term == "treatmentNPK") %>%
  mutate(trait = "δ13C (C3)", source = "Individual")


c4_result <- tidy(cd13_c4_raw_lmer, effects = "fixed", conf.int = TRUE) %>%
  filter(term == "treatmentNPK") %>%
  mutate(trait = "δ13C (C4)", source = "Individual")




model_results <- bind_rows(model_results, c3_result, c4_result)

model_results <- model_results %>%
  mutate(significant = !(conf.low <= 0 & conf.high >= 0))


trait_order <- model_results %>%
  filter(source == "Individual") %>%
  arrange(estimate) %>%
  pull(trait) %>%
  unique()

model_results <- model_results %>%
  mutate(trait = factor(trait, levels = trait_order))




ggplot(model_results, aes(x = estimate, y = trait, color = significant)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("gray30", "red3")) +
  labs(
    title = "Effect of NPK Treatment on Plant Traits",
    x = "Estimate (NPK vs Control)",
    y = NULL,
    color = "Significant\n(CI excludes 0)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0, size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )
