###############################################################################
############## AFEC 2022 Dung Beetle Thermal Tolerance  project  ##############
###############################################################################

## Install packages
library(ggplot2)
library(tidyverse) 
library(dbplyr) 
library(readxl) 
library(iNEXT) # to create sampling sufficiency curves
library(vegan)  # for multivariate analysis
library(lme4) # to run the lmer function (linear mixed effects model)
library(lmerTest) # to perform automatic backward model selection of fixed and random parts of the linear mixed effect model

############################################
# SETUP
############################################

## Set working directory
getwd()
setwd("/Volumes/GoogleDrive/My Drive/afec-x/0_miniproject/dung-beetle-thermal-tolerance")

## Import data
beetles_raw <- read_excel("data/thermal_tolerance.xlsx")
head(beetles_raw)
colnames(beetles_raw)

## Clean data
beetles <- rename(beetles_raw,
                  beetleID = 'number',
                  wet_weight = 'wet_body_weight(g)',
                  after_weight = 'after_body_weight(g)',
                  set_temp = 'set(°C)',
                  actual_temp = 'actual(°C)',
                  bodysize_a = 'BSa(mm)',
                  bodysize_b = 'BSb(mm)',
                  wingsize_a = 'WSa(mm)', 
                  bodysize_a = 'BSa(mm)',
                  bodysize_b = 'BSb(mm)',
                  wingsize_a = 'WSa(mm)',
                  wingsize_b = 'WSb(mm)',
                  dry_weight = 'dry_body_weight(g)',
                  subfamily = 'Subfamily'
                  )

beetles$beetleID<- as.factor(beetles$beetleID)
beetles$site <- as.factor(beetles$site)
beetles$species <- as.factor(beetles$species)
beetles$site <- as.factor(beetles$site)
beetles$metric <- as.factor(beetles$metric)
beetles$humidity <- as.factor(beetles$humidity)
beetles$treatment <- as.factor(beetles$treatment)
beetles$subfamily <- as.factor(beetles$subfamily)

## Create CTmax and CTmin tibbles
beetles_max <- beetles %>% filter(str_detect(metric,'CTmax'))
beetles_min <- beetles %>% filter(str_detect(metric,'CTmin'))

## Check histograms by species
beetles %>% filter(str_detect(species,'sp1')) %>%
  ggplot() +
  aes(x=wet_weight)+
  geom_histogram()+ 
  geom_histogram(fill="darkgreen",color="darkgreen")

beetles %>% filter(str_detect(species,'sp2')) %>%
  ggplot() +
  aes(x=wet_weight)+
  geom_histogram()+
  geom_histogram(fill="darkgreen",color="darkgreen")

beetles %>% filter(str_detect(species,'sp3')) %>%
  ggplot() +
  aes(x=wet_weight)+
  geom_histogram()+
  geom_histogram(fill="darkgreen",color="darkgreen")

beetles %>% filter(str_detect(species,'sp3')) %>%
  ggplot() +
  aes(x=actual_temp)+
  geom_histogram()+
  geom_histogram(fill="darkgreen",color="darkgreen")


############################################
# ANALYSIS
############################################

####################
## Hypothesis 1: subtropical vs tropical

####################
## Hypothesis 2: humidity impact on CTmax and CTmin

# Analysis method: multi-comparison t-test and anova

# re-factor and test for significant difference between groups

# Box plot

####################
## Hypothesis 3: rainforest vs rubber plantation

####################
## Hypothesis 4: trait impact on CTmax and CTmin

# Analysis method: linear model with species as a random effect

# CTmax
lm <- lmer(actual_temp~wet_weight*humidity+(1|species),beetles_max)
summary(lm)
anova(lm)  # marginally insignificant

# CTmin
lm2 <- lmer(actual_temp~wet_weight*humidity+(1|species),beetles_min)
summary(lm2)
anova(lm2) # insignificant

# Line plot


