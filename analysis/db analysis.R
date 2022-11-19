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
library(emmeans) # for post-hoc pairwise tests
library(car) # to run the Anova function that calculates p values for each factor (ANOVA table)


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

## Check histograms

# check for normality of response (CTmax and CTmin temperature)
# it looks ok, but  will look better once we have more data points
beetles_max %>%
  ggplot() +
  aes(x=actual_temp, fill=species, color=species)+
  geom_histogram()+
  theme_classic()

beetles_min %>%
  ggplot() +
  aes(x=actual_temp, fill=species, color=species)+
  geom_histogram()+
  theme_classic()

# check for normality of predictors
beetles_long <- pivot_longer(beetles, c(wet_weight:after_weight,bodysize_a:dry_weight), names_to = "trait")

ggplot(beetles_long, 
  aes(x = value, fill=species, color=species))+
  geom_histogram(position = "identity") +
  facet_wrap(~ trait, scale = "free")+
  theme_classic()

# again, this should look better once we have more data

# check wet weight histogram
beetles %>% 
  ggplot() +
  aes(x=wet_weight, fill=species, color=species)+
  geom_histogram()+
  facet_wrap(~ species, scale="free")+
  theme_classic()


############################################
# ANALYSIS
############################################

####################
## Hypothesis 2 and 4: humidity and trait impact on CTmax and CTmin

# Analysis method: linear model with species as a random effect

# TODO / FOR DISCUSSION WITH GROUP: how shall we calculate wing size and body size? since we have a and b measurements, should we calculate the area of an ellipse?

# CTmax
lm_max <- lmer(actual_temp~wet_weight*humidity+(1|species),beetles_max)
summary(lm_max)
Anova(lm_max)   # TODO: I am still not sure what the difference is between Anova (from car package) and regular anova... need to look this up
step(lm_max) # eliminated the interaction effect and wet_weight

lm_max2 <- lmer(actual_temp~humidity+(1|species),beetles_max)
summary(lm_max2)
anova(lm_max2) 

# post-hoc test (Tukey) to determine which pairs are significantly different 
pairs(emmeans(lm_max2, "humidity")) # 30 vs 50 and 30 vs 90significant, 50 vs 90 not significant

# box plot for species 1 only
beetles_max %>% filter(str_detect(species,'sp1')) %>%
  ggplot() +
  aes(x = humidity, y = actual_temp, fill = humidity)+
  geom_jitter(aes(color = humidity), size = 3, alpha = 0.8, width = 0.1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) + 
  labs(y = "Critical thermal maximum (CTmax)", x = "Humidity (%RH)") +
  ggtitle("Species 1 CTmax")+
  theme_classic() +
  theme(legend.position = "none") +
  annotate("text", x = 1, y = 40, size=6, label = "a")+
  annotate("text", x = 2, y = 40, size=6, label = "b")+ 
  annotate("text", x = 3, y = 40, size=6, label = "b")

# line plot for CTmax vs wet body weight
beetles_max %>%
  arrange(wet_weight) %>%
  ggplot() +
  aes(x=wet_weight,y=actual_temp, color=species)+
  geom_point()+
  geom_smooth(method="lm")+
  labs(y = "Critical thermal maximum (CTmax)", x = "Wet body weight (g)") +  
  ggtitle("CTmax vs Body weight by species and humidity treatment ")+
  scale_x_log10()+
  facet_grid(humidity~species, scale = "free")

# CTmin
lm_min <- lmer(actual_temp~wet_weight*humidity+(1|species),beetles_min)
summary(lm_min)
Anova(lm_min) # all insignificant
step(lm_min) # gets rid of everything
pairs(emmeans(lm_min, "humidity")) # indeed, all insignificant

# box plot for species 1 only
beetles_min %>% filter(str_detect(species,'sp1')) %>%
  ggplot() +
  aes(x = humidity, y = actual_temp, fill = humidity)+
  geom_jitter(aes(color = humidity), size = 3, alpha = 0.8, width = 0.1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) + 
  labs(y = "Critical thermal minimum (CTmin)", x = "Humidity (%RH)") +
  ggtitle("Species 1 CTmin")+
  theme_classic() +
  theme(legend.position = "none") +
  annotate("text", x = 1, y = 2, size=6, label = "a")+
  annotate("text", x = 2, y = 2, size=6, label = "a")+ 
  annotate("text", x = 3, y = 2, size=6, label = "a")

# line plot for CTmin vs wet body weight
beetles_min %>%
  arrange(wet_weight) %>%
  ggplot() +
  aes(x=wet_weight,y=actual_temp, color=species)+
  geom_point()+
  geom_smooth(method="lm")+
  labs(y = "Critical thermal minimum (CTmin)", x = "Wet body weight (g)") +  
  ggtitle("CTmin vs Body weight by species and humidity treatment ")+
  scale_x_log10()+
  facet_grid(humidity~species, scale = "free")

####################
## Hypothesis 1: subtropical vs tropical

# Awaiting data

####################
## Hypothesis 3: rainforest vs rubber plantation

# Awaiting data
