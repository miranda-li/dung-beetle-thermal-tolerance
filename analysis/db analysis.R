###############################################################################
############## AFEC 2022 Dung Beetle Thermal Tolerance  project  ##############
###############################################################################

## Install packages
library(ggplot2)
library(ggpubr)
library(ggbreak)
library(tidyverse) 
library(dbplyr) 
library(readxl) 
library(iNEXT) # to create sampling sufficiency curves
library(vegan)  # for multivariate analysis
library(lme4) # to run the lmer function (linear mixed effects model)
library(lmerTest) # to perform automatic backward model selection of fixed and random parts of the linear mixed effect model
library(emmeans) # for post-hoc pairwise tests
library(car) # to run the Anova function that calculates p values for each factor

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
                  wet_weight = 'wet_body_weight(g)',
                  after_weight = 'after_body_weight(g)',
                  water_loss = 'water_loss(g)',
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

beetles$ID<- as.factor(beetles$ID)
beetles$site <- as.factor(beetles$site)
beetles$species <- as.factor(beetles$species)
beetles$site <- as.factor(beetles$site)
beetles$metric <- as.factor(beetles$metric)
beetles$humidity <- as.factor(beetles$humidity)
beetles$treatment <- as.factor(beetles$treatment)
beetles$subfamily <- as.factor(beetles$subfamily)

# Calculate traits
# TODO: discuss with team how to calculate these traits! For now I'm just using the equation for an ellipse

# Body area = use formula for ellipse ?
beetles$bodysize <- pi * (beetles$bodysize_a/2) * (beetles$bodysize_b/2)

# Wing area = use formula for ellipse ?
beetles$wingsize <- pi * (beetles$wingsize_a/2) * (beetles$wingsize_b/2)

# Water loss over dry weight
# TODO: discuss: should we make water loss 0 or higher?
beetles$water_loss_pos <- ifelse(beetles$water_loss<=0, 0, beetles$water_loss)

# TODO: subsitute wet weight for dry weight once we've measured it
beetles$water_loss_prop <- beetles$water_loss_pos / beetles$wet_weight

# check for normality of predictors
beetles_long <- pivot_longer(beetles, c(wet_weight:water_loss,bodysize_a:dry_weight,bodysize:water_loss_prop), names_to = "trait")

ggplot(beetles_long, 
  aes(x = value, fill=species, color=species))+
  geom_histogram(position = "identity") +
  facet_wrap(~ trait, scale = "free")+
  ggtitle("Beetle traits histograms")+
  theme_classic()


# Check histograms of individual predictors that we're interested in

# body size 
beetles %>% 
  ggplot() +
  aes(x=bodysize, fill=species, color=species)+
  geom_histogram()+
  facet_wrap(~ species, scale="free")+
  theme_classic()

# wing size 
beetles %>% 
  ggplot() +
  aes(x=wingsize, fill=species, color=species)+
  geom_histogram()+
  facet_wrap(~ species, scale="free")+
  theme_classic()

# water loss proportion
beetles %>% 
  ggplot() +
  aes(x=water_loss_prop, fill=species, color=species)+
  geom_histogram()+
  facet_wrap(~ species, scale="free")+
  theme_classic()
# TODO: maybe we should  transform water loss proportion

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


# shades of blue that we can use for our three different humidity levels
humidity_colors <- c("aliceblue", "cadetblue2", "deepskyblue4")

# p-value comparisons for ggplot
humidity_comparisons <- list( c("30", "90"), c("30", "50"), c("50", "90") )


############################################
# ANALYSIS
############################################

####################
## Question 1: 

# Hypothesis 1: Beetles from the subtropics have lower CTmax and CTmin, but larger thermal tolerance range, because of higher temperature range in subtropics.

# Analysis method: linear model: temp ~ site, with species and humidity as random effects

lm_site <- lmer(actual_temp~site+(1|species:humidity),beetles_max,na.action=na.omit)
summary(lm_site)
Anova(lm_site) 
step(lm_site) 

beetles_max %>%
  ggplot() +
  aes(x = site, y = actual_temp, fill = humidity)+
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  labs(y = "Critical thermal maximum (CTmax) (\u00B0C)", x = element_blank()) +
  ggtitle("CTmax by site (all species)")+
  theme_classic() +
  scale_fill_manual(values=humidity_colors)
#  stat_compare_means(method = "anova", comparisons = humidity_comparisons)

####################
## Question 2: How do dung beetle thermal tolerances respond under different humidity gradients?

# Hypothesis: Increasing humidity will decrease CTmax of dung beetles, because it reduces the ability of the beetles to use evapotranspiration to thermoregulate.

# Analysis method: linear model
# CTmax ~ water loss proportion * humidity  + (1|species)
# where water loss proportion = water loss / dry body mass

# CTmin ~ humidity + (1|species) ; it does not makes sense to test water loss proportion for CTmin (most of the water loss values are negative)

# CTmax

lm_max <- lmer(actual_temp~water_loss_prop*humidity+(1|species),beetles_max,na.action=na.omit)
summary(lm_max)
Anova(lm_max)   # TODO: I am still not sure what the difference is between Anova (from car package) and regular anova... need to look this up
step(lm_max) # eliminated the interaction effect and water loss

lm_max2 <- lmer(actual_temp~humidity+(1|species),beetles_max) 
summary(lm_max2)
anova(lm_max2) 

# post-hoc test (Tukey) to determine which pairs are significantly different 
pairs(emmeans(lm_max2, "humidity")) # all significant when we consider all species

# box plot for species 1 only
beetles_max %>% filter(str_detect(species,'sp1')) %>%
  ggplot() +
  aes(x = humidity, y = actual_temp, fill = humidity)+
  #geom_jitter(aes(color = humidity), size = 3, alpha = 0.8, width = 0.1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  labs(y = "Critical thermal maximum (CTmax) (\u00B0C)", x = "Humidity (%RH)") +
  ggtitle("Species 1 CTmax")+
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values=humidity_colors)+
  annotate("text", x = 1, y = 40, size=6, label = "a")+
  annotate("text", x = 2, y = 40, size=6, label = "b")+ 
  annotate("text", x = 3, y = 40, size=6, label = "b")
#  stat_compare_means(method = "anova", comparisons = humidity_comparisons)

# box plot for all species
beetles_max %>%
  ggplot() +
  aes(x = humidity, y = actual_temp, fill = humidity)+
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  labs(y = "Critical thermal maximum (CTmax) (\u00B0C)", x = "Humidity (%RH)") +
  ggtitle("CTmax by humidity")+
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values=humidity_colors)+
  #stat_compare_means(comparisons = humidity_comparisons)+ # Add pairwise comparisons p-value
  facet_wrap(~ species, scale = "free")

# need to do some kind of box plot for many different species

# CTmin
lm_min <- lmer(actual_temp~water_loss_prop*humidity+(1|species),beetles_min, na.action=na.omit) # water loss for CTmin does not make sense...most of them are negative anyways. Ignore
summary(lm_min)
Anova(lm_min) # all insignificant
step(lm_min) # gets rid of everything
pairs(emmeans(lm_min, "humidity")) # indeed, all insignificant

# box plot for species 1 only
beetles_min %>% filter(str_detect(species,'sp1')) %>%
  ggplot() +
  aes(x = humidity, y = actual_temp, fill = humidity)+
  #geom_jitter(aes(color = humidity), size = 3, alpha = 0.8, width = 0.1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  labs(y = "Critical thermal minimum (CTmin) (\u00B0C)", x = "Humidity (%RH)") +
  ggtitle("Species 1 CTmin")+
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values=humidity_colors)+
  annotate("text", x = 1, y = 2, size=6, label = "a")+
  annotate("text", x = 2, y = 2, size=6, label = "a")+ 
  annotate("text", x = 3, y = 2, size=6, label = "a")

# box plot combining CTmax and CTmin for species 1
beetles %>% filter(str_detect(species,'sp1')) %>%
  ggplot() +
  aes(x = humidity, y = actual_temp, group = treatment, fill=humidity)+
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  labs(y = "Temperature (\u00B0C)", x = "Humidity (%RH)") +
  ggtitle("Species 1 CTmax and CTmin")+
  theme_classic() +
  scale_fill_manual(values=humidity_colors)+
  scale_y_break(c(3, 32))+
  annotate("text", x = 1, y = 40, size=6, label = "a")+
  annotate("text", x = 2, y = 40, size=6, label = "b")+ 
  annotate("text", x = 3, y = 40, size=6, label = "b")+
  annotate("text", x = 1, y = 2, size=6, label = "a")+
  annotate("text", x = 2, y = 2, size=6, label = "a")+ 
  annotate("text", x = 3, y = 2, size=6, label = "a")


# Water loss line plot
# If the interaction is significant, we expect to see significantly different slopes for different humidity treatments

beetles_max %>% filter(str_detect(species,'sp1')) %>%
  arrange(water_loss_prop) %>%
  ggplot() +
  aes(x=water_loss_prop,y=actual_temp, color=species)+
  geom_point()+
  geom_smooth(method="lm")+
  labs(y = "Critical thermal maximum (CTmax) (\u00B0C)", x = "Water loss / wet body weight") +  
  ggtitle("CTmax vs water loss proportion by humidity")+
  scale_x_log10()+
  facet_grid(humidity~species, scale = "free")+
  theme_classic()


####################
## Question 3: Which morphological traits have significant impact on the thermal tolerance of dung beetles?

# Analysis method: linear model:
# CTmax ~ body size + (1 | species:humidity)
# CTmin ~ wing size + (1 | species:humidity)

# FIRST FILTER OUT THE NA (temporary as we haven't finished measuring the traits)
# TODO: remove this line of code for the final version
beetles_max_na<-subset(beetles_max, select = -c(dry_weight)) # first remove dry weight because they are all NA
beetles_max_na<-na.omit(beetles_max_na) # then remove the rows with NA (haven't measured bodysize and wingsize yet)
beetles_min_na<-subset(beetles_min, select = -c(dry_weight)) # first remove dry weight because they are all NA
beetles_min_na<-na.omit(beetles_min_na) # then remove the rows with NA  (haven't measured bodysize and wingsize yet)

# CTmax and body size
lm_body <- lmer(actual_temp~bodysize+(1|species:humidity),beetles_max_na,na.action=na.omit)
summary(lm_body)
Anova(lm_body) 
step(lm) # nothing is significant

# CTmin and wingsize
lm_wing <- lmer(actual_temp~wingsize+(1|species:humidity),beetles_min_na,na.action=na.omit)
summary(lm_wing)
Anova(lm_wing) 
step(lm_wing) # nothing is significant

# line plot for CTmax vs bodysize
beetles_max_na %>%
  arrange(bodysize) %>%
  ggplot() +
  aes(x=bodysize,y=actual_temp, color=species)+
  geom_point()+
  geom_smooth(method="lm")+
  labs(y = "Critical thermal maximum (CTmax) (\u00B0C)", x = expression(Body~size~(mm^2)))+
  ggtitle("CTmax vs Body size by species and humidity treatment ")+
  scale_x_log10()+
  facet_grid(humidity~species, scale = "free")+
  theme_classic()

# line plot for CTmin vs wingsize
beetles_min_na %>%
  arrange(wingsize_a) %>%
  ggplot() +
  aes(x=wingsize,y=actual_temp, color=species)+
  geom_point()+
  geom_smooth(method="lm")+
  labs(y = "Critical thermal minimum (CTmin (\u00B0C)", x = expression(Wing~size~(mm^2)))+
  ggtitle("CTmin vs Wing size by species and humidity treatment ")+
  scale_x_log10()+
  facet_grid(humidity~species, scale = "free")+
  theme_classic()
