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
library(MuMIn)

###########################################################################
# SETUP
###########################################################################

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

# remove species 3 ; too few samples!
beetles <- beetles[beetles$species != "sp3",]  

# Convert to factors
beetles$site <- as.factor(beetles$site)
beetles$species <- as.factor(beetles$species)
beetles$ID<- as.factor(beetles$ID)
beetles$site <- as.factor(beetles$site)
beetles$metric <- as.factor(beetles$metric)
beetles$humidity <- as.factor(beetles$humidity)
beetles$treatment <- as.factor(beetles$treatment)
beetles$subfamily <- as.factor(beetles$subfamily)

# Rename "rainforest" to "tropics" because no "rubber plantation" for now
# Relevel site so subtropic is first
levels(beetles$site)
levels(beetles$site) <- c("tropics", "subtropics")
beetles$site <- factor(beetles$site, levels = c("subtropics", "tropics"))

# Relevel species to be in order
levels(beetles$species)
beetles$species <- factor(beetles$species, levels = c("sp1", "sp2", "sp6", "sp7", "sp8", "sp9", "sp10"))
levels(beetles$species) <- c("s1", "s2", "t1", "t2", "t3", "t4", "t5")

# Create new factors for use in ggplot
beetles$site.species <- paste(beetles$site, beetles$species)
beetles$site.humidity <- paste(beetles$site, beetles$humidity)
beetles$site.metric <- paste(beetles$site, beetles$metric)
beetles$site.treatment <- paste(beetles$site, beetles$treatment)
beetles$species.metric <- paste(beetles$species, beetles$metric)
beetles$site.humidity <- as.factor(beetles$site.humidity)
beetles$site.metric <- as.factor(beetles$site.metric)
beetles$site.treatment <- as.factor(beetles$site.treatment)
beetles$species.metric <- as.factor(beetles$species.metric)


## See what data we have
beetles %>%
  count(species)

beetles %>%
  count(species, humidity)

beetles %>%
  count(species, metric, humidity)

# shades of blue that we can use for our three different humidity levels
humidity_colors <- c("aliceblue", "cadetblue2", "deepskyblue4")

############################################
## Calculate traits

# Body area = use formula for ellipse 
beetles$bodysize <- pi * (beetles$bodysize_a/2) * (beetles$bodysize_b/2)

# TODO: discuss which wingsize trait to use
beetles$wingsize <- pi * (beetles$wingsize_a/2) * (beetles$wingsize_b/2)

# Water loss over dry weight
# Make water loss 0 or higher
# Justification: we assume losing water is the method to cool down their body. Water gain is probably due to two reasons: 1) condensation for CTmin measurements; 2) maybe human error in measuring CTmax body weight?
beetles$water_loss_pos <- ifelse(beetles$water_loss<=0, 0, beetles$water_loss)

# Water loss proportion = water loss / dry weight
# TODO: subsitute wet weight for dry weight once we've measured it
beetles$water_loss_prop <- beetles$water_loss_pos / beetles$wet_weight

############################################
## Check for normality

# First: visually

# check for normality of response (CTmax and CTmin temperature)
# it looks ok, but  will look better once we have more data points
beetles %>% filter(str_detect(metric,'CTmax')) %>%
  ggplot() +
  aes(x=actual_temp, fill=species, color=species)+
  geom_histogram()+
  ggtitle("CTmax")+
  theme_classic()

beetles %>% filter(str_detect(metric,'CTmin')) %>%
  ggplot() +
  aes(x=actual_temp, fill=species, color=species)+
  geom_histogram()+
  ggtitle("CTmin")+
  theme_classic()

# Check for normality of predictors

# Create a histogram
beetles_long <- pivot_longer(beetles, c(wet_weight:water_loss,bodysize_a:dry_weight,bodysize:water_loss_prop), names_to = "trait")

ggplot(beetles_long, 
  aes(x = value, fill=species, color=species))+
  geom_histogram(position = "identity") +
  facet_wrap(~ trait, scale = "free")+
  ggtitle("Beetle traits histograms")+
  theme_classic()

# Check individual histograms

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
# maybe cube root transformaiton
beetles %>% 
  ggplot() +
  aes(x=water_loss_pos, fill=species, color=species)+
  geom_histogram()+
  facet_wrap(~ species, scale="free")+
  theme_classic()
# TODO: maybe we should  transform water loss proportion


# Second: Shapiro test for normality by each species
# This doesn't work if the data is incomplete... need to fix this

num_values <- colnames(select_if(beetles, is.numeric)) # gives all the numeric column names
shapiro_result <- data.frame(matrix(ncol = length(num_values), nrow = length(levels(beetles$species)))) # create new results matrix by species and numerical trait
colnames(shapiro_result) <- num_values
rownames(shapiro_result) <- levels(beetles$species)

# define a function to apply shapiro test and extract p-value
shapiro_test_p <- function(values) {
  return(tryCatch(shapiro.test(values)$p.value, error=function(e) NULL))
}

for(i in 1:nrow(shapiro_result)) {
  species_i <- filter(beetles,str_detect(species,rownames(shapiro_result)[i])) # create tibble of only this species
  species_i <- species_i[which(sapply(species_i, is.numeric))] # remove all non-numeric values
  shapiro_result[i,] <- apply(species_i,2,shapiro_test_p) # extract pvalue
}

shapiro_result # check results
# water prop is the most not normal

############################################
## Create new data tables

# Create CTmax and CTmin tibbles
beetles_max <- beetles %>% filter(str_detect(metric,'CTmax'))
beetles_min <- beetles %>% filter(str_detect(metric,'CTmin'))

# Create tibbles for humidity 50
beetles_max_50 <- beetles_max %>% filter(str_detect(humidity,'50'))
beetles_min_50 <- beetles_min %>% filter(str_detect(humidity,'50'))
beetles_50 <- beetles %>% filter(str_detect(humidity,'50'))

############################################
## Create species summary data tables

species_temp <- 
  beetles %>%
  group_by(species, site, metric, humidity) %>%
  summarise(temp = mean(actual_temp))
species_temp
    
species_temp <- pivot_wider(species_temp, names_from=metric, values_from=temp) # make different columns for average ctmin and ctmax

species_temp$range = species_temp$CTmax - species_temp$CTmin # calculate range from CTmax to CTmin
species_temp

species_traits <- 
  beetles %>%
  group_by(species, site) %>%
  summarise(across(everything(), mean, na.rm=TRUE))

species_traits <- species_traits[,colSums(is.na(species_traits))<nrow(species_traits)] # remove columns that are all NA (these are non-summarizable columns)
species_traits <- select(species_traits, -c(set_temp,actual_temp)) #remove the temperature columns
species_traits # this table summarizes the average traits of each species

###########################################################################
# ANALYSIS
###########################################################################

############################################
## Question 1: 

# Hypothesis 1: Beetles from the subtropics have lower CTmax and CTmin, but larger thermal tolerance range, because of higher temperature range in subtropics.

# Analysis method: linear model: temp ~ site, with species as a random effect

lm_site_max <- lmer(actual_temp~site+(1|species),beetles_max_50,na.action=na.omit)
summary(lm_site_max) # HERE: p-value for tropics vs subtropics CTmax
r.squaredGLMM(lm_site_max)
Anova(lm_site_max) 

lm_site_min <- lmer(actual_temp~site+(1|species),beetles_min_50,na.action=na.omit)
summary(lm_site_min) # HERE: p-value for tropics vs subtropics CTmin
r.squaredGLMM(lm_site_min)
Anova(lm_site_min) 

species_temp_50 <- species_temp %>% filter(str_detect(humidity,'50')) # should we use 50 only?
lm_site_range <- lm(range~site,species_temp_50, na.action=na.omit)
summary(lm_site_range) # the range is not significantly different
r.squaredGLMM(lm_site_range) 
anova(lm_site_range) 

lm_site_range2 <- lmer(range~site+(1|species), species_temp, na.action=na.omit)
summary(lm_site_range2)
r.squaredGLMM(lm_site_range2) 
Anova(lm_site_range2) 

#take the range and present it? 
# 


# PLOTS

# Plot all CTmax and CTmin for each species separately
beetles_50 %>%
  ggplot() +
  aes(x = site.species, y = actual_temp, group = species.metric, fill = species)+
  geom_boxplot(position="identity", outlier.shape = NA, alpha = 0.8) + 
  labs(y = "Temperature (\u00B0C)", x = element_blank()) +
  ggtitle("CTmax and CTmin by site at 50% humidity")+
  scale_y_break(c(20, 34))+
  theme_classic()

# CTmax 
beetles_max_50 %>%
  ggplot() +
  aes(x = site, y = actual_temp)+
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  labs(y = "Critical thermal maximum (CTmax) (\u00B0C)", x = element_blank()) +
  ggtitle("CTmax for 50% humidity by site (all species)")+
  theme_classic()

# CTmin 
beetles_min_50 %>%
  ggplot() +
  aes(x = site, y = actual_temp)+
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  labs(y = "Critical thermal minimum (CTmin) (\u00B0C)", x = element_blank()) +
  ggtitle("CTmin for 50% humidity by site (all species)")+
  theme_classic()

# Everything
beetles_50 %>%
  ggplot() +
  aes(x = site, y = actual_temp, group = site.metric, fill=site)+
  geom_boxplot(position="identity",outlier.shape = NA, alpha = 0.8) + 
  labs(y = "Temperature (\u00B0C)", x = element_blank()) +
  ggtitle("CTmax and CTmin for 50% humidity by site (all species)")+
  theme(legend.position = "none") +
  scale_y_break(c(15, 30))+
  theme_classic()

############################################
## Question 2: How do dung beetle thermal tolerances respond under different humidity gradients?

# Hypothesis: Increasing humidity will decrease CTmax of dung beetles, because it reduces the ability of the beetles to use evapotranspiration to thermoregulate.

# Analysis method: linear model
# CTmax ~ water loss proportion * humidity  + (1|species)
# where water loss proportion = water loss / dry body mass

# third hypothesis: size as an additional predictor in the model? but there's a correlation between size and site
# for size analysis, just exclude the subtropic beetles

# 

# CTmax
lm_humid_max <- lmer(actual_temp~site*humidity+(1|species),beetles_max,  na.action = na.fail)
summary(lm_humid_max)
Anova(lm_humid_max) 
r.squaredGLMM(lm_humid_max)
step(lm_humid_max) # stepwise selection eliminated the interaction effect and water loss
dredge(lm_humid_max) # dredge kept everything based on AIC

lm_humid_max2 <- lmer(actual_temp~humidity+(1|species),beetles_max,na.action=na.fail)
summary(lm_humid_max2)
anova(lm_humid_max2) 
r.squaredGLMM(lm_humid_max2)

AIC(lm_humid_max, lm_humid_max2)

# post-hoc test (Tukey) to determine which pairs are significantly different 
pairs(emmeans(lm_humid_max, "humidity")) 
pairs(emmeans(lm_humid_max2, "humidity")) 

# CTmin
lm_humid_min <- lmer(actual_temp~humidity+(1|species),beetles_min, na.action=na.omit)
summary(lm_humid_min)
Anova(lm_humid_min) 
step(lm_humid_min) 
pairs(emmeans(lm_humid_min, "humidity"))


# QUESTIONS FOR KYLE:
# how to think about species as a random effect?
# how do we graph the results when combining between different species makes it look not significant when the result actually is? 
# should we just plot species by species?


# PLOTS

# make a prediction plot with confidence intervals instead
# that will show overall effect corrected for species

# Overall
beetles_max %>%
  ggplot() +
  aes(x = humidity, y = actual_temp, group = treatment, fill = humidity)+
  geom_boxplot(position="identity",outlier.shape = NA, alpha = 0.8) +
  labs(y = "Temperature (\u00B0C)", x = "Humidity (%RH)") +
  ggtitle("CTmax and min by humidity level")+
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values=humidity_colors)

# box plot for species 1 only
beetles_max %>% filter(str_detect(species,'s1')) %>%
  ggplot() +
  aes(x = humidity, y = actual_temp, fill = humidity)+
  #geom_jitter(aes(color = humidity), size = 3, alpha = 0.8, width = 0.1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  labs(y = "Critical thermal maximum (CTmax) (\u00B0C)", x = "Humidity (%RH)") +
  ggtitle("Species 1 CTmax")+
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values=humidity_colors)+
  annotate("text", x = 1, y = 38, size=6, label = "a")+
  annotate("text", x = 2, y = 38, size=6, label = "b")+ 
  annotate("text", x = 3, y = 38, size=6, label = "b")

# box plot for species 1 only
beetles_min %>% filter(str_detect(species,'s1')) %>%
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
beetles %>% filter(str_detect(species,'s1')) %>%
  ggplot() +
  aes(x = humidity, y = actual_temp, group = treatment, fill=humidity)+
  geom_boxplot(position = "identity", outlier.shape = NA, alpha = 0.8) + 
  labs(y = "Temperature (\u00B0C)", x = "Humidity (%RH)") +
  ggtitle("Species 1 CTmax and CTmin")+
  theme_classic() +
  scale_fill_manual(values=humidity_colors)+
  scale_y_break(c(3, 32))+
  annotate("text", x = 1, y = 38, size=6, label = "a")+
  annotate("text", x = 2, y = 38, size=6, label = "b")+ 
  annotate("text", x = 3, y = 38, size=6, label = "b")+
  annotate("text", x = 1, y = 2, size=6, label = "a")+
  annotate("text", x = 2, y = 2, size=6, label = "a")+ 
  annotate("text", x = 3, y = 2, size=6, label = "a")


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
  facet_wrap(~ species, scale = "free")



# Water loss line plot
# If the interaction is significant, we expect to see significantly different slopes for different humidity treatments?

beetles_max %>% filter(str_detect(species,'s1')) %>%
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


############################################
## Question 3: Which morphological traits have significant impact on the thermal tolerance of dung beetles?

# species as a ranef, delta ~ traits ?  traits*humidity
# when we go for individual, species as ranef may take all expl power


species_t <- 
  beetles %>%
  filter(str_detect(metric,'CTmax')) %>%
  group_by(species, site, treatment) %>%
  summarise(temp = mean(actual_temp))

species_t
species_t <- pivot_wider(species_t, names_from=treatment, values_from=temp) # make different columns for average ctmin and ctmax
species_t$CTmax_90v50 <- species_t$CTmax_90 - species_t$CTmax_50
species_t$CTmax_50v30 <- species_t$CTmax_50 - species_t$CTmax_30
species_t$CTmax_90v30 <- species_t$CTmax_90 - species_t$CTmax_30
# we should really only use 90v50 , the other data is not very complete

species_traits

# need to get rid of the NaN values 
rda.beetles <-dbrda(species_t~., species_traits, distance="bray", scale=T, na.action = na.fail)

# Analysis method: linear model:
# CTmax ~ body size + (1 | species:humidity)
# CTmin ~ wing size + (1 | species:humidity)

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
beetles %>%
  arrange(bodysize) %>%
  ggplot() +
  aes(x=bodysize,y=actual_temp)+
  geom_point()+
  geom_smooth(method="lm")+
  labs(y = "Critical thermal maximum (CTmax) (\u00B0C)", x = expression(Body~size~(mm^2)))+
  ggtitle("CTmax vs Body size by species and humidity treatment ")+
  facet_grid(~humidity, scale = "free")+
  theme_classic()

# line plot for CTmin vs wingsize
beetles %>%
  arrange(wingsize) %>%
  ggplot() +
  aes(x=wingsize,y=actual_temp)+
  geom_point()+
  geom_smooth(method="lm")+
  labs(y = "Critical thermal minimum (CTmin (\u00B0C)", x = expression(Wing~size~(mm^2)))+
  ggtitle("CTmin vs Wing size by species and humidity treatment ")+
  scale_x_log10()+
  facet_grid(~humidity, scale = "free")+
  theme_classic()


# check correlations between traits
