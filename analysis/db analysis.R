###############################################################################
############## AFEC 2022 Dung Beetle Thermal Tolerance  project  ##############
###############################################################################

# Fred Munyao, Xiaoyu Yu, Miranda Li, and Hanchen Shuai
# Advised by Thilina Nimalrathna and Akihiro Nakamura
# China Academy of Sciences - Xishuangbanna Tropical Botanical Garden 
# AFEC-X 2022

# To find outputs for graphs and p-values, search for "!HERE"
# If you have any questions about the code, feel free to contact Miranda at mirandali1995@gmail.com

## Install packages
library(ggplot2)
library(ggpubr)
library(ggbreak)
library(MASS) # for boxcox
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
library(ciTools) # for confidence intervals
library(performance) # for checking diagnostics
library(RColorBrewer) # for beautiful colors
library(grid)

###########################################################################
#### SETUP ####
###########################################################################

#### Set working directory####
getwd()
setwd("/Volumes/GoogleDrive/My Drive/afec-x/0_miniproject/dung-beetle-thermal-tolerance")

## Import data
beetles_raw <- read_excel("data/thermal_tolerance.xlsx")
head(beetles_raw)
colnames(beetles_raw)

## Clean data
beetles <- 
  rename(beetles_raw,
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

# remove species with fewer than 5 samples; too few! 
beetles <- beetles%>% group_by(species) %>% filter(n()>5)
beetles <- droplevels(beetles)

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
# Relevel site so subtropics is first
levels(beetles$site)
levels(beetles$site) <- c("tropics", "subtropics")
beetles$site <- factor(beetles$site, levels = c("subtropics", "tropics"))

# Relevel and rename species
levels(beetles$species)
beetles$species <- factor(beetles$species, levels = c("sp1", "sp2", "sp6", "sp7", "sp8", "sp9", "sp10"))
levels(beetles$species) <- c("s1", "s2", "t1", "t2", "t3", "t4", "t5")

# Create new variables for use in ggplot
beetles$site.humidity <- paste(beetles$site, beetles$humidity)
beetles$site.metric <- paste(beetles$site, beetles$metric)
beetles$site.treatment <- paste(beetles$site, beetles$treatment)
beetles$species.metric <- paste(beetles$species, beetles$metric)
beetles$top_temp <- ifelse(beetles$metric == "CTmax", beetles$actual_temp, -999) # just taking CTmax, for graphing later

beetles$site.humidity <- as.factor(beetles$site.humidity)
beetles$site.metric <- as.factor(beetles$site.metric)
beetles$site.treatment <- as.factor(beetles$site.treatment)
beetles$species.metric <- as.factor(beetles$species.metric)

# shades of blue that we can use for our three different humidity levels
humidity_colors <- c("aliceblue", "cadetblue2", "deepskyblue4")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888") # for Aki-san :)


#####Calculate traits#####

# Body area = use formula for ellipse 
beetles$bodysize <- pi * (beetles$bodysize_a/2) * (beetles$bodysize_b/2)

# Wing size = formula for ellipse
# Relative wing size = longer length of wing / longer length of body
beetles$wingsize <- pi * (beetles$wingsize_a/2) * (beetles$wingsize_b/2)
beetles$relative_wingsize <- beetles$wingsize_a / beetles$bodysize_a

# Water loss over dry weight
# First make water loss 0 or higher...Justification: losing water is how beetles cool down their body. Water gain is probably due to two reasons: 1) condensation for CTmin measurements; 2) maybe human error in measuring CTmax body weight?
beetles$water_loss_pos <- ifelse(beetles$water_loss<=0, 0, beetles$water_loss)

# Water loss proportion = water loss / dry weight
beetles$water_loss_prop <- beetles$water_loss_pos / beetles$dry_weight

#####Check for normality####### 

# First: visually

# Check for normality of response (CTmax and CTmin temperature)
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
  aes(x=relative_wingsize, fill=species, color=species)+
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


# Second: Shapiro test for normality by each species
# This doesn't work if the data is incomplete

num_values <- colnames(select_if(beetles, is.numeric)) # gives all the numeric column names
num_values <- num_values[num_values!="species"] # for some reason species is getting selected as a numerical column...why?
shapiro_result <- data.frame(matrix(ncol = length(num_values), nrow = length(levels(beetles$species)))) # create new results matrix by species and numerical trait
colnames(shapiro_result) <- num_values
rownames(shapiro_result) <- levels(beetles$species)

# define a function to apply shapiro test and extract p-value

shapiro_test_p <- function(values) {
  shapiro.test(values)$p.value
}

for(i in 1:nrow(shapiro_result)) {
  species_i <- filter(beetles,str_detect(species,rownames(shapiro_result)[i])) # create tibble of only this species
  species_i <- species_i[which(sapply(species_i, is.numeric))] # remove all non-numeric values
  shapiro_result[i,] <- apply(species_i,2,shapiro_test_p) # extract pvalue
}

shapiro_result # check results

# Transform water loss proportion
# Note in the final version, we don't use the transformation because when we use dry body weight instead of wet body weight to calculate water loss, the results are more normal

# use boxcox to find best transformation 
par(mfrow=c(1,1))
min(beetles$actual_temp, na.rm = TRUE) # find the minimum temperature
bc <- boxcox(actual_temp+5.5~water_loss_prop, data=beetles) # use the basic model for boxcox

beetles$water_loss_prop_sqrt <- beetles$water_loss_prop^(1/2)

# check the graph. let's not worry about the subtropics, because we'll only use tropics for question 3 (about traits)
beetles %>% 
  ggplot() +
  aes(x=water_loss_prop_sqrt, fill=species, color=species)+
  geom_histogram()+
  facet_wrap(~ species, scale="free")+
  theme_classic()

beetles %>% 
  ggplot() +
  aes(x=water_loss_prop, fill=species, color=species)+
  geom_histogram()+
  facet_wrap(~ species, scale="free")+
  theme_classic()

#####Create new data tables####

# Create CTmax and CTmin tibbles
beetles_max <- beetles %>% filter(str_detect(metric,'CTmax'))
beetles_min <- beetles %>% filter(str_detect(metric,'CTmin'))

# Create tibbles by site
beetles_subtropics <- beetles %>% filter(str_detect(site,'subtropics'))
beetles_subtropics_max <- beetles_subtropics %>% filter(str_detect(metric,'CTmax'))
beetles_subtropics_min <- beetles_subtropics %>% filter(str_detect(metric,'CTmin'))

beetles_tropics <- beetles %>% filter(!str_detect(site,'subtropics')) # for some reason, having "tropics" here doesn't work. So we take what's NOT the subtropics
beetles_tropics_max <- beetles_tropics%>% filter(str_detect(metric,'CTmax'))
beetles_tropics_min <- beetles_tropics%>% filter(str_detect(metric,'CTmin'))

# Create tibbles for by humidity level
beetles_50 <- beetles %>% filter(str_detect(humidity,'50'))
beetles_max_50 <- beetles_max %>% filter(str_detect(humidity,'50'))
beetles_min_50 <- beetles_min %>% filter(str_detect(humidity,'50'))

beetles_30 <- beetles %>% filter(str_detect(humidity,'30'))
beetles_max_30 <- beetles_max %>% filter(str_detect(humidity,'30'))
beetles_min_30 <- beetles_min %>% filter(str_detect(humidity,'30'))

beetles_90 <- beetles %>% filter(str_detect(humidity,'90'))
beetles_max_90 <- beetles_max %>% filter(str_detect(humidity,'90'))
beetles_min_90 <- beetles_min %>% filter(str_detect(humidity,'90'))

beetles_tropics_min_50 <- beetles_tropics_min %>% filter(str_detect(humidity,'50'))
beetles_tropics_min_30 <- beetles_tropics_min %>% filter(str_detect(humidity,'30'))
beetles_tropics_min_90 <- beetles_tropics_min %>% filter(str_detect(humidity,'90'))

beetles_tropics_max_50 <- beetles_tropics_max %>% filter(str_detect(humidity,'50'))
beetles_tropics_max_30 <- beetles_tropics_max %>% filter(str_detect(humidity,'30'))
beetles_tropics_max_90 <- beetles_tropics_max %>% filter(str_detect(humidity,'90'))


#####Create species summary data tables####

# species_temp  table summarizes CTmin and CTmax by species and treatment
# this will be used for some of our linear models later that need averages by species
species_temp <- 
  beetles %>%
  group_by(species, site, metric, humidity) %>%
  summarise(temp = mean(actual_temp,na.rm=TRUE))
species_temp <- pivot_wider(species_temp, names_from=metric, values_from=temp) # make different columns for average ctmin and ctmax
species_temp$range = species_temp$CTmax - species_temp$CTmin # calculate range from CTmax to CTmin
species_temp

# species_traits table summarizes mean traits for each species
species_traits <- 
  beetles %>%
  group_by(species, site) %>%
  summarise(across(everything(), mean, na.rm=TRUE))
species_traits <- species_traits[,colSums(is.na(species_traits))<nrow(species_traits)] # remove columns that are all NA (these are non-summarizable columns)
species_traits <- select(species_traits, -c(set_temp,actual_temp,top_temp)) #remove the temperature columns
species_traits 

#### Summary of data ####

# See what data we have

beetles %>%
  count(species)

beetles %>%
  count(site)

beetles %>%
  group_by(site, treatment) %>%
  count()

beetles_50 %>%
  group_by(site, treatment) %>%
  count()

beetles %>%
  group_by(species, metric, humidity) %>%
  count()

beetles_tropics %>%
  group_by(species, metric, humidity) %>%
  count() %>% print(n=100)

beetles_tropics_max %>%
  group_by(species, metric, humidity) %>%
  count() %>% print(n=100)

beetles_tropics_min %>%
  group_by(species, metric, humidity) %>%
  count() %>% print(n=100)


  

###########################################################################
#### ANALYSIS ####
###########################################################################

#### Question 1 ####

# Question: How do dung beetle thermal tolerances respond in different climate zones?

# Hypothesis 1: Beetles from the subtropics have lower CTmax and CTmin, but larger thermal tolerance range, because of higher temperature range in subtropics.

## ANALYSIS 

## Analysis for CTmin and CTmax
# Linear model: temp ~ site, with species as a random effect, for humidity=50% only
lm_site_max50 <- lmer(actual_temp~site+(1|species),beetles_max_50,na.action=na.omit)
summary(lm_site_max50) # !HERE: p-value for tropics vs subtropics CTmax
r.squaredGLMM(lm_site_max50)
Anova(lm_site_max50) 
check_model(lm_site_max50)

lm_site_min50 <- lmer(actual_temp~site+(1|species),beetles_min_50,na.action=na.omit)
summary(lm_site_min50) # !HERE: p-value for tropics vs subtropics CTmin
r.squaredGLMM(lm_site_min50)
Anova(lm_site_min50) 
check_model(lm_site_min50)

# We used 50% humidity, but let's check if the results are the same for overall humidity
# Overall linear model: temp ~ site, with species and humidity as a random effect
lm_site_max <- lmer(actual_temp~site+(1|species:humidity),beetles_max,na.action=na.omit)
summary(lm_site_max) 
r.squaredGLMM(lm_site_max)
Anova(lm_site_max) 

lm_site_min <- lmer(actual_temp~site+(1|species:humidity),beetles_min,na.action=na.omit)
summary(lm_site_min) 
r.squaredGLMM(lm_site_min)
Anova(lm_site_min) 

# Run analyses for 30 and 90
lm_site_max90 <- lmer(actual_temp~site+(1|species),beetles_max_90,na.action=na.omit)
summary(lm_site_max90) 
r.squaredGLMM(lm_site_max90)
Anova(lm_site_max90) 

lm_site_min90 <- lmer(actual_temp~site+(1|species),beetles_min_90,na.action=na.omit)
summary(lm_site_min90) 
r.squaredGLMM(lm_site_min90)
Anova(lm_site_min90) 

lm_site_max30 <- lmer(actual_temp~site+(1|species),beetles_max_30,na.action=na.omit)
summary(lm_site_max30)  # this is not significantly different! hmm
r.squaredGLMM(lm_site_max30)
Anova(lm_site_max30) 

lm_site_min30 <- lmer(actual_temp~site+(1|species),beetles_min_30,na.action=na.omit)
summary(lm_site_min30) 
r.squaredGLMM(lm_site_min30)
Anova(lm_site_min30) 


## Analysis for range difference

# Method 1: by species averages
# Species average ranges by humidity are in the species_temp tibble
species_temp
species_temp_50 <- species_temp %>% filter(str_detect(humidity,'50')) 

# Now check if they are statistically significantly different
# Linear model: temp_range ~ site

# use 50% humidity only
lm_site_range50 <- lm(range~site,species_temp_50, na.action=na.omit)
summary(lm_site_range50) # !HERE: the range is not significantly different
r.squaredGLMM(lm_site_range50) 
anova(lm_site_range50) 
check_model(lm_site_range50)

# Calculate the average ranges by site
(site_avg_range <- 
    species_temp_50 %>%
    group_by(site) %>%
    summarise(avg_range = mean(range)))

# all humidity?
lm_site_range <- lmer(range~site+(1|species), species_temp, na.action=na.omit)
summary(lm_site_range) # the range is not significantly different here either
r.squaredGLMM(lm_site_range) 
Anova(lm_site_range) 


# Method 2: Random pairs method for range calculation
# Create random pairs for each species, so each species has 5 data points 

range_random_pairs <- 
  beetles %>% filter(str_detect(humidity,'50')) %>%
  group_by(species, site, metric, ID, humidity) %>%
  summarise(temp = mean(actual_temp,na.rm=TRUE))

range_random_pairs_max <- range_random_pairs %>% filter(str_detect(metric,'CTmax'))
range_random_pairs_min <- range_random_pairs %>% filter(str_detect(metric,'CTmin'))

n <- as.integer(length(levels(range_random_pairs$species)))
range_random_pairs_max$ID <- rep(c("1","2","3","4","5"),times=n)
range_random_pairs_max$speciesID <- paste(range_random_pairs_max$species, range_random_pairs_max$ID)
range_random_pairs_max <- rename(range_random_pairs_max, CTmax='temp')
range_random_pairs_min$ID <- rep(c("1","2","3","4","5"),times=n)
range_random_pairs_min$speciesID <- paste(range_random_pairs_min$species, range_random_pairs_min$ID)
range_random_pairs_min <- rename(range_random_pairs_min, CTmin='temp')

range_random_pairs <- left_join(range_random_pairs_max, range_random_pairs_min, by=c("species","site","ID","humidity","speciesID")) # reassign random pairs by merging the CTmax and CTmin
range_random_pairs$range <- range_random_pairs$CTmax - range_random_pairs$CTmin

lmer_range_pairs <- lmer(range~site+(1|species),range_random_pairs)
summary(lmer_range_pairs) # still not significant
Anova(lmer_range_pairs)
r.squaredGLMM(lmer_range_pairs)

## Analysis for difference between species for 50% humidity (so we can get the letters for the plot)
lm_max_species <- lm(actual_temp~species,beetles_max_50,na.action=na.fail)
summary(lm_max_species)
pairs(emmeans(lm_max_species, "species")) # post-hoc test (Tukey) to determine which pairs are significantly different 

lm_min_species <- lm(actual_temp~species,beetles_min_50,na.action=na.fail)
summary(lm_min_species)
pairs(emmeans(lm_min_species, "species")) # post-hoc test (Tukey) to determine which pairs are significantly different 

## PLOTS

# CTmax 
beetles_max_50 %>%
  ggplot() +
  aes(x = site, y = actual_temp, fill=site)+
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  labs(y = "Critical thermal maximum (CTmax) (\u00B0C)", x = element_blank()) +
  ggtitle("CTmax for 50% humidity by site (all species)")+
  theme_classic()

# CTmin 
beetles_min_50 %>%
  ggplot() +
  aes(x = site, y = actual_temp, fill=site)+
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  labs(y = "Critical thermal minimum (CTmin) (\u00B0C)", x = element_blank()) +
  ggtitle("CTmin for 50% humidity by site (all species)")+
  theme_classic()

# !HERE: Output: aggregated 
beetles_50 %>%
  ggplot() +  
  theme_classic()+
  theme(text = element_text(size = 18))+
  aes(x = site, y = actual_temp, group = site.metric, fill=metric)+
  geom_boxplot(position="identity",outlier.shape = NA, alpha = 0.8) + 
  labs(y = "Temperature (\u00B0C)", x = element_blank()) +
  scale_y_break(c(15, 30))

# !HERE: Output: range, using species level 
species_temp_50 %>%
  ggplot() +
  theme_classic()+
  theme(text = element_text(size = 18))+
  aes(x = site, y = range)+
  geom_boxplot(position="identity",outlier.shape = NA, fill="darkseagreen",alpha = 0.8) + 
  labs(y = "CTmax - CTmin (\u00B0C)", x = element_blank())

# Output: plot for species separately, using CTmax to reorder them 

# first, reorder species using CTmax
beetles_50_reordered <- beetles_50 
beetles_50_reordered$reordered_species <- fct_reorder(beetles_50$species, beetles_50$top_temp, mean)
levels(beetles_50_reordered$reordered_species) <- c("s1","s2","t1","t2","t3","t4","t5") # so the x axis can be in order and beautiful
beetles_50_reordered %>% count(species,reordered_species) # note how we reordered them!

# !HERE then plot
beetles_50_reordered %>%
  ggplot() +
  theme_classic()+
  theme(text=element_text(size=18))+
  aes(x = reordered_species, y = actual_temp, group = species.metric, fill = species)+
  geom_boxplot(position="identity", outlier.shape = NA, alpha = 0.8) + 
  labs(y = "Temperature (\u00B0C)", x = "Species") +
  scale_y_break(c(20, 34))+
  scale_fill_manual(values=safe_colorblind_palette)

#### Question 2 ###### 

# Question: How do dung beetle thermal tolerances respond under different humidity gradients?

# Hypothesis: Increasing humidity will decrease CTmax of dung beetles, because it reduces the ability of the beetles to use evapotranspiration to thermoregulate.

## ANALYSIS

# Linear model: CTmax ~ site*humidity+ (1|species) 
# We include site because we want to understand if there is an interaction between humidity and site

# CTmax
lm_humid_max <- lmer(actual_temp~site*humidity+(1|species),beetles_max,  na.action = na.fail)
summary(lm_humid_max)
Anova(lm_humid_max) 
r.squaredGLMM(lm_humid_max)
step(lm_humid_max) # stepwise selection eliminated the interaction effect
dredge(lm_humid_max) # dredge kept everything based on AIC, but the interaction is very insignificant
# post-hoc test (Tukey) to determine which pairs are significantly different 
pairs(emmeans(lm_humid_max, "humidity"))  # 30-50 not significantly different, all others are

# Simplified model
lm_humid_max2 <- lmer(actual_temp~site+humidity+(1|species),beetles_max,  na.action = na.fail)
summary(lm_humid_max2)
Anova(lm_humid_max2) 
r.squaredGLMM(lm_humid_max2)
pairs(emmeans(lm_humid_max2, "humidity"))  # !HERE: pairwise comparison between 50 and 90

# Subtropics only CTmax
lm_humid_max_s <- lmer(actual_temp~humidity+(1|species),beetles_subtropics_max,  na.action = na.fail)
summary(lm_humid_max_s)
Anova(lm_humid_max_s) 
r.squaredGLMM(lm_humid_max_s)
pairs(emmeans(lm_humid_max_s, "humidity")) # all significantly different

# Tropics only CTmax
lm_humid_max_t <- lmer(actual_temp~humidity+(1|species),beetles_tropics_max,  na.action = na.fail)
summary(lm_humid_max_t)
Anova(lm_humid_max_t) 
r.squaredGLMM(lm_humid_max_t)
pairs(emmeans(lm_humid_max_t, "humidity"))

# CTmin
lm_humid_min <- lmer(actual_temp~site*humidity+(1|species),beetles_min, na.action=na.omit)
summary(lm_humid_min)
Anova(lm_humid_min) 
step(lm_humid_min) 
dredge(lm_humid_min) 
pairs(emmeans(lm_humid_min, "humidity"))

#Simplified model
lm_humid_min2 <- lmer(actual_temp~site+humidity+(1|species),beetles_min, na.action=na.omit)
summary(lm_humid_min2)
Anova(lm_humid_min2) 
step(lm_humid_min2) 
dredge(lm_humid_min2) 
pairs(emmeans(lm_humid_min2, "humidity"))

# Subtropics only CTmin
lm_humid_min_s <- lmer(actual_temp~humidity+(1|species),beetles_subtropics_min,  na.action = na.omit)
summary(lm_humid_min_s)
Anova(lm_humid_min_s) 
r.squaredGLMM(lm_humid_min_s)
pairs(emmeans(lm_humid_min_s, "humidity"))

# Tropics only CTmin
lm_humid_min_t <- lmer(actual_temp~humidity+(1|species),beetles_tropics_min,na.action = na.omit)
summary(lm_humid_min_t)
Anova(lm_humid_min_t) 
r.squaredGLMM(lm_humid_min_t)
pairs(emmeans(lm_humid_min_t, "humidity"))


## PLOTS

# Plots with prediction data and confidence intervals (ignoring effects by species)

# Create prediction for CTmax by site
(new_data_max <- 
    beetles_max %>%
    group_by(site, metric, humidity) %>%
    summarise(temp = mean(actual_temp, na.rm=TRUE)))

(pred_max <- add_ci(
  new_data_max,
  lm_humid_max,
  alpha = 0.05,
  includeRanef = FALSE,
  nSims = 500,
))

# Create prediction for CTmin by site
(new_data_min <- 
    beetles_min %>%
    group_by(site, metric, humidity) %>%
    summarise(temp = mean(actual_temp, na.rm=TRUE)))

(pred_min <- add_ci(
  new_data_min,
  lm_humid_min,
  alpha = 0.05,
  includeRanef = FALSE,
  nSims = 500,
))

# join data together and plot
pred_all <- union(pred_max,pred_min)
pred_all

# plotting all on one slide
pred_all %>% 
  ggplot() +
  aes(x = humidity, group = metric,  linetype = site, y = pred, ymin = LCB0.025, ymax = UCB0.975, color=metric)+
  geom_errorbar(width = 0.1) +
  geom_point(size = 1.5)+
  theme_classic() +
  theme(legend.title=element_blank())+
  labs(y = "Temperature (\u00B0C)", x = "Humidity (%RH)")+
  theme(text = element_text(size = 14))
  

# !HERE: plot results
pred_all %>% 
  ggplot() +
  aes(x = humidity, group = metric,  y = pred, ymin = LCB0.025, ymax = UCB0.975, color=metric)+
  geom_errorbar(width = 0.1) +
  geom_point(size = 1.5)+
  theme_classic() +
  theme(legend.title=element_blank())+
  theme(text = element_text(size = 14))+
  labs(y = "Temperature (\u00B0C)", x = "Humidity (%RH)")+
  facet_wrap(metric~site,scale="free")+
  theme(strip.text.x = element_blank(), text=element_text(size=14), panel.spacing = unit(2, "lines"))


# plotting individually
pred_all %>% filter(str_detect(site,'subtropics')) %>%
  ggplot() +
  aes(x = humidity, group = metric,  y = pred, ymin = LCB0.025, ymax = UCB0.975, color=metric)+
  geom_errorbar(width = 0.1) +
  geom_point(size = 1.5)+
  theme_classic() +
  theme(legend.title=element_blank())+
  theme(text = element_text(size = 14))+
  labs(y = "Temperature (\u00B0C)", x = "Humidity (%RH)")+
  scale_y_break(c(3, 33))

pred_all %>% filter(!str_detect(site,'subtropics')) %>% # we want to filter by tropics
  ggplot() +
  aes(x = humidity, group = metric, y = pred, ymin = LCB0.025, ymax = UCB0.975, color=metric)+
  geom_errorbar(width = 0.1) +
  geom_point(size = 1.5)+
  theme_classic() +
  theme(legend.title=element_blank())+
  theme(text = element_text(size = 14))+
  labs(y = "Temperature (\u00B0C)", x = "Humidity (%RH)")+
  scale_y_break(c(12, 40))

# Overall plots with real data
beetles %>%
  ggplot() +
  aes(x = humidity, y = actual_temp, group = treatment, fill = humidity)+
  geom_boxplot(position="identity",outlier.shape = NA, alpha = 0.8) +
  labs(y = "Temperature (\u00B0C)", x = "Humidity (%RH)") +
  ggtitle("CTmax and min by humidity level")+
  theme_classic() +
  scale_y_break(c(20, 30))+
  theme(legend.position = "none") +
  scale_fill_manual(values=humidity_colors)

# Plotting actual values
beetles_max %>%
  ggplot() +
  aes(x = humidity, y = actual_temp, group = treatment, fill = humidity)+
  geom_boxplot(position="identity",outlier.shape = NA, alpha = 0.8) +
  labs(y = "Temperature (\u00B0C)", x = "Humidity (%RH)") +
  ggtitle("CTmax by humidity level")+
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values=humidity_colors)

beetles_min %>%
  ggplot() +
  aes(x = humidity, y = actual_temp, group = treatment, fill = humidity)+
  geom_boxplot(position="identity",outlier.shape = NA, alpha = 0.8) +
  labs(y = "Temperature (\u00B0C)", x = "Humidity (%RH)") +
  ggtitle("CTmin by humidity level")+
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values=humidity_colors)

# Plot just for species 1

# Species 1
# CTmax
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

# CTmin
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

# Both
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

# All species, CTmax
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



#### Question 3 ###### 
# Question: Which morphological traits have significant impact on the thermal tolerance of dung beetles?

# Check for correlation between traits
cor.test(beetles_tropics$wingsize, beetles_tropics$bodysize) # of course this is very high, that's why we use relative wingsize
cor.test(beetles_tropics$relative_wingsize, beetles_tropics$bodysize)
cor.test(beetles_tropics$relative_wingsize, beetles_tropics$water_loss_prop_sqrt)
cor.test(beetles_tropics$bodysize, beetles_tropics$water_loss_prop_sqrt)
# all the other correlations are acceptable

# OPTION 1: analysis at a species level
# The problem with this is that when we condense to a species level, we only get 4 means that we can use, leading to loss of degrees of freedom to estimate parameters in our model

# Calculate CTmax90-CTmax50 for our 4 species
species_delta_t <- 
  beetles %>%
  filter(str_detect(metric,'CTmax')) %>%
  filter(!str_detect(site,'subtropics')) %>%
  group_by(species, site, treatment) %>%
  summarise(temp = mean(actual_temp))

species_delta_t 
species_delta_t <- pivot_wider(species_delta_t, names_from=treatment, values_from=temp) # make different columns for average ctmin and ctmax
species_delta_t$CTmax_90v50 <- species_delta_t$CTmax_50 - species_delta_t$CTmax_90 # we should really only use 90v50 , the other data is not very complete

# Filter out what we don't have the delta for
species_delta_t <- species_delta_t %>% drop_na(CTmax_90v50)

# Combine with traits table
species_delta_t <- left_join(species_delta_t,species_traits)
species_delta_t # check final table

# Run statistical analysis
lm_1 <- lm(CTmax_90v50~water_loss_prop_sqrt,species_delta_t)
summary(lm_1) 
lm_2 <- lm(CTmax_90v50~bodysize,species_delta_t)
summary(lm_2)
lm_3 <- lm(CTmax_90v50~relative_wingsize,species_delta_t)
summary(lm_3)

# not significant, really not enough data points


# OPTION 2: generate random pairs to assess CTmax50 - CTmax90

# we'll use beetles from the tropics, CTmax, humidity 50 and 90 to generate our random pairs
pairs_50 <- beetles_tropics_max_50
pairs_90 <- beetles_tropics_max_90

pairs_50 <- pairs_50 %>% filter(!str_detect(species,"t2")) # this manually removes t2...I'm too lazy to do an automated version where we count n < 10 and remove... :)

# check that we have the same number of observations for each
str(pairs_50)
str(pairs_90)

# let's drop the levels we don't need
pairs_50 <- droplevels(pairs_50)
pairs_90 <- droplevels(pairs_90)

# In order to generate  random pairs, in each of pairs_50 and pairs_90, we reassign the ID column to 1-5 for each species. Then, we'll have two species t1 & ID 1 beetles - one for CTmax_50 and one for CTmax_90 - and so on.
n <- as.integer(length(levels(pairs_50$species)))
pairs_50$ID <- rep(c("1","2","3","4","5"),times=n) 
pairs_50$speciesID <- paste(pairs_50$species, pairs_50$ID)
pairs_90$ID <- rep(c("1","2","3","4","5"),times=n)
pairs_90$speciesID <- paste(pairs_90$species, pairs_90$ID)


# Then we do a left-join of these two tibbles, matching by speciesID. Each row of the final tibble becomes an aggregate of two beetles of the same species: one that underwent CTmax50, one that underwent CTmax90.
pairs <- left_join(pairs_50, pairs_90, by=c("species","site","metric","ID","speciesID")) 

# We calculate the difference between the CTmax at the two humidities
pairs$delta <- pairs$actual_temp.x - pairs$actual_temp.y
hist(pairs$delta) 

# then, we create new traits we care about (bodysize, wingsize, relative wingsize, and water loss prop) by taking the average trait data of the two beetles
pairs$bodysize <- (pairs$bodysize.x+pairs$bodysize.y)/2
pairs$wingsize <- (pairs$wingsize.x+pairs$wingsize.y)/2
pairs$relative_wingsize <- (pairs$relative_wingsize.x+pairs$relative_wingsize.y)/2
pairs$water_loss_prop <- (pairs$water_loss_prop.x+pairs$water_loss_prop.y)/2
pairs$water_loss_prop_sqrt <- (pairs$relative_wingsize.x+pairs$relative_wingsize.y)/2

# check the histograms of these traits
# the previous normality checks were by species, now that we are aggregating everything, we should check it again

pairs_long <- pivot_longer(pairs, c(bodysize:water_loss_prop_sqrt), names_to = "trait")

pairs_long %>%
  ggplot()+
  aes(x = value)+
  geom_histogram(position = "identity") +
  facet_wrap(~ trait, scale = "free")+
  ggtitle("Beetle traits histograms")+
  theme_classic()

# bodysize not that normal, let's transform
hist(pairs$bodysize^(1/3)) # cube root looks ok
pairs$bodysize_cub <- pairs$bodysize^(1/3)

# Now check for significance!
random_water_lm <- lm(delta~water_loss_prop,pairs)
summary(random_water_lm)
anova(random_water_lm)
r.squaredGLMM(random_water_lm)

random_body_lm <- lm(delta~bodysize_cub,pairs)
summary(random_body_lm)
anova(random_body_lm)
r.squaredGLMM(random_body_lm)

random_wing_lm <- lm(delta~relative_wingsize,pairs)
summary(random_wing_lm)
anova(random_wing_lm)
r.squaredGLMM(random_water_lm)

# nothing is significant... let's plot line graphs to understand what's happening
# we see a lot of spread by species in each plot
# the overall trend 

pairs %>%
  arrange(water_loss_prop) %>%
  ggplot() +
  theme_classic()+
  theme(text = element_text(size = 14))+
  aes(x=water_loss_prop,y=delta)+
  geom_point(aes(color=species))+
  scale_color_manual(values=safe_colorblind_palette)+
  labs(y = "CTmax at 50% vs 90% humidity (\u00B0C)", x = "Water loss / Dry body weight")+
  geom_smooth(method="lm")+
  theme(legend.position = "none") +

pairs %>%
  arrange(bodysize) %>%
  ggplot() +
  theme_classic()+
  theme(text = element_text(size = 14))+
  aes(x=bodysize,y=delta)+
  geom_point(aes(color=species))+
  scale_color_manual(values=safe_colorblind_palette)+
  geom_smooth(method="lm")+
  labs(y = element_blank(), x = expression(Body~size~(mm^2)~(cubed)))+
  theme(legend.position = "none") +


pairs %>%
  arrange(relative_wingsize) %>%
  ggplot() +
  theme_classic()+
  theme(text = element_text(size = 14))+
  aes(x=relative_wingsize,y=delta)+
  geom_point(aes(color=species))+
  scale_color_manual(values=safe_colorblind_palette)+
  geom_smooth(method="lm")+
  labs(y = element_blank(), x = "Relative wing size")



# OPTION 3: Analysis at an individual level looking at interaction between traits*humidity
# CTmax ~ traits * humidity + (1|species)

# there is an out outlier with very low CTmax at humidity 30, let's take it out
beetles_tropics_max_clean <- beetles_tropics_max %>% filter(!str_detect(ID,'624'))

# for these, it makes sense to run without species as a random effect because we care about differences BETWEEN species (which have different traits), rather than the differences WITHIN species


##  WATER LOSS

# note, we original used the sqrt transformed water loss metric when we were calculating with wet body weight, but when we use dry body weight, we can just use the un-transformed metric 

# water loss CTmax

lm_water_max <- lm(actual_temp~water_loss_prop*humidity,beetles_tropics_max_clean,na.action=na.omit)
summary(lm_water_max)
Anova(lm_water_max)   #!HERE: humidity is significant, water prop is not
r.squaredGLMM(lm_water_max) 

emtrends(lm_water_max, specs=pairwise~humidity, var="water_loss_prop")
# slopes are not significantly different

# just because we're curious, we saw what happens when we put species as a ranef too :)
lmer_water_max <- lmer(actual_temp~water_loss_prop*humidity+(1|species),beetles_tropics_max_clean,na.action=na.fail)
summary(lmer_water_max)
anova(lmer_water_max) 
r.squaredGLMM(lmer_water_max)
dredge(lmer_water_max) # keeps everything
step(lmer_water_max)
pairs(emmeans(lmer_water_max, "humidity")) 
AIC(lmer_water_max,lm_water_max) 


# BODY SIZE

# Body size CTmax

lm_body_max <- lm(actual_temp~bodysize*humidity,beetles_tropics_max_clean,na.action=na.omit)
summary(lm_body_max) 
anova(lm_body_max) #!HERE
r.squaredGLMM(lm_body_max)
emtrends(lm_body_max, specs=pairwise~humidity, var="bodysize") # slopes are not significantly different

lmer_body_max <- lmer(actual_temp~bodysize*humidity+(1|species),beetles_tropics_max_clean,na.action=na.omit)
summary(lmer_body_max)
Anova(lmer_body_max)  
r.squaredGLMM(lmer_body_max)
dredge(lmer_body_max)
step(lmer_body_max)

# Body size CTmin
lm_body_min <- lm(actual_temp~bodysize*humidity,beetles_tropics_min,na.action=na.omit)
summary(lm_body_min)
Anova(lm_body_min) 
r.squaredGLMM(lm_body_min)

lmer_body_min <- lmer(actual_temp~bodysize*humidity+(1|species),beetles_tropics_min,na.action=na.fail)
summary(lmer_body_min)
Anova(lmer_body_min)  
r.squaredGLMM(lmer_body_min)
dredge(lmer_body_min)
step(lmer_body_min)



# WING SIZE

# Wing size CTmax
lm_wing_max <- lm(actual_temp~relative_wingsize*humidity,data=beetles_tropics_max_clean,na.action=na.omit)
summary(lm_wing_max)
Anova(lm_wing_max)  
r.squaredGLMM(lm_wing_max)

lmer_wing_max <- lmer(actual_temp~relative_wingsize*humidity+(1|species),beetles_tropics_max_clean,na.action=na.fail)
summary(lmer_wing_max)
Anova(lmer_wing_max)  
r.squaredGLMM(lmer_wing_max)


# Wing size CTmin
lm_wing_min <- lm(actual_temp~relative_wingsize*humidity,beetles_tropics_min,na.action=na.omit)
summary(lm_wing_min)
anova(lm_wing_min)   #!HERE
r.squaredGLMM(lm_wing_min)
emtrends(lm_wing_min, specs=pairwise~humidity, var="relative_wingsize") # 50vs90 significant

lmer_wing_min <- lmer(actual_temp~relative_wingsize*humidity+(1|species),beetles_tropics_min,na.action=na.fail)
summary(lmer_wing_min)
Anova(lmer_wing_min)  
r.squaredGLMM(lmer_wing_min)
dredge(lmer_wing_min)
step(lmer_wing_min) 

# is wing size significant for different humidities?
lm_wing_min_30 <- lm(actual_temp~relative_wingsize,beetles_tropics_min_30,na.action=na.omit)
summary(lm_wing_min_30)
Anova(lm_wing_min_30)  
r.squaredGLMM(lm_wing_min_30)

lm_wing_min_50 <- lm(actual_temp~relative_wingsize,beetles_tropics_min_50,na.action=na.omit)
summary(lm_wing_min_50)
Anova(lm_wing_min_50)  
r.squaredGLMM(lm_wing_min_50)

lm_wing_min_90 <- lm(actual_temp~relative_wingsize,beetles_tropics_min_90,na.action=na.omit)
summary(lm_wing_min_90)
Anova(lm_wing_min_90)  
r.squaredGLMM(lm_wing_min_90)




# PLOTS # 

# Check relationships between predictors

beetles_tropics_max %>% 
  arrange(bodysize) %>%
  ggplot() +
  aes(x=bodysize,y=water_loss)+
  geom_point(aes(color=species))+
  geom_smooth(method="lm")+
  scale_x_log10()+
  facet_grid(~humidity, scale = "free")+
  theme_classic()

beetles_tropics_max %>% 
  arrange(bodysize) %>%
  ggplot() +
  aes(x=bodysize,y=water_loss_prop)+
  geom_point(aes(color=species))+
  geom_smooth(method="lm")+
  scale_x_log10()+
  facet_grid(~humidity, scale = "free")+
  theme_classic()

beetles_tropics_max %>% 
  arrange(bodysize) %>%
  ggplot() +
  aes(x=bodysize,y=relative_wingsize)+
  geom_point(aes(color=species))+
  geom_smooth(method="lm")+
  scale_x_log10()+
  theme_classic()


# Water loss

# We do see different slopes for different humidity treatments! At 90, most beetles are unable to lose water at all

#!HERE
beetles_tropics_max %>% 
  filter(!str_detect(ID,'624')) %>% # take out the outlier
  arrange(water_loss_prop) %>%
  ggplot() +
  theme_classic()+
  theme(text = element_text(size = 18))+   
  aes(x=water_loss_prop,y=actual_temp)+
  geom_point(aes(color=species))+
  scale_color_manual(values=safe_colorblind_palette)+
  geom_smooth(method="lm")+
  scale_fill_brewer(palette="Spectral")+
  labs(y = "CTmax (\u00B0C)", x = "Water loss / Dry body weight") +
  facet_grid(~humidity, scale = "free")

# It does not make sense to calculate or plot water loss vs CTmin becuase of so many negative values for water loss for CTmin due to condensation

# Body size
# !HERE CTmax vs body size

beetles_tropics_max %>%
  filter(!str_detect(ID,'624')) %>% # take out the outlier
  arrange(bodysize) %>%
  ggplot() +
  theme_classic()+
  theme(text = element_text(size = 18))+   
  aes(x=bodysize,y=actual_temp)+
  geom_point(aes(color=species))+
  scale_color_manual(values=safe_colorblind_palette)+
  geom_smooth(method="lm")+
  labs(y = "CTmax (\u00B0C)", x = expression(Body~size~(mm^2)))+
  facet_grid(~humidity, scale = "free")

# CTmin vs body size
beetles_tropics_min %>%
  arrange(bodysize) %>%
  ggplot() +
  theme_classic()+
  theme(text = element_text(size = 18))+   
  aes(x=bodysize,y=actual_temp)+
  geom_point(aes(color=species))+
  geom_smooth(method="lm")+
  labs(y = "CTmin (\u00B0C)", x = expression(Body~size~(mm^2)))+
  facet_grid(~humidity, scale = "free")

# Wing size

# CTmax vs wingsize
beetles_tropics_max %>%
  arrange(relative_wingsize) %>%
  ggplot() +
  theme_classic()+
  theme(text = element_text(size = 18))+
  aes(x=wingsize,y=actual_temp)+
  geom_point(aes(color=species))+
  geom_smooth(method="lm")+
  labs(y = "CTmax (\u00B0C)", x = "Wing size / body size")+
  facet_grid(~humidity, scale = "free")

# !HERE CTmin vs wingsize
beetles_tropics_min %>%
  arrange(relative_wingsize) %>%
  ggplot() +
  theme_classic()+
  theme(text = element_text(size = 18))+
  aes(x=relative_wingsize,y=actual_temp)+
  geom_point(aes(color=species))+
  scale_color_manual(values=safe_colorblind_palette)+
  geom_smooth(method="lm")+
  labs(y = "CTmin (\u00B0C)", x = "Relative wing size")+
  facet_grid(~humidity, scale = "free")

#### end ####

