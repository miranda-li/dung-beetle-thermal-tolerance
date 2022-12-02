

library(ggplot2)
library(tidyverse) 
library(dbplyr) 

getwd()
setwd("/Volumes/GoogleDrive/My Drive/afec-x/0_miniproject/dung-beetle-thermal-tolerance")

#### Plot temperature in machine vs box ####

temp_raw<-read_csv(file="data/dungbeetle.csv",skip=1)

colnames(temp_raw)
temp <- 
  rename(temp_raw,
         num = "#",
         time = "Date Time, GMT+08:00",
         box_temperature = "air temperature in the small box",
         air_temperature = "air temperature",
         body_temperature = "body temperature"
  )

colnames(temp)
str(temp)
temp


temp[which.max(temp$air_temperature),] # find out which one the max is


temp %>%
  ggplot()+
  theme_bw()+
  theme(text= element_text(size=14),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_point(aes(x=num,y=air_temperature,color="#1e7099"),size=1)+
  geom_point(aes(x=num,y=box_temperature,color="#afcfd4"),size=1)+
  geom_point(aes(x=num,y=body_temperature,color="#b45560"),size=1)+
  labs(y='temperature (\u00B0C)',x='time',title='Comparison between different temperatures')+
  xlim(934,3687)+
  scale_color_manual(values=c("#1e7099", "#afcfd4", "#b45560"),
                     name=" ",labels=c('air temperature','box temperature','beetle body temperature'))
