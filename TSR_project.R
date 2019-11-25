rm(list=ls())

#libraries
library(tidyverse)
library(RColorBrewer)
library(devtools)
library(readxl)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
source_url("https://raw.githubusercontent.com/VFugere/LakePulse/master/custom_plots.R")

#### temp data ####

load('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/basic_data.RData')
kestrel <- read.csv2('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/2017 not merged/LakePulse2017_kestrel_QC.csv')
colnames(kestrel)[1] <- 'Lake_ID'
rbr <- read.csv2('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/2017 not merged/LakePulse2017_RBR_top_bottom_QC.csv')
colnames(rbr)[1] <- 'Lake_ID'
alex <- read.csv('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/2017 not merged/Alex_weather_station_data.csv')
colnames(alex)[1] <- 'Lake_ID'

kestrel <- kestrel %>% select(Lake_ID,temp) %>% rename(airtemp = temp)
rbr <- rbr %>% filter (depth == 'Tube lenght') %>% select(Lake_ID,temperature) %>% rename(watertemp = temperature)
alex <- alex %>% filter(Time >= 1950) %>% group_by(Lake_ID) %>% summarize(histtemp = mean(Mean_Ann_Temp, na.rm = T))
lat <- basic.data %>% select(Lake_ID, latitude)

merge <- inner_join(rbr,kestrel) %>% inner_join(alex) %>% inner_join(lat)
#corrgram::corrgram(merge)
plot(merge[,2:5])

#air temp and water temp correlate well on sampling day. So I can use air temp as proxy for water temp.
#mean annual air temp since 1950 correlate almost perfectly with latitude. So can use latitude as a proxy for air temp
#latitude should provide a proxy for water temp for all lakes.

#### plankton data ####

bad.zoo.samples <- c('07-057','17-050')#,'08-205','07-029') # remove these two if anything weird (see email Cindy 30-Sept-2019)

zoo <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton2017/all raw data 2017.xlsx', sheet='raw') %>%
  filter(!(ID_lakepulse %in% bad.zoo.samples)) %>%
  rename(Lake_ID = ID_lakepulse)



zoo.abund$name <- str_replace(zoo.abund$name, '  ', ' ')
zoo.abund$name <- str_replace(zoo.abund$name, 'spp', 'sp')


#phyto

phyto <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/phytoplankton2017/Phytoplankton LakePulse 2017_Jelena.xlsx', sheet='Longform')

phyto <- phyto %>% rename(Lake_ID = lake_id, bv = Biomass.mgm3, name = totalbinomial) %>%
  select(Lake_ID, name, bv) %>%
  filter(bv != 0) %>% 
  group_by(Lake_ID, name) %>% 
  summarize(bv = sum(bv, na.rm=T)) %>%
  ungroup %>%
  spread(name, bv) %>%
  as.data.frame

phyto[is.na(phyto)] <- 0

