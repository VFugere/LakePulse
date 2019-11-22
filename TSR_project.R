rm(list=ls())

#libraries
library(tidyverse)
library(RColorBrewer)
library(devtools)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
source_url("https://raw.githubusercontent.com/VFugere/LakePulse/master/custom_plots.R")

#cols
cols <- brewer.pal(7, 'Dark2')[1:3] 
cols2 <- brewer.pal(8, 'Dark2')[c(4,5,6,8)] 

#data
load('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/basic_data.RData')

#
#zoo Cindy

bad.zoo.samples <- c('07-057','17-050')#,'08-205','07-029') # remove these two if anything weird (see email Cindy 30-Sept-2019)

zoo.abund <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton2017/all raw data 2017.xlsx', sheet='abundances') %>%
  filter(!(ID_lakepulse %in% bad.zoo.samples)) %>%
  rename(Lake_ID = ID_lakepulse, abund = '#individuals counted') %>%
  select(Lake_ID, name, abund)

zoo.abund$name <- str_replace(zoo.abund$name, '  ', ' ')
zoo.abund$name <- str_replace(zoo.abund$name, 'spp', 'sp')

zoo.abund <- zoo.abund %>% 
  group_by(Lake_ID,name) %>% 
  summarize(abund = sum(abund, na.rm=T)) %>%
  ungroup %>%
  spread(name,abund)

#

zoo.biomass <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton2017/all raw data 2017.xlsx', sheet='clean biomass') %>%
  filter(!(ID_lakepulse %in% bad.zoo.samples)) %>%
  rename(Lake_ID = ID_lakepulse, biomass = 'species biomass (µg d.w./L)') %>%
  select(Lake_ID, name, biomass)

zoo.biomass$name <- str_replace(zoo.biomass$name, '  ', ' ')
zoo.biomass$name <- str_replace(zoo.biomass$name, 'spp', 'sp')

zoo.biomass <- spread(zoo.biomass,name, biomass) #ug DM L^-1, copepodite not grouped with adults

#

zoo.biomass.grouped <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton2017/ALLfinal_grouping2017.csv', stringsAsFactors = F) %>%
  filter(!(Lake_ID %in% bad.zoo.samples)) %>%
  mutate_at(vars(Ergasilus.spp.:Simocephalus..spp.), as.numeric)

colnames(zoo.biomass.grouped) <- str_replace(colnames(zoo.biomass.grouped), '\\.', ' ')
colnames(zoo.biomass.grouped) <- str_replace(colnames(zoo.biomass.grouped), '\\.spp\\.', 'sp\\.')
colnames(zoo.biomass.grouped) <- str_replace(colnames(zoo.biomass.grouped), 'spp\\.', 'sp\\.')
colnames(zoo.biomass.grouped) <- str_replace(colnames(zoo.biomass.grouped), 'Daphnia galeata\\.mendotae', 'Daphnia galeata mendotae')
colnames(zoo.biomass.grouped) <- str_replace(colnames(zoo.biomass.grouped), ' \\.', ' ')

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

