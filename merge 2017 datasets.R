rm(list=ls())

library(tidyverse)
library(readxl)

#env data

env.data <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/env_all_sites_May2018.csv', sep = ';', stringsAsFactors = F)

# #MP's zoo trait database
# zooT <- read.csv('~/Desktop/zooTtraits.csv', sep=';', stringsAsFactors = F)
# zooT$Body.length <- str_replace(zooT$Body.length, ',','.')
# zooT$Body.length <- as.numeric(zooT$Body.length)
# 
# zooT <- zooT %>%
#   filter(Habitat == 'Freshwater') %>%
#   select(Genus:Body.length) %>%
#   rename(length = Body.length) %>%
#   select(-Ref.tg, -Habitat)
# 
# zooT$sp <- paste(zooT$Genus,zooT$Species,sep=' ')

#zoo Cindy

bad.zoo.samples <- c('07-057','17-050')#,'08-205','07-029') # remove these two if anything weird (see email Cindy 30-Sept-2019)

zoo.abund <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton2017/all raw data 2017.xlsx', sheet='abundances') %>%
  filter(!(ID_lakepulse %in% bad.zoo.samples)) %>%
  rename(Lake_ID = ID_lakepulse, abund = '#individuals counted') %>%
  mutate(species = paste(genus,species)) %>%
  select(Lake_ID, species, abund) %>%
  group_by(Lake_ID, species) %>% 
  summarize(abund = sum(abund)) %>%
  ungroup %>%
  spread(species, abund)

zoo.biomass <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton2017/all raw data 2017.xlsx', sheet='clean biomass') %>%
  filter(!(ID_lakepulse %in% bad.zoo.samples)) %>%
  rename(Lake_ID = ID_lakepulse, biomass = 'species biomass (µg d.w./L)') %>%
  select(Lake_ID, name, biomass) %>%
  spread(name, biomass) #ug DM L^-1, copepodite not grouped with adults

#phyto Bruno

phyto.biov <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/phytoplankton2017/Phytoplankton LakePulse 2017_Bruno.xlsx', sheet='Longform')

phyto.biov <- phyto.biov %>% rename(Lake_ID = lake_id, bv = BV.mm3L.corrected, name = totalbinomial) %>%
  select(Lake_ID, name, bv) %>%
  group_by(Lake_ID, name) %>% 
  summarize(bv = sum(bv, na.rm=T)) %>%
  ungroup %>%
  spread(name, bv) %>%
  as.data.frame

phyto.biov[is.na(phyto.biov)] <- 0

colnames(phyto.biov)
str(phyto.biov)

# ASV data Susanne

bacterio <- read.table('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/bacterioplankton2017/otu_table_ASV_level.txt', header=T, stringsAsFactors = F)
colnames(bacterio)[1] <- 'Lake_ID'

# Mario M. recommended deleting global singleton but not global doubleton or more. I.e. no need for rarefaction
sum(colSums(bacterio[,2:13435]) == 0)
sum(is.na(colSums(bacterio[,2:13435])))
#no NA, no empty columns. 
sum(colSums(bacterio[,2:13435]) == 1)
sum(colSums(bacterio[,2:13435]) == 2)
#5 global singletons, 557 global doubletons

# master dataframe

#check that all names are ok and found in metadata file
zoo.biomass$Lake_ID[!(zoo.biomass$Lake_ID %in% env.data$Lake_ID)]
phyto.biov$Lake_ID[!(phyto.biov$Lake_ID %in% env.data$Lake_ID)]
bacterio$Lake_ID[!(bacterio$Lake_ID %in% env.data$Lake_ID)]

data <- inner_join(env.data, zoo.biomass, by = 'Lake_ID')
data <- inner_join(data, phyto.biov, by = 'Lake_ID') #we loose 13 lakes

# in phyto but not zoo database
phyto.biov$Lake_ID[!(phyto.biov$Lake_ID %in% zoo.biomass$Lake_ID)] #indeed, 13 lakes

#in zoo database but not phyto database
zoo.biomass$Lake_ID[!(zoo.biomass$Lake_ID %in% phyto.biov$Lake_ID)]

merged.data <- inner_join(data, bacterio, by = 'Lake_ID') #we loose 4 additional lakes

rm(data,bad.zoo.samples)
save.image('~/Desktop/2017data.RData')
