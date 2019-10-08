rm(list=ls())

library(tidyverse)
library(readxl)

#basic data data

d2017 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2017_basic_info.csv', sep = ';', stringsAsFactors = F)
d2017$year <- 2017
d2018 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2018_basic_info.csv', sep = ';', stringsAsFactors = F)
d2018$year <- 2018
d2019 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2019_basic_info.csv', sep = ';', stringsAsFactors = F)
d2019$year <- 2019

basic.data <- bind_rows(d2017,d2018) %>%
  bind_rows(d2019) %>%
  rename(Lake_ID = lakepulse_id, area = size_km2, HI = hi_index)

rm(d2017,d2018,d2019)
  
#MP's zoo trait database

comma2dot <- function(x){as.numeric(str_replace(x,',','.'))}
return1stval <- function(x){x[1]}

zooT <- read.csv('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/traits/zootraits.csv', sep=';', stringsAsFactors = F)
zooT <- zooT %>%
  mutate_at(c('Body.length','Min.bl','Max.bl','Dry.mass','Min.dm','Max.dm'), comma2dot) %>% 
  filter(Habitat == 'Freshwater') %>%
  select(-Ref.tg, -Ref.bl, -Habitat, -Ref.dm) %>%
  select(Genus:Max.dm)
                  
zooT$taxon <- paste(zooT$Genus,zooT$Species,sep=' ')
colnames(zooT)[c(3,6:11)] <- c('rep','length','min.length','max.length','dry.mass','min.dry.mass','max.dry.mass')

zooT <- select(zooT, taxon, everything())

zooT <- zooT %>% group_by(taxon) %>% 
  summarize(genus = return1stval(Genus),
            species = return1stval(Species),
            class = return1stval(Group),
            TL = return1stval(Trophic.group),
            length = mean(length, na.rm=T),
            max.length = max(max.length, na.rm=T),
            min.length = min(min.length, na.rm=T),
            dry.mass = mean(dry.mass, na.rm=T),
            max.dry.mass = max(max.dry.mass, na.rm=T),
            min.dry.mass = min(min.dry.mass, na.rm=T)) %>%
  as.data.frame()

for(i in 6:11){
  zooT[!(is.finite(zooT[,i])),i] <- NA
}

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

#phyto Bruno

phyto <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/phytoplankton2017/Phytoplankton LakePulse 2017_Bruno.xlsx', sheet='Longform')

phyto <- phyto %>% rename(Lake_ID = lake_id, bv = BV.mm3L.corrected, name = totalbinomial) %>%
  select(Lake_ID, name, bv) %>%
  group_by(Lake_ID, name) %>% 
  summarize(bv = sum(bv, na.rm=T)) %>%
  ungroup %>%
  spread(name, bv) %>%
  filter(Lake_ID != '06-174') %>% #empty lake! biovolume of 0
  as.data.frame

phyto[is.na(phyto)] <- 0

# colnames(phyto) <- str_replace(colnames(phyto), '\\(', '')
# colnames(phyto) <- str_replace(colnames(phyto), '\\)', '')

# phyto traits and taxonomy

phytoT <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/phytoplankton2017/Phytoplankton LakePulse 2017_Bruno.xlsx', sheet='Unique taxa entries')

# ASV data Susanne

bacterio <- read.table('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/bacterioplankton2017/otu_table_ASV_level.txt', header=T, stringsAsFactors = F)
colnames(bacterio)[1] <- 'Lake_ID'

# Mario M. recommended deleting global singleton but not global doubleton or more. I.e. no need for rarefaction
sum(colSums(bacterio[,2:ncol(bacterio)]) == 0)
sum(is.na(colSums(bacterio[,2:ncol(bacterio)])))
#no NA, no empty columns. 
sum(colSums(bacterio[,2:ncol(bacterio)]) == 1)
sum(colSums(bacterio[,2:ncol(bacterio)]) == 2)
#5 global singletons, 557 global doubletons

#removing two sites with very poor sequencing depth (< 1100 reads, while all other sites are > 10K)
to.rm <- which(rowSums(bacterio[,2:ncol(bacterio)]) < 1500)
bacterio <- bacterio[-(to.rm),]

#removing global single and doubletons
to.rm <- as.numeric(which(colSums(bacterio[,2:ncol(bacterio)]) < 3))+1
bacterio <- bacterio[,-(to.rm)]

#bacterio taxonomy

bacterioT <- read.table('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/bacterioplankton2017/taxa_ASV.txt', header=T, stringsAsFactors = F)

#check that all names are ok and found in metadata file
zoo.abund$Lake_ID[!(zoo.abund$Lake_ID %in% basic.data$Lake_ID)]
zoo.biomass$Lake_ID[!(zoo.biomass$Lake_ID %in% basic.data$Lake_ID)]
zoo.biomass.grouped$Lake_ID[!(zoo.biomass.grouped$Lake_ID %in% basic.data$Lake_ID)]
phyto$Lake_ID[!(phyto$Lake_ID %in% basic.data$Lake_ID)]
bacterio$Lake_ID[!(bacterio$Lake_ID %in% basic.data$Lake_ID)]

# saving

rm(bad.zoo.samples,i,to.rm,comma2dot,return1stval)
save.image('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/2017data.RData')
