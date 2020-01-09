### LOAD LAKE PULSE 2017 PLANKTON DATASETS
### ADD SOME TRAITS AND OTHER INFO

# code by Vincent Fugère (2019)

rm(list=ls())

library(tidyverse)
library(readxl)
library(vegan)

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

zoo.abund <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/all raw data 2017.xlsx', sheet='abundances') %>%
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

zoo.biomass <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/all raw data 2017.xlsx', sheet='clean biomass') %>%
  filter(!(ID_lakepulse %in% bad.zoo.samples)) %>%
  rename(Lake_ID = ID_lakepulse, biomass = 'species biomass (µg d.w./L)') %>%
  select(Lake_ID, name, biomass)

zoo.biomass$name <- str_replace(zoo.biomass$name, '  ', ' ')
zoo.biomass$name <- str_replace(zoo.biomass$name, 'spp', 'sp')

zoo.biomass <- spread(zoo.biomass,name, biomass) #ug DM L^-1, copepodite not grouped with adults

#

zoo.biomass.grouped <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/ALLfinal_grouping2017.csv', stringsAsFactors = F) %>%
  filter(!(Lake_ID %in% bad.zoo.samples)) %>%
  mutate_at(vars(Ergasilus.spp.:Simocephalus..spp.), as.numeric)

colnames(zoo.biomass.grouped) <- str_replace(colnames(zoo.biomass.grouped), '\\.', ' ')
colnames(zoo.biomass.grouped) <- str_replace(colnames(zoo.biomass.grouped), '\\.spp\\.', 'sp\\.')
colnames(zoo.biomass.grouped) <- str_replace(colnames(zoo.biomass.grouped), 'spp\\.', 'sp\\.')
colnames(zoo.biomass.grouped) <- str_replace(colnames(zoo.biomass.grouped), 'Daphnia galeata\\.mendotae', 'Daphnia galeata mendotae')
colnames(zoo.biomass.grouped) <- str_replace(colnames(zoo.biomass.grouped), ' \\.', ' ')

#phyto

phyto <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/phytoplankton/phyto2017_clean.xlsx', sheet='Longform')

# phyto traits and taxonomy
phytoT <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/phytoplankton/phyto2017_clean.xlsx', sheet='taxonomy_clean')
phytoT$group <- tolower(phytoT$group)

#use clean taxonomy sheet to replace erroneous species name in community matrix
phyto$totalbinomial <- phytoT$clean.name[match(phyto$totalbinomial,phytoT$totalbinomial)]

#format into community matrix
phyto <- phyto %>% rename(Lake_ID = lake_id, bv = Biomass.mgm3, name = totalbinomial) %>%
  select(Lake_ID, name, bv) %>%
  filter(bv != 0) %>% 
  group_by(Lake_ID, name) %>% 
  summarize(bv = sum(bv, na.rm=T)) %>%
  ungroup %>%
  spread(name, bv) %>%
  as.data.frame

phyto[is.na(phyto)] <- 0

# colnames(phyto) <- str_replace(colnames(phyto), '\\(', '')
# colnames(phyto) <- str_replace(colnames(phyto), '\\)', '')

# ASV data Susanne

bacterio <- read.table('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/bacterioplankton/otu_table_ASV_level.txt', header=T, stringsAsFactors = F)
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

#rarefying to 15K to be consistent with Susanne
bacterio[,2:ncol(bacterio)] <- rrarefy(bacterio[,2:ncol(bacterio)], 15000)
                                       
#removing global single and doubletons
to.rm <- as.numeric(which(colSums(bacterio[,2:ncol(bacterio)]) < 3))+1
bacterio <- bacterio[,-(to.rm)]

#bacterio taxonomy

bacterioT <- read.table('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/bacterioplankton/taxa_ASV.txt', header=T, stringsAsFactors = F)

# saving

rm(bad.zoo.samples,i,to.rm,comma2dot,return1stval)
save.image('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/2017data.RData')
