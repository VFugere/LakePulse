rm(list=ls())

#libraries
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(rworldmap)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
devtools::source_url("https://raw.githubusercontent.com/VFugere/LakePulse/master/custom_plots.R")

#cols
cols <- brewer.pal(3, 'Dark2')
cols2 <- brewer.pal(8, 'Dark2')[c(4,5,6,8)] 

#data
load('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/2017data.RData')
zoo <- zoo.biomass.grouped 
basic.data <- filter(basic.data, year == 2017)
basic.data$ecozone <- as.factor(basic.data$ecozone)
levels(basic.data$ecozone) <- c('AH','AM','BS','MP')

#map
map <- getMap(resolution = "low")

#land use data
lulc <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2017_LULC_QC.csv', sep = ';', stringsAsFactors = F)
lulc <- select(lulc, -comment,-flag)
basic.data <- left_join(basic.data,lulc, by = c('Lake_ID' = 'lakepulse_id'))
rm(lulc)

#### Creating a merged dataset with diversity of prokaryotic, eukaryotic phyto, and eukaryotic zoo

hist(log(apply(phyto[,2:ncol(phyto)],1,FUN='sum')),breaks=50)
hist(log(apply(zoo[,2:ncol(zoo)],1,FUN='sum')),breaks=50)

phyto.long <- gather(phyto, taxon, biov ,-Lake_ID)
phyto.long$class <- phytoT$KINDGOM[match(phyto.long$taxon, phytoT$totalbinomial)]
  
phyto.euk <- filter(phyto.long, class != 'CYANOBACTERIA') %>%
  select(-class) %>%
  spread(taxon, biov)

b.div <- data.frame('Lake_ID' = bacterio$Lake_ID, 'bdiv' = specnumber(bacterio[,2:ncol(bacterio)]))
p.div <- data.frame('Lake_ID' = phyto.euk$Lake_ID, 'pdiv' = specnumber(phyto.euk[,2:ncol(phyto.euk)]))
z.div <- data.frame('Lake_ID' = zoo$Lake_ID, 'zdiv' = specnumber(zoo[,2:ncol(zoo)]))

prop.func <- function(x){x/.data$tdiv}
div.merge <- inner_join(b.div, p.div) %>% inner_join(z.div) %>%
  mutate(tdiv = bdiv+pdiv+zdiv) %>%
  left_join(basic.data)
div.merge$brel <- with(div.merge, bdiv/tdiv)
div.merge$prel <- with(div.merge, pdiv/tdiv)
div.merge$zrel <- with(div.merge, zdiv/tdiv)
div.merge$pz.ratio <- with(div.merge, zdiv/pdiv)

mybubble(div.merge$area,div.merge$HI,div.merge$pdiv,name='phyto div',ez=div.merge$ecozone)
mybubble(div.merge$area,div.merge$HI,div.merge$bdiv,name='bacterio div',ez=div.merge$ecozone)
mybubble(div.merge$area,div.merge$HI,div.merge$zdiv,name='zoo div',ez=div.merge$ecozone)
mybubble(div.merge$area,div.merge$HI,div.merge$prel,name='rel phyto div',ez=div.merge$ecozone)
mybubble(div.merge$area,div.merge$HI,div.merge$brel,name='rel bacterio div',ez=div.merge$ecozone)
mybubble(div.merge$area,div.merge$HI,div.merge$zrel,name='rel zoo div',ez=div.merge$ecozone)
mybubble(div.merge$area,div.merge$HI,div.merge$pz.ratio,name='zoo/phyto div',ez=div.merge$ecozone)

mapplot(div.merge$longitude,div.merge$latitude,log(div.merge$bdiv),'bacterial diversity (logged)')
mapplot(div.merge$longitude,div.merge$latitude,div.merge$pdiv,'phytoplankton diversity')
mapplot(div.merge$longitude,div.merge$latitude,div.merge$zdiv,'zooplankton diversity')
mapplot(div.merge$longitude,div.merge$latitude,div.merge$prel,'phytoplankton relative diversity')
mapplot(div.merge$longitude,div.merge$latitude,div.merge$zrel,'zooplankton relative diversity')
mapplot(div.merge$longitude,div.merge$latitude,div.merge$pz.ratio,'zooplankton:phytoplankton diversity')

div.merge %>% select(depth_m,watershed_km2,latitude,longitude,fraction_agriculture:fraction_water,pdiv,bdiv,zdiv,prel,brel,zrel,pz.ratio) %>% plot
div.merge %>% select(depth_m,pdiv,bdiv,zdiv,prel,brel,zrel,pz.ratio) %>% mutate_all(log) %>% plot

#### same thing with biomass/biovolume ####

ep.biov <- data.frame('Lake_ID' = phyto.euk$Lake_ID, 'epbiov' = apply(phyto.euk[,2:ncol(phyto.euk)],1,FUN='sum'))
p.biov <- data.frame('Lake_ID' = phyto$Lake_ID, 'pbiov' = apply(phyto[,2:ncol(phyto)],1,FUN='sum'))
z.biom <- data.frame('Lake_ID' = zoo$Lake_ID, 'zbiom' = apply(zoo[,2:ncol(zoo)],1,FUN='sum'))

bio.merge <- inner_join(p.biov, z.biom) %>%
  left_join(basic.data) %>%
  mutate(prop = zbiom/pbiov)

mybubble(bio.merge$area,bio.merge$HI,log(bio.merge$pbiov))
mybubble(bio.merge$area,bio.merge$HI,bio.merge$zbiom)
mybubble(bio.merge$area,bio.merge$HI,bio.merge$prop)

mapplot(bio.merge$longitude,bio.merge$latitude,log(bio.merge$pbiov),'phyto biovolume (logged)')
mapplot(bio.merge$longitude,bio.merge$latitude,log(bio.merge$zbiom),'zoo biomass (logged)')
mapplot(bio.merge$longitude,bio.merge$latitude,bio.merge$prop,'zoo:phyto')

bio.merge %>% select(depth_m,pbiov,zbiom,prop) %>% plot
bio.merge %>% select(depth_m,pbiov,zbiom,prop) %>% mutate_all(log) %>% plot

#

bio.merge <- inner_join(ep.biov, z.biom) %>%
  left_join(basic.data) %>%
  mutate(prop = zbiom/epbiov)

mybubble(bio.merge$area,bio.merge$HI,log(bio.merge$epbiov))
mybubble(bio.merge$area,bio.merge$HI,bio.merge$zbiom)
mybubble(bio.merge$area,bio.merge$HI,bio.merge$prop,'zoo/phyto',ez=bio.merge$ecozone)

mapplot(bio.merge$longitude,bio.merge$latitude,log(bio.merge$epbiov),'phyto biovolume (logged)')
mapplot(bio.merge$longitude,bio.merge$latitude,log(bio.merge$zbiom),'zoo biomass (logged)')
mapplot(bio.merge$longitude,bio.merge$latitude,log(bio.merge$prop),'zoo:phyto')

bio.merge %>% select(depth_m,epbiov,zbiom,prop) %>% plot
bio.merge %>% select(depth_m,epbiov,zbiom,prop) %>% mutate_all(log) %>% plot
