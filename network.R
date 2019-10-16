rm(list=ls())

#libraries
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(scales)
library(netassoc)


?make_netassoc_network

https://cran.r-project.org/web/packages/netassoc/index.html

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
devtools::source_url("https://raw.githubusercontent.com/VFugere/LakePulse/master/custom_plots.R")

#cols
cols <- brewer.pal(3, 'Dark2')
cols2 <- brewer.pal(8, 'Dark2')[c(4,5,6,8)] 

#data
load('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/2017data.RData')
basic.data <- filter(basic.data, year == 2017)
basic.data$ecozone <- as.factor(basic.data$ecozone)
levels(basic.data$ecozone) <- c('AH','AM','BS','MP')

#land use data
lulc <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2017_LULC_QC.csv', sep = ';', stringsAsFactors = F)
lulc <- select(lulc, -comment,-flag)
basic.data <- left_join(basic.data,lulc, by = c('Lake_ID' = 'lakepulse_id'))
rm(lulc)

## formatting and adding some taxonomic info

phyto.long <- gather(phyto, taxon, biov ,-Lake_ID)
phyto.long$class <- phytoT$KINDGOM[match(phyto.long$taxon, phytoT$totalbinomial)]

phyto.euk <- filter(phyto.long, class != 'CYANOBACTERIA') %>%
  select(-class) %>%
  spread(taxon, biov)

zoo <- zoo.biomass.grouped 
zooT <- readxl::read_xlsx('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/traits/zoo_trait_summary.xlsx')

zoo.long <- gather(zoo.biomass.grouped, taxon, biom ,-Lake_ID)
zoo.long$genus <- sub('\\s.*', '', zoo.long$taxon)
zoo.long$class <- zooT$class[match(zoo.long$taxon,zooT$taxon)]
zoo.long$TG <- zooT$TG[match(zoo.long$taxon,zooT$taxon)]
zoo.long$TGS <- zooT$TGS[match(zoo.long$taxon,zooT$taxon)]

zoo.herb <- filter(zoo.long, TG %in% c('herb','omni')) %>%
  select(-TG) %>%
  spread(taxon, biom)

# dataframe for SEM

data <- inner_join(zoo.herb,)