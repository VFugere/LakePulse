# Analysis script to track invasive species Bosmina Coregoni

rm(list=ls())

#libraries
library(tidyverse)
library(readxl)
library(vegan)
library(RColorBrewer)
library(rworldmap)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
devtools::source_url("https://raw.githubusercontent.com/VFugere/LakePulse/master/custom_plots.R")

## Load and format zoo data

# 2017

d2017 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/ALLfinal_grouping2017.csv', stringsAsFactors = F) %>%
  mutate_at(vars(Ergasilus.spp.:Simocephalus..spp.), as.numeric)

colnames(d2017) <- str_replace(colnames(d2017), '\\.', ' ')
colnames(d2017) <- str_replace(colnames(d2017), '\\.spp\\.', 'sp\\.')
colnames(d2017) <- str_replace(colnames(d2017), 'spp\\.', 'sp\\.')
colnames(d2017) <- str_replace(colnames(d2017), 'Daphnia galeata\\.mendotae', 'Daphnia galeata mendotae')
colnames(d2017) <- str_replace(colnames(d2017), ' \\.', ' ')

colSums(d2017[,2:69]) %>% sort

# 2018

d2018 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/ALLfinal_grouping2018.csv', stringsAsFactors = F) %>%
  mutate_at(vars(Ergasilus.spp.:Simocephalus..spp.), as.numeric)

colnames(d2018) <- str_replace(colnames(d2018), '\\.', ' ')
colnames(d2018) <- str_replace(colnames(d2018), '\\.spp\\.', 'sp\\.')
colnames(d2018) <- str_replace(colnames(d2018), 'spp\\.', 'sp\\.')
colnames(d2018) <- str_replace(colnames(d2018), 'Daphnia galeata\\.mendotae', 'Daphnia galeata mendotae')
colnames(d2018) <- str_replace(colnames(d2018), ' \\.', ' ')

colSums(d2018[,2:93]) %>% sort

#bind all

zooLP <- bind_rows(d2017,d2018)
zooLP[is.na(zooLP)] <- 0
zooLP <- zooLP[,c(1,order(names(zooLP[,2:ncol(zooLP)]))+1)]

# species total biomass

colSums(zooLP[,2:100]) %>% sort

# number of lakes occupied

occupancy <- zooLP[,2:ncol(zooLP)]
occupancy[occupancy > 0] <- 1
occupancy <- select(zooLP, Lake_ID) %>% bind_cols(occupancy)

sort((colSums(occupancy[,2:ncol(occupancy)])/nrow(occupancy))*100)

#make a map

library(rgdal)
map <- getMap(resolution = "low")

map <- getMap(resolution = 'low')[which(getMap()$ADMIN %in% c('Canada','United States of America')),]
map <- spTransform(map, CRS=CRS("+proj=aea +lat_1=26.605547828888216 +lat_2=62.451651032950174 +lon_0=-96.85546875"))
plot(map)  

load('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/basic_data.RData')

invaded <- zooLP %>% select(Lake_ID, `Eubosmina coregoni`) %>% filter(`Eubosmina coregoni` > 0)
merged <- inner_join(basic.data,invaded, by = 'Lake_ID')

col <- 'coral3'

x <- basic.data$longitude
y <- basic.data$latitude
xrange <- range(x)+c(-2,2)
yrange <- range(y)+c(-1,1)
plot(map, xlim = xrange, ylim = yrange,col='light gray',border=0,asp=1,axes=F,cex.lab=0.5)
points(x=x,y=y,pch=1,col=1,cex=1.2)
points(x=merged$longitude,y=merged$latitude,pch=16,col=alpha(col,0.5),cex=1.1)
