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

d2019 <- read.csv('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/ALLfinal_grouping2019.csv', stringsAsFactors = F)

colnames(d2019) <- str_replace(colnames(d2019), '\\.\\.', ' ')
colnames(d2019) <- str_replace(colnames(d2019), 'spp\\.', 'sp\\.')
colnames(d2019) <- str_replace(colnames(d2019), 'sp\\.', 'sp')
colnames(d2019) <- str_replace(colnames(d2019), '\\.', ' ')

d2017 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/ALLfinal_grouping2017.csv', stringsAsFactors = F) %>%
  mutate_at(vars(Ergasilus.spp.:Simocephalus..spp.), as.numeric)

colnames(d2017) <- str_replace(colnames(d2017), '\\.\\.', ' ')
colnames(d2017) <- str_replace(colnames(d2017), 'spp\\.', 'sp\\.')
colnames(d2017) <- str_replace(colnames(d2017), 'sp\\.', 'sp')
colnames(d2017) <- str_replace(colnames(d2017), '\\.', ' ')

d2018 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/ALLfinal_grouping2018.csv', stringsAsFactors = F) %>%
  mutate_at(vars(Ergasilus.spp.:Simocephalus..spp.), as.numeric)

colnames(d2018) <- str_replace(colnames(d2018), '\\.\\.', ' ')
colnames(d2018) <- str_replace(colnames(d2018), 'spp\\.', 'sp\\.')
colnames(d2018) <- str_replace(colnames(d2018), 'sp\\.', 'sp')
colnames(d2018) <- str_replace(colnames(d2018), '\\.', ' ')

#bind all

zooLP <- bind_rows(d2017,d2018,d2019)
zooLP[is.na(zooLP)] <- 0
zooLP <- zooLP[,c(1,order(names(zooLP[,2:ncol(zooLP)]))+1)]

rm(d2017,d2018,d2019)

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

#map <- spTransform(map, CRS=CRS("+proj=aea +lat_1=26.605547828888216 +lat_2=62.451651032950174 +lon_0=-96.85546875"))
#plot(map)  

#https://projectionwizard.org/#
#+proj=aea +lat_1=33.99791731855025 +lat_2=63.00714506395516 +lon_0=-96.85546875

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

# only quebec

qc.lakes <- basic.data %>% filter(province == 'QUEBEC')

zoo <- filter(zooLP, Lake_ID %in% qc.lakes$Lake_ID)

invaded <- zoo %>% select(Lake_ID, `Eubosmina coregoni`) %>% filter(`Eubosmina coregoni` > 0)
merged <- inner_join(qc.lakes,invaded, by = 'Lake_ID')

col <- 'coral3'

#map

x <- qc.lakes$longitude
y <- qc.lakes$latitude
xrange <- range(x)+c(-2,2)
yrange <- range(y)+c(-1,1)
plot(map, xlim = xrange, ylim = yrange,col='light gray',border=0,asp=1,axes=F,cex.lab=0.5)
points(x=x,y=y,pch=1,col=1,cex=1.2)
points(x=merged$longitude,y=merged$latitude,pch=16,col=alpha(col,0.5),cex=1.1)

colSums(zoo[,2:100]) %>% sort

# number of lakes occupied

occupancy <- zoo[,2:ncol(zoo)]
occupancy[occupancy > 0] <- 1
occupancy <- select(zoo, Lake_ID) %>% bind_cols(occupancy)

sort((colSums(occupancy[,2:ncol(occupancy)])/nrow(occupancy))*100)

