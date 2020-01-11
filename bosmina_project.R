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

# 2019

#first do what Cindy did for other years, splitting biomass of un-Id'd species across species, and grouping copepodids with adults

dat <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/all raw data 2019.xlsx', sheet = 'clean biomass') 

# first group all copepodids id'd to species level with adults
taxo <- dat %>% filter(division == 'Copepoda') %>%
  distinct(`species name`) %>%
  filter(`species name` %!in% c("calanoid copepodid", "cyclopoid copepodid")) %>% 
  rename(name = `species name`) %>%
  arrange(name)
taxo$out <- taxo$name
taxo$out <- str_remove(taxo$out, pattern = ' copepodid')
taxo <- as.data.frame(taxo)
taxo$out[taxo$name == 'Heterocope septentrionalis'] <- 'Heterocope  septentrionalis'

dat[dat$`species name` %in% taxo$name,'species name'] <- taxo$out[match(dat$`species name`[dat$`species name` %in% taxo$name],taxo$name)]
colnames(dat)[2:4] <- c('name','division','biomass')
dat <- dat %>% group_by(ID_lakepulse,name,division) %>% summarize(biomass = sum(biomass)) %>% ungroup

#making a taxonomy data frame to distinguish calanoids from cyclopoids
taxo <- dat %>% filter(division == 'Copepoda') %>% distinct(name)
taxo$order <- 'calanoid'
taxo <- arrange(taxo, name)
taxo$order[12] <- 'poecilostomatoid'
#finding cyclopoids
taxo$order[c(1,5,6,7,8,13,14,15,28,29,30,31,33,36)] <- 'cyclopoid'

dat$genus <- str_remove(dat$name, '\\ .*')
dat$cyc_order <- taxo$order[match(dat$name,taxo$name)]

writexl::write_xlsx(dat, '~/Desktop/2019.xlsx')

out.dat <- data.frame

lakes <- unique(dat$ID_lakepulse)

i<-1

sub <- dat %>% filter(ID_lakepulse == lakes[i])

#split copepodids based on species, as long as one species was id'd to genus or species
sub$copepopid <- 'no'
sub$copepopid[str_detect(sub$name, 'copepodid')] <- 'yes'
sub.adults <- filter(sub, copepopid == 'no', !is.na(cyc_order))
sub.copepodid <- filter(sub, copepopid == 'yes')
if(nrow(sub.adults)>0){
  for(c in unique(sub.copepodid$cyc_order)){
    order <- c
    babybiomass <- sub.copepodid %>% filter(cyc_order == c) %>% pull(biomass)
    prop.dat <- sub.adults %>% filter(cyc_order == c)
    prop.dat$props <- prop.dat$biomass/sum(prop.dat$biomass)
    prop.dat$bm.to.add <- prop.dat$props*babybiomass
    prop.dat$final.bm <- prop.dat$biomass+prop.dat$bm.to.add
    sub$biomass[match(prop.dat$name,sub$name)] <- prop.dat$final.bm
  }
  sub <- filter(sub, copepopid == 'no')
}  
#### need to add scenario where one order has an adult id, and one doesnt, but there are copepodids of both orders


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
