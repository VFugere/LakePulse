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

d2017 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton2017/ALLfinal_grouping2017.csv', stringsAsFactors = F) %>%
  mutate_at(vars(Ergasilus.spp.:Simocephalus..spp.), as.numeric)

colnames(d2017) <- str_replace(colnames(d2017), '\\.', ' ')
colnames(d2017) <- str_replace(colnames(d2017), '\\.spp\\.', 'sp\\.')
colnames(d2017) <- str_replace(colnames(d2017), 'spp\\.', 'sp\\.')
colnames(d2017) <- str_replace(colnames(d2017), 'Daphnia galeata\\.mendotae', 'Daphnia galeata mendotae')
colnames(d2017) <- str_replace(colnames(d2017), ' \\.', ' ')

colSums(d2017[,2:69]) %>% sort

# 2018

d2018 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton2018/ALLfinal_grouping2018.csv', stringsAsFactors = F) %>%
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

#make a map

colSums(zooLP[,2:100]) %>% sort
