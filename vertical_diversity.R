rm(list=ls())

#libraries
library(tidyverse)
library(vegan)
library(RColorBrewer)
#library(party)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))

#cols
cols <- brewer.pal(3, 'Dark2')

#data
load('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/2017data.RData')

plot(HumanImpac~Lake.area_km2,env.data) #no large impacted lakes. 
# Does not make sense to look at variation within groups

#richness across gradients
rich <- merged.data %>%
  mutate(rich.z = specnumber())
  