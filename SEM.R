rm(list=ls())

#libraries
library(tidyverse)
library(RColorBrewer)
library(performance)
library(nlme)
library(lme4)
library(lavaan)
library(piecewiseSEM)
library(scales)

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

phyto.ffg <- phyto.long %>% group_by(Lake_ID, class) %>%
  summarize(biov = sum(biov)) %>%
  ungroup %>%
  spread(class, biov)
colnames(phyto.ffg)[2:ncol(phyto.ffg)] <- str_to_lower(colnames(phyto.ffg)[2:ncol(phyto.ffg)])
phyto.ffg$others <- with(phyto.ffg, chrysophyceae+cryptophyceae+dinophyceae+euglenophyceae+haptophyte)

zoo <- zoo.biomass.grouped 
zooT <- readxl::read_xlsx('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/traits/zoo_trait_summary.xlsx')

zoo.long <- gather(zoo.biomass.grouped, taxon, biom ,-Lake_ID)
zoo.long$genus <- sub('\\s.*', '', zoo.long$taxon)
zoo.long$class <- zooT$class[match(zoo.long$taxon,zooT$taxon)]
zoo.long$TG <- zooT$TG[match(zoo.long$taxon,zooT$taxon)]
zoo.long$TGS <- zooT$TGS[match(zoo.long$taxon,zooT$taxon)]

zoo.ffg <- zoo.long %>% group_by(Lake_ID, TGS) %>%
  summarize(biom = sum(biom)) %>%
  ungroup %>%
  spread(TGS, biom)

# dataframe for SEM

data <- inner_join(phyto.ffg, zoo.ffg) %>%
  left_join(basic.data)

data$omni <- with(data, omni+pred)

data <- data %>% rename('greens' = chlorophyceae, 'cyanos' = cyanobacteria, 'depth' = depth_m) %>%
  select(greens,diatoms,cyanos,others,small_herb,large_herb,omni,area,depth,HI,ecozone,HI,fraction_agriculture:fraction_urban)

data$euk.phyto <- with(data, greens+diatoms+others)

data <- mutate_at(data, vars(greens:others,area,depth,euk.phyto), log)
data <- mutate_at(data, vars(small_herb,large_herb,omni,HI,fraction_agriculture:fraction_urban), log1p)

data$greens[which(data$greens == -Inf)] <- NA
data$cyanos[which(data$cyanos == -Inf)] <- NA
data$others[which(data$others == -Inf)] <- NA
data$diatoms[which(data$diatoms == -Inf)] <- NA

dat <- select(data, greens:omni, euk.phyto, depth, HI, ecozone) %>% as.data.frame
dat <- na.omit(dat)
dat <- mutate_if(dat, is.numeric, scale)

mod1 <- psem(
  #cyanos
  lme(cyanos ~ depth + HI, random = ~ 1|ecozone, data = dat),
  #greens
  lme(greens ~ depth + HI, random = ~ 1|ecozone, data = dat),
  #diatoms
  lme(diatoms ~ depth + HI, random = ~ 1|ecozone, data = dat),
  #others
  lme(others ~ depth + HI, random = ~ 1|ecozone, data = dat),
  #small herb
  lme(small_herb ~ depth + HI + cyanos + greens + diatoms + others, random = ~ 1|ecozone, data = dat),
  #large herb
  lme(large_herb ~ depth + HI + cyanos + greens + diatoms + others, random = ~ 1|ecozone, data = dat),
  #omni
  lme(omni ~ depth + HI + cyanos + greens + diatoms + others + small_herb, random = ~ 1|ecozone, data = dat),
  data = dat)

summary(mod1,conserve=T)

mod2 <- psem(
  #cyanos
  lm(cyanos ~ depth + HI, data = dat),
  #greens
  lm(greens ~ depth + HI + cyanos, data = dat),
  #diatoms
  lm(diatoms ~ depth + HI + cyanos, data = dat),
  #others
  lm(others ~ depth + HI, data = dat),
  #herb
  lm(herb ~ depth + HI + cyanos + greens + diatoms + others, data = dat),
  #omni
  lm(omni ~ depth + HI + cyanos + greens + diatoms + others + herb, data = dat)
)
summary(mod2)

mod3 <- psem(

  lm(cyanos ~ depth + HI, data = dat),

  lm(diatoms ~ depth + HI, data = dat),
  lm(greens ~ depth + HI, data = dat),
  lm(others ~ depth + HI, data = dat),

  lm(small_herb ~ depth + HI + cyanos + diatoms + greens + others, data = dat),

  lm(large_herb ~ depth + HI + diatoms + greens + others + small_herb, data = dat),

  lm(omni ~ depth + HI + diatoms + greens + others + small_herb, data = dat),
  
  cyanos %~~% diatoms,
  cyanos %~~% greens,
  cyanos %~~% others,
  diatoms %~~% greens,
  others %~~% greens,
  others %~~% diatoms,
  large_herb %~~% omni
  
)
summary(mod3)

coefs$est.sc <- rescale(coefs$Std.Estimate, to=c(0.5,5))
coefs %>% select(Predictor, Response, Estimate, P.Value, est.sc) %>% arrange(Predictor, P.Value)
