# Analysis script to model max trophic level of fish communities in Lake Pulse sites,
# using either only abiotic parameters or using also zooplankton. Does including zoopplankton
# increase predictive power?

rm(list=ls())

#libraries
library(tidyverse)
library(readxl)
library(vegan)
library(devtools)
library(party)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
source_url("https://raw.githubusercontent.com/VFugere/LakePulse/master/custom_plots.R")
mean.nona <- function(x){mean(x, na.rm=T)}

### load data ###

## lake pulse abiotic data

load('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/basic_data.RData')
basic.data <- basic.data %>% arrange(Lake_ID) %>% select(Lake_ID,latitude,longitude,area,depth_m,ecozone,Shore_dev,Res_time,feow,cont.watershed)
altitude <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_altitude.xlsx') %>% rename(Lake_ID = lakepulse_id) %>% arrange(Lake_ID)
basic.data <- left_join(basic.data, altitude)
rm(altitude)
basic.data <- basic.data %>% mutate_at(vars(ecozone,cont.watershed), as.factor)

lulc.watershed <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LULC_QC20200323.xlsx') %>% arrange(idLakePulse)
#NOTE THAT WATERSHED AREA AND HII ARE NOT EXACTLY EQUAL IN BASIC AND LULC DATASETS. LULC data more up-to-date so I use this instead.
colnames(lulc.watershed )[c(1,2,3)] <- c('Lake_ID','watershed_area_km2','HII')
lulc.watershed  <- select(lulc.watershed , Lake_ID:HII,fraction_agriculture:fraction_urban)

# kestrel <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LPkestrel.xlsx') %>% arrange(Lake_ID)
# kestrel <- kestrel %>% mutate_at(vars(temperature,rel_humid,atm_press,wind_max,wind_avg), as.numeric) %>% group_by(Lake_ID) %>%
#   summarize(air.temp = mean.nona(temperature), rel_humid = mean.nona(rel_humid), atm_press = mean.nona(atm_press), 
#             wind.max = mean.nona(wind_max), wind.avg = mean.nona(wind_avg)) %>% ungroup
# #kestrel will not be very useful
# kestrel <- replace(kestrel, is.na(kestrel), NA)

rbr <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_RBR_QC.csv', stringsAsFactors = F) %>% 
  filter(depth == 'Tube length')

rbr <- rbr %>% select(lakepulse_id, conductivity_mean, temp_mean, do_mean, do_concentration_mean, chla_mean, ph_mean, salinity_mean, specific_conductance_mean)
colnames(rbr)[1] <- 'Lake_ID'
rbr <- rbr %>% mutate_at(vars(conductivity_mean:specific_conductance_mean), as.numeric)
rbr2 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_RBR_QC.csv', stringsAsFactors = F) %>% 
  filter(depth == 'Bottom meter') %>% select(lakepulse_id, do_concentration_mean)
colnames(rbr2)[c(1,2)] <- c('Lake_ID','bottom_do')
rbr2$bottom_do <- as.numeric(rbr2$bottom_do)
rbr <- left_join(rbr,rbr2)
rm(rbr2)
colnames(rbr)[c(2,3,4,5,6,7,8,9)] <- c('conductivity','rbr_temp','DO_sat','DO_mg','rbr_chla','pH','salinity','SPC')

strati <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_stratification.xlsx') %>% rename(Lake_ID = lakepulse_id) %>% arrange(Lake_ID)
colnames(strati)[2:5] <- c('epilimnion_depth','thermocline_depth','hypolimnion_depth','stratified')
strati <- strati %>% mutate_at(vars(epilimnion_depth:hypolimnion_depth), as.numeric)
strati <- replace(strati, is.na(strati), NA)
strati$stratified[strati$stratified == 'NA'] <- NA
strati$stratified <- as.factor(strati$stratified)

chla <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_chla_QC.csv', stringsAsFactors = F)
colnames(chla)[c(1,3)] <- c('Lake_ID','chla')
chla$chla[str_detect(chla$chla_flags,'-1')]
chla <- chla %>% mutate_at('chla', as.numeric) %>% group_by(Lake_ID) %>% summarize(chla = mean.nona(chla)) %>% ungroup %>% arrange(Lake_ID)
# test <- left_join(chla, rbr)
# plot(rbr_chla~chla,test,log='xy') #not that terrible

secchi <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_secchi_QC.csv', stringsAsFactors = F)
colnames(secchi)[c(1,6)] <- c('Lake_ID','secchi_depth')
secchi <- filter(secchi, site_name == 'Index Site') %>% select(Lake_ID, secchi_depth) %>% arrange(Lake_ID) %>% mutate_at('secchi_depth', as.numeric)

TN <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_TN_QC.csv', stringsAsFactors = F)
colnames(TN)[c(1,3)] <- c('Lake_ID','TN')
TN <- filter(TN, depth == 'Epilimnion') %>% select(Lake_ID, TN) %>% arrange(Lake_ID) %>% mutate_at('TN', as.numeric)

allenv <- left_join(basic.data, chla) %>% left_join(lulc.watershed) %>% left_join(rbr) %>%
  left_join(secchi) %>% left_join(strati) %>% left_join(TN)

## fish

load(file='~/Google Drive/Recherche/Lake Pulse Postdoc/R/FisHab_Git/formatted_open_data/species_codes.RData')
load(file='~/Google Drive/Recherche/Lake Pulse Postdoc/R/FisHab_Git/formatted_open_data/fishbase.RData')
fish <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/LP_fish_communitymatrix.xlsx')
fish.long <- fish %>% pivot_longer(FS005:FS364, names_to = 'species', values_to = 'presence')

## zooplankton biomass

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
zooLP <- bind_rows(d2017,d2018,d2019)
zooLP[is.na(zooLP)] <- 0
zooLP <- zooLP[,c(1,order(names(zooLP[,2:ncol(zooLP)]))+1)]
rm(d2017,d2018,d2019)
  
#### calculating max TL ####

fish.long$species.name <- species_codes$clean.species.name[match(fish.long$species,species_codes$fish_species_ID)]
#only keeping stuff id'd to species
fish.long <- drop_na(fish.long)
#merging species codes that are the same species by grouping by name instead of FS code
fish.long <- fish.long %>% group_by(id_lakepulse,species.name) %>% summarize(pres = max(presence)) %>% ungroup
#to wide
fish.wide <- fish.long %>% spread(species.name,pres)
fish.wide <- fish.wide %>% replace(is.na(.), 0)
#adding trophic level
for(co in 2:ncol(fish.wide)){
  sp <- colnames(fish.wide)[co]
  sp <- str_replace(sp, '_', ' ')
  tl <- estimate.df %>% filter(Species == sp) %>% pull(Troph)
  fish.wide[,co] <- fish.wide[,co]*tl
}

fish.wide$max.trophic.level <- apply(fish.wide[,2:ncol(fish.wide)], 1, max)

# send to Cindy
# cindy <- select(fish.wide, id_lakepulse, max.trophic.level) %>% filter(max.trophic.level > 0)
# writexl::write_xlsx(cindy, '~/Desktop/Cindy.xlsx')

#### machine learning ####

library(randomForest)
library(party)
library(caret)
library(gbm)
library(skimr)

TL <- fish.wide %>% select(id_lakepulse, max.trophic.level) %>% rename(Lake_ID = id_lakepulse, fishTL = max.trophic.level) %>% arrange(Lake_ID)

subenv <- allenv
subenv <- allenv %>% select(Lake_ID,latitude,area:ecozone,altitude:fraction_urban,stratified,TN)
#removing continuous vars with lots of NAs (pH, rbr chla, secchi, all the stratification vars)

dat <- left_join(TL, subenv) %>% select(-Lake_ID) %>% filter(fishTL > 0)
skim(dat)
dat <- drop_na(dat)
dat <- filter(dat, fishTL > 3)
  
sample <- createDataPartition(dat$fishTL,p=0.75,list=F)
training.set <- dat[sample,]
test.set <- dat[-sample,]

## simplest of random forest with randomForest

forest<-randomForest(fishTL ~ ., training.set, localImp = TRUE)
importance(forest, scale=F)
results <- predict(forest, newdata = test.set[,-1])
forest

## conditional random forest with party

cforest1 <-cforest(fishTL ~ ., training.set, controls = cforest_unbiased(mtry=6)) # mtry optimizied with a loop
train.results <- predict(cforest1, OOB=T)
test.results <- predict(cforest1, newdata = test.set)
postResample(pred = train.results, obs = training.set$fishTL)
postResample(pred = test.results, obs = test.set$fishTL)
plot(test.results~test.set$fishTL)
vars <- varimp(cforest1, conditional = T)
sort(vars)

## GBM

gbm1 <- gbm(fishTL ~ ., data=training.set, n.trees = 500, interaction.depth = 3)
gbm.perf(gbm1)
plot(predict(gbm1, n.trees = 500)~training.set$fishTL)
plot(density(training.set$fishTL))
plot(predict(gbm1, newdata = test.set, n.trees = 500)~test.set$fishTL)

#trying to optimize with caret

fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

gbmGrid <-  expand.grid(interaction.depth = 1:9, 
                       n.trees = 500, 
                       shrinkage = c(0.1,0.5,0.9),
                       n.minobsinnode = 20)


gbm2 <- train(fishTL ~ ., data=training.set, 
                 method = "gbm", 
                 trControl = fitControl, 
                 verbose = FALSE, 
                 tuneGrid = gbmGrid)

gbm2
plot(gbm2)

plot(predict(gbm2, newdata = test.set, n.trees = 500)~test.set$fishTL) #exact same thing

## Bayesian Additive Regression Trees (BART)

library(BayesTree)

bart1 <- bart(x.train = as.data.frame(training.set[,-1]), y.train = training.set$fishTL, x.test=as.data.frame(test.set[,-1]))
bart1
plot(bart1) #awful

## BART with caret & bartMachine

library(bartMachine)

