# Script to put together all environmental parameters for Lake Pulse sites

rm(list=ls())

#libraries
library(tidyverse)
library(readxl)
library(devtools)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
source_url("https://raw.githubusercontent.com/VFugere/LakePulse/master/custom_plots.R")
mean.nona <- function(x){mean(x, na.rm=T)}
sum.nona <- function(x){sum(x, na.rm=T)}

### LOAD AND PASTE TOGETHER ALL LAKE PULSE ENVIRONMENTAL VARIABLES ###

# basic data

d2017 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2017_basic_info.csv', sep = ';', stringsAsFactors = F)
d2017$sampling_date <- as.Date(d2017$sampling_date, '%Y-%m-%d')

d2018 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2018_basic_info.csv', sep = ';', stringsAsFactors = F)
d2018$sampling_date <- as.Date(d2018$sampling_date, '%Y-%m-%d')

d2019 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2019_basic_info.csv', sep = ';', stringsAsFactors = F)
d2019$sampling_date <- str_replace(d2019$sampling_date, 'juil', '07')
d2019$sampling_date <- str_replace(d2019$sampling_date, 'ao.t', '08')
d2019$sampling_date <- str_replace(d2019$sampling_date, 'juin', '06')
d2019$sampling_date <- as.Date(d2019$sampling_date, '%d-%m-%y')

basic.data <- bind_rows(d2017,d2018) %>%
  bind_rows(d2019) %>%
  rename(Lake_ID = lakepulse_id, area_km2 = size_km2, HII = hi_index)
rm(d2017,d2018,d2019)

basic.data$year <- format(basic.data$sampling_date, '%Y') %>% as.numeric
basic.data$julian.day <- format(basic.data$sampling_date, '%j') %>% as.numeric

basic.data <- basic.data %>% arrange(Lake_ID) %>% 
  select(-HII,-comment,-flags,-utc) %>% 
  select(Lake_ID, lake_name, province, year, julian.day, sampling_date, everything())

altitude <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_altitude.xlsx') %>% rename(Lake_ID = lakepulse_id)

#land use data

#entire watershed
lulc.watershed <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LULC_QC20200323.xlsx') %>% arrange(idLakePulse)
#NOTE THAT WATERSHED AREA AND HII ARE NOT EXACTLY EQUAL IN BASIC AND LULC DATASETS. LULC data more up-to-date so I use this instead.
colnames(lulc.watershed )[c(1,2,3)] <- c('Lake_ID','watershed_area_km2','HII')
lulc.watershed  <- select(lulc.watershed , Lake_ID:HII,fraction_agriculture:fraction_urban)

#land use proportions within 1.5 km of the lakeshore
lulc.buffer1 <- read_xls('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LU2017ringsBuffer_finalQC20200323.xls')
lulc.buffer2 <- read_xls('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LU2018ringsBuffer_finalQC20200323.xls')
lulc.buffer3 <- read_xls('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LU2019ringsBuffer_finalQC20200323.xls')
lulc.buffer <- bind_rows(lulc.buffer1,lulc.buffer2,lulc.buffer3)
rm(lulc.buffer1,lulc.buffer2,lulc.buffer3)
lulc.buffer <- select(lulc.buffer, -...22, -Flags, -Comments)
colnames(lulc.buffer)[3] <- 'watershed_area'
lulc.buffer <- lulc.buffer %>% mutate_at(vars(watershed_area:fraction_water),as.numeric)
lulc.buffer <- lulc.buffer %>% filter(Ring %!in% c('Ring 1500-3100 m','Ring 3100-6300 m','Ring 6300 m-end')) %>%
  select(idLakePulse, area_km2_NA:total_km2_area) %>% group_by(idLakePulse) %>% summarize_if(is.numeric, sum.nona)
lulc.buffer$prop.ag.buffer <- lulc.buffer$area_km2_agriculture/lulc.buffer$total_km2_area
lulc.buffer$prop.forestry.buffer <- lulc.buffer$area_km2_forestry/lulc.buffer$total_km2_area
lulc.buffer$prop.mines.buffer <- lulc.buffer$area_km2_mines/lulc.buffer$total_km2_area
lulc.buffer$prop.natural.buffer <- lulc.buffer$area_km2_naturalLandscapes/lulc.buffer$total_km2_area
lulc.buffer$prop.pasture.buffer <- lulc.buffer$area_km2_pasture/lulc.buffer$total_km2_area
lulc.buffer$prop.urban.buffer <- lulc.buffer$area_km2_urban/lulc.buffer$total_km2_area
lulc.buffer <- lulc.buffer %>% rename(Lake_ID = idLakePulse) %>% select(Lake_ID, prop.ag.buffer:prop.urban.buffer)

#kestrel

kestrel <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LPkestrel.xlsx') %>% arrange(Lake_ID)
kestrel <- kestrel %>% mutate_at(vars(temperature,rel_humid,atm_press,wind_max,wind_avg), as.numeric) %>% group_by(Lake_ID) %>%
  summarize(kestrel.air.temp = mean.nona(temperature), rel.humid = mean.nona(rel_humid), atm.press = mean.nona(atm_press),
            wind.max = mean.nona(wind_max), wind.avg = mean.nona(wind_avg)) %>% ungroup
kestrel <- replace(kestrel, is.na(kestrel), NA)

#RBR

#top-bottom
rbr <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_RBR_QC.csv', stringsAsFactors = F) %>%
  filter(depth == 'Tube length')
rbr <- rbr %>% select(lakepulse_id, conductivity_mean, temp_mean, do_mean, do_concentration_mean, chla_mean, ph_mean, salinity_mean, specific_conductance_mean)
colnames(rbr)[1:9] <- c('Lake_ID','conductivity_top','temp_top','DOsat_top','DOmg_top','Chla_RBR_top','pH_top','salinity_top','SPC_top')
rbr <- rbr %>% mutate_at(vars(conductivity_top:SPC_top), as.numeric)
rbr2 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_RBR_QC.csv', stringsAsFactors = F) %>%
  filter(depth == 'Bottom meter') %>% select(lakepulse_id, conductivity_mean, temp_mean, do_mean, do_concentration_mean, chla_mean, ph_mean, salinity_mean, specific_conductance_mean)
colnames(rbr2)[1:9] <- c('Lake_ID','conductivity_bottom','temp_bottom','DOsat_bottom','DOmg_bottom','Chla_RBR_bottom','pH_bottom','salinity_bottom','SPC_bottom')
rbr2 <- rbr2 %>% mutate_at(vars(conductivity_bottom:SPC_bottom), as.numeric)
rbr.top.bottom <- left_join(rbr,rbr2)
rm(rbr, rbr2)

#average of entire depth profile
rbr <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_RBR_whole_water_column.xlsx')
rbr <- select(rbr, -depth_averaged_m, -pres24, -Flags, -Comments, -pressure_dbar, -sound.velocity_m.per.s)
colnames(rbr)[1:9] <- c('Lake_ID','conductivity_Zavg','temp_Zavg','DOsat_Zavg','DOmg_Zavg','Chla_RBR_Zavg','pH_Zavg','salinity_Zavg','SPC_Zavg')

#presence of stratification
strati <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_stratification.xlsx') %>% rename(Lake_ID = lakepulse_id) %>% arrange(Lake_ID)
colnames(strati)[2:5] <- c('epilimnion_depth','thermocline_depth','hypolimnion_depth','stratified')
strati <- strati %>% mutate_at(vars(epilimnion_depth:hypolimnion_depth), as.numeric)
strati <- replace(strati, is.na(strati), NA)
strati$stratified[strati$stratified == 'NA'] <- NA
strati$stratified <- as.factor(strati$stratified)

# Chlorophyll a (extraction + spectro)

chla <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_chla_QC.csv', stringsAsFactors = F)
colnames(chla)[c(1,3)] <- c('Lake_ID','chla')
chla <- chla %>% mutate_at('chla', as.numeric) %>% group_by(Lake_ID) %>% summarize(chla = mean.nona(chla)) %>% ungroup %>% arrange(Lake_ID)
# test <- left_join(chla, rbr)
# plot(rbr_chla~chla,test,log='xy') #not that terrible
colnames(chla)[2] <- 'Chla_spectro'

# Secchi

secchi <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_secchi_QC.csv', stringsAsFactors = F)
colnames(secchi)[c(1,6)] <- c('Lake_ID','secchi_depth')
secchi <- filter(secchi, site_name == 'Index Site') %>% select(Lake_ID, secchi_depth) %>% arrange(Lake_ID) %>% mutate_at('secchi_depth', as.numeric)

# Nutrients

TN <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LP_TN_QC.csv', stringsAsFactors = F)
colnames(TN)[c(1,3)] <- c('Lake_ID','TN')
TN <- filter(TN, depth == 'Epilimnion') %>% select(Lake_ID, TN) %>% arrange(Lake_ID) %>% mutate_at('TN', as.numeric)

#### ANCIALLRY DATA FROM OTHER SOURCES

#climate data from weather canada
climate <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/other sources/LP_climate_Cindy.csv', stringsAsFactors = F) %>% 
  mutate_at(vars(distance:sampling_year), as.numeric) %>% select(Lake_ID,mean_temp:total_precip)
colnames(climate)[2:5] <- c('mean.air.temp','min.air.temp','max.air.temp','total.precip')

#hydrolakes, freshwater ecoregion, continental watershed
someGISdata <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/other sources/LP_GISaddons.xlsx')

#### PASTING EVERYTHING TOGETHER

LP.env.data <- left_join(basic.data, altitude) %>% 
  left_join(someGISdata) %>% left_join(lulc.watershed) %>% left_join(lulc.buffer) %>%
  left_join(kestrel) %>% left_join(climate) %>% left_join(strati) %>%
  left_join(secchi) %>%  left_join(TN) %>%  left_join(rbr.top.bottom) %>%
  left_join(rbr) %>%  left_join(chla)

save(LP.env.data, file='~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/allenvdata.RData')

