### LOAD LAKE PULSE BASIC INFO
### ADD HYDROLAKES ATTRIBUTES
### ADD FRESHWATER ECOREGION

# code by Vincent Fugère (2019)

rm(list=ls())
library(tidyverse)
library(skimr)
library(readxl)
library(sp)

#basic data

d2017 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2017_basic_info.csv', sep = ';', stringsAsFactors = F)
d2017$year <- 2017
d2018 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2018_basic_info.csv', sep = ';', stringsAsFactors = F)
d2018$year <- 2018
d2019 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2019_basic_info.csv', sep = ';', stringsAsFactors = F)
d2019$year <- 2019

basic.data <- bind_rows(d2017,d2018) %>%
  bind_rows(d2019) %>%
  rename(Lake_ID = lakepulse_id, area = size_km2, HI = hi_index)

rm(d2017,d2018,d2019)

# ### save Qc lakes as shapefile
# 
# LP_Qc <- filter(basic.data, province  == 'QUEBEC')
# LP_shp <-  LP_Qc
# sp::coordinates(LP_shp) <- ~longitude+latitude
# sp::proj4string(LP_shp) <- CRS("+proj=longlat +datum=WGS84")
# raster::shapefile(LP_shp, "~/Desktop/LakePulse_QC.shp")
# write_csv(LP_Qc, '~/Desktop/LP_Qc.csv')

### get info from hydrolakes

LP_spatial <-  basic.data
coordinates(LP_spatial) <- ~longitude+latitude

hydro <- rgdal::readOGR(dsn = "/Users/vincentfugere/Desktop/GIS/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10_shp/", layer = "HydroLAKES_polys_v10")
proj4string(hydro)

proj4string(LP_spatial) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

overlay <- over(LP_spatial, hydro)
overlay <- droplevels(overlay)

plot(Pour_lat~Pour_long, overlay, pch=1)
points(latitude~longitude, basic.data, pch=16, col='red',cex=0.3)

hydroLakes_LPsites <- overlay
save(hydroLakes_LPsites, file='~/Desktop/HydroLakes_LP.RData')

basic.data <- bind_cols(basic.data,overlay)
basic.data <- rename(basic.data, HL_name = Lake_name, HL_area = Lake_area, HL_watershed_area = Wshd_area)
basic.data <- dplyr::select(basic.data, -(HL_name:Poly_src), -Grand_id, -Lake_type, -HL_watershed_area)

### adding FW ecoregion

feow <- rgdal::readOGR(dsn = "/Users/vincentfugere/Desktop/FEOW-TNC/", layer = "FEOWv1_TNC")
proj4string(feow)

LP_spatial <- spTransform(LP_spatial, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj4string(LP_spatial)

overlay <- over(LP_spatial, feow)
overlay <- droplevels(overlay)
overlay <- dplyr::select(overlay, ECOREGION) %>% rename(feow = ECOREGION)

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
plot(latitude~longitude, basic.data, pch=16, col=sample(color, 17)[overlay$feow])

basic.data <- bind_cols(basic.data,overlay)

save(basic.data, file='/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/basic_data.RData')

