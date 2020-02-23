rm(list=ls())

#libraries
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(adespatial)
library(devtools)
library(sp)
library(party)
library(rworldmap)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
source_url("https://raw.githubusercontent.com/VFugere/LakePulse/master/custom_plots.R")

#cols
cols <- brewer.pal(7, 'Dark2')[1:3] 
cols2 <- brewer.pal(8, 'Dark2')[c(4,5,6,8)] 

#map
map <- getMap(resolution = "low")

#data
load('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/2017_community_data.RData')
load('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/basic_data.RData')

# #land use data
# lulc <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2017_LULC_QC.csv', sep = ';', stringsAsFactors = F)
# lulc <- select(lulc, -comment,-flag)
# basic.data <- left_join(basic.data,lulc, by = c('Lake_ID' = 'lakepulse_id'))
# rm(lulc)

# #fish richness
# fish$richness <- apply(fish[,2:62], 1, sum)

#remove PNN sites from QC
bad.fish <- c('06-085','06-096','06-131','17-093','17-017','17-012','17-105','06-125','06-095','06-141')
fish <- filter(fish, Lake_ID %!in% bad.fish)

# fish <- select(fish, Lake_ID, richness, everything())
# 
# subd <- select(fish, Lake_ID, richness)
# submeta <- basic.data %>% select(-team, -sampling_date, -lake_name, -hi_class, -size_class, -province,
#                                  -comment,-flags,-type,-utc,-year,-Hylak_id,-HL_area,-Vol_res,-Vol_src)
# submeta$ecozone <- as.factor(submeta$ecozone)
# subd <- left_join(subd,submeta, by = 'Lake_ID') %>% select(-Lake_ID)
# subd$area_ratio <- with(subd, area/watershed_km2)
# subd$LUI <- with(subd,area_ratio*HI)
# plot(ctree(richness ~ ., subd))
# mapplot(subd$longitude,subd$latitude,subd$richness,'fish richness')
# 
# pdf('~/Desktop/plots.pdf', width=20,height = 20,pointsize = 8)
# plot(subd)
# dev.off()
# 
# subd$d_log <- scale(log(subd$depth_m))
# subd$hi_log <- scale(log(subd$HI+0.001))
# #subd$hi_log <- scale(subd$LUI)
# 
# colfunc <- colorRampPalette(RColorBrewer::brewer.pal(11,'RdYlBu'))
# mapcols <- colfunc(100)[100:1]
# 
# library(mgcv)
# library(itsadug)
# library(performance)
# mod1 <- gam(richness ~ s(latitude,longitude, bs='gp',k=6) + s(d_log, k = 6) + s(hi_log, k = 6) + ti(d_log,hi_log,k=6), data=subd)
# plot(mod1)
# gam.check(mod1)
# anova(mod1) #depth has significant impact but not HI
# model_performance(mod1)
# fvisgam(mod1, view = c('d_log','hi_log'), dec=1,xlab='max depth (log)',ylab='human impact (log)',lwd=1.5,color = mapcols, main = NULL)

basic.data$hi_class2 <- 'low'
basic.data$hi_class2[basic.data$HI > 0.05] <- 'moderate'
basic.data$hi_class2[basic.data$HI > 0.2] <- 'high'

trt <- basic.data %>% select(Lake_ID, hi_class2, size_class, ecozone, HI, area, latitude, longitude) %>%
  rename(hi_cont = HI, area_cont = area, lat=latitude,long=longitude) %>%
  rename(area = size_class, HI = hi_class2)

#remove lakes not in qc

lakes.to.keep <- basic.data %>% filter(province %in% c('QUEBEC','ONTARIO')) %>% pull(Lake_ID)

basic.data <- filter(basic.data, Lake_ID %in% lakes.to.keep)
phyto <- filter(phyto, Lake_ID %in% lakes.to.keep)
bacterio <- filter(bacterio, Lake_ID %in% lakes.to.keep)
zoo <- zoo.biomass.grouped %>% filter(Lake_ID %in% lakes.to.keep)

fish <- as.data.frame(fish)

#### figure ####

pdf('~/Desktop/betadiv1.pdf',width=13,height = 7,pointsize = 12)
par(cex=1,mar=c(4,4,1,1),mfrow=c(4,3))

mybxp <- function(x,y,colvec,ylab,ylims,xlab){
  boxplot(y~x, lty=1, staplewex=0, whisklwd=1,boxwex=0.6, xlab=xlab, outline=F,boxlwd=1, medlwd=1, col=colvec, ylab=ylab, ylim=ylims, cex.lab=1.2)
  points(x=jitter(as.numeric(x)),y=y, pch=16,cex=0.5, col=alpha(1,0.3))
}

for(i in 1:4){
  
  if(i == 1){com <- fish}
  if(i == 2){com <- zoo.biomass.grouped}
  if(i == 3){com <- phyto}
  if(i == 4){com <- bacterio}
  
  meta <- com %>% select(Lake_ID)
  meta <- left_join(meta, trt, by = 'Lake_ID')
  meta <- mutate_at(meta, vars(HI:ecozone), as.factor)
  meta$HI <- factor(meta$HI, levels=c('low','moderate','high'))
  com <- com[,2:ncol(com)] %>% as.matrix
  com[is.na(com)] <- 0
  row.names(com) <- 1:nrow(com)
  row.names(meta) <- 1:nrow(meta)
  com.sub <- com
  meta.sub <- meta
  com.sub <- decostand(com.sub,method='pa')
  dm <- vegdist(com.sub,method = 'jaccard')
  
  richness <- apply(com.sub, 1, sum)
  mybxp(y=richness,x=meta.sub$HI,colvec=cols,ylab='richness',ylims=range(richness),xlab='land use intensity')
  
  bd <- betadisper(d=dm, group=meta.sub$HI, type='centroid')
  if(i==1){
    plot(betadisper(d=dm, group=meta.sub$HI,type='centroid'),ellipse=F,hull=F,label=F,seg.col=alpha(cols,0.1),segments=T,col=cols,pch=rep(1,3),cex=0.3,main=NULL,sub=NULL,ylab='PCoA2',xlab='PCoA1',cex.lab=1.2)
    legend('topleft',bty='n',legend=levels(meta.sub$HI),pch=16,col=cols)
    
  }else{
    plot(betadisper(d=dm, group=meta.sub$HI,type='centroid'),ellipse=F,hull=F,label=F,seg.col=alpha(cols,0.1),segments=T,col=cols,pch=rep(1,3),cex=0.1,main=NULL,sub=NULL,ylab='PCoA2',xlab='PCoA1',cex.lab=1.2)
  }
  
  distances <- (betadisper(d=dm, group=meta.sub$HI, type='centroid'))$distances
  mybxp(y=distances,x=meta.sub$HI,colvec=cols,ylab='distance to centroid',ylims=range(distances),xlab='land use intensity')
  
}

dev.off()

