rm(list=ls())

options(tibble.print_max = 100, scipen = 999)

#libraries
library(tidyverse)
library(RColorBrewer)
library(devtools)
library(readxl)
library(lme4)
library(mgcv)
library(performance)
library(itsadug)
library(scales)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
source_url("https://raw.githubusercontent.com/VFugere/LakePulse/master/custom_plots.R")
mean.nona <- function(x){mean(x, na.rm=T)}

cols <- brewer.pal(3, 'Dark2')

#### temp data ####

load('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/basic_data.RData')
kestrel <- read.csv2('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/2017 not merged/LakePulse2017_kestrel_QC.csv')
colnames(kestrel)[1] <- 'Lake_ID'
rbr <- read.csv2('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/2017 not merged/LakePulse2017_RBR_top_bottom_QC.csv')
colnames(rbr)[1] <- 'Lake_ID'
alex <- read.csv('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/2017 not merged/Alex_weather_station_data.csv')
colnames(alex)[1] <- 'Lake_ID'

kestrel <- kestrel %>% select(Lake_ID,temp) %>% rename(airtemp = temp)
rbr <- rbr %>% filter (depth == 'Tube lenght') %>% select(Lake_ID,temperature) %>% rename(watertemp = temperature)
alex <- alex %>% filter(Time >= 1950) %>% group_by(Lake_ID) %>% summarize(histtemp = mean(Mean_Ann_Temp, na.rm = T))
lat <- basic.data %>% select(Lake_ID, latitude)

merge <- inner_join(rbr,kestrel) %>% inner_join(alex) %>% inner_join(lat)
# corrgram::corrgram(merge)
# plot(merge[,2:5])

#air temp and water temp correlate well on sampling day. So I can use air temp as proxy for water temp.
#mean annual air temp since 1950 correlate almost perfectly with latitude. So can use latitude as a proxy for air temp
#latitude should provide a proxy for water temp for all lakes.

pdf('~/Desktop/tempgraph.pdf',width=10,height = 5,pointsize = 12)
par(mfrow=c(1,2), cex=1)

merge <- merge %>% drop_na

scatter.smooth(y=merge$watertemp,x=merge$airtemp,lpars = list(lwd = 5,col=alpha(1,0.5)),pch=16,col=alpha(cols[3],0.5),xlab='air temperature (°C)',ylab='water column  temperature (°C)',bty='l')
r <- round(cor(y=merge$watertemp,x=merge$airtemp),2)
legend('topright',bty='n',legend=bquote(R^2 == .(r)))

scatter.smooth(y=merge$histtemp,x=merge$latitude,lpars = list(lwd = 5,col=alpha(1,0.5)),pch=16,col=alpha(cols[3],0.5),xlab='latitude (degrees)',ylab='mean annual temperature (°C)',bty='l')
r <- round(cor(y=merge$histtemp,x=merge$latitude),2)
legend('topright',bty='n',legend=bquote(R^2 == .(r)))

dev.off()

#### plankton data ####

bad.zoo.samples <- c('07-057','17-050')#,'08-205','07-029') # remove these two if anything weird (see email Cindy 30-Sept-2019)

zoo <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton2017/all raw data 2017.xlsx', sheet='raw') %>%
  filter(!(ID_lakepulse %in% bad.zoo.samples)) %>%
  rename(Lake_ID = ID_lakepulse)

#calculate mean body size, excluding 0s and NAs
sizes <- select(zoo, D1:D10)
sizes[sizes == 0] <- NA
sizes$mean <- apply(sizes, MARGIN = 1, FUN = mean.nona)
zoo$mean.size <- sizes$mean
rm(sizes)

#removing all taxa that aren't id'd to species
zoo$species <- paste(zoo$genus,zoo$species,sep='_')
distinct(zoo, species) %>% pull(species) -> taxlist
to.rm <- taxlist[str_detect(taxlist, '.spp')]
to.rm <- c(to.rm, taxlist[str_detect(taxlist, 'copepodid')])
to.rm <- c(to.rm, taxlist[str_detect(taxlist, '_NA')])
to.rm <- c(to.rm, c('Alona_sp.','Skistodiaptomus_sp.','immature_cladoceran','Macrothrix_sp.','Acanthocyclops_sp.'))

#clean up, remove species with a single site, and add latitude
zoo <- filter(zoo, species %!in% to.rm)
zoo$division[zoo$species == 'Leptodora_kindtii'] <- 'Cladocera'
zoo <- zoo %>% select(Lake_ID,species,division,mean.size) %>% drop_na %>% filter(is.numeric(mean.size))
zoo <- inner_join(zoo, lat, by = 'Lake_ID')
zoo <- zoo %>% add_count(species)
zoo <- zoo %>% filter(n >= 5)

#add scaled variables within species
taxlist <- zoo %>% distinct(species) %>% pull(species)
zoo$scl.size <- 0
zoo$scl.lat <- 0
for(i in 1:length(taxlist)){
  tmp <- filter(zoo, species == taxlist[i])
  scl.size <- scale(tmp$mean.size, scale = F)[,1]
  scl.lat <- scale(tmp$latitude, scale = F)[,1]
  zoo[zoo$species == taxlist[i], 'scl.size'] <- scl.size
  zoo[zoo$species == taxlist[i], 'scl.lat'] <- scl.lat
}
zoo$species <- as.factor(zoo$species)

### phytoplankton

phyto <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/phytoplankton2017/Phytoplankton LakePulse 2017_Jelena.xlsx', sheet='Longform')

phyto <- phyto %>% rename(Lake_ID = lake_id, bv = Biomass.mgm3, name = totalbinomial) %>%
  select(Lake_ID, name, bv) %>%
  filter(bv != 0) %>% 
  group_by(Lake_ID, name) %>% 
  summarize(bv = sum(bv, na.rm=T)) %>%
  ungroup %>%
  spread(name, bv) %>%
  as.data.frame

phyto[is.na(phyto)] <- 0

#### analysis ####

## approach 1: fit a linear regression to each species with 5+ data points, and check distribution of slopes

reg.results <- data.frame('focalsp' = character(0),
                          'group' = character(0),
                          'min.lat' = numeric(0),
                          'max.lat' = numeric(0),
                          'n' = numeric(0),
                          'slope' = numeric(0),
                          'SE' = numeric(0),
                          'p' = numeric(0),
                          'r2' = numeric(0), stringsAsFactors = F)

for(i in 1:nlevels(zoo$species)){
  focalsp <- levels(zoo$species)[i]
  spdat <- zoo %>% filter(species == focalsp)
  spdat <- arrange(spdat, latitude)
  lmmod <- lm(mean.size~latitude,spdat)
  r2 <- summary(lmmod)$r.squared
  slope <- coef(lmmod)[2]
  p <- summary(lmmod)$coefficients[2,4]
  SE <- summary(lmmod)$coefficients[2,2]
  n <- nrow(spdat)
  min.lat <- min(spdat$latitude)
  max.lat <- max(spdat$latitude)
  group <- spdat$division[1]
  results1 <- c(focalsp,group)
  results2 <- c(min.lat,max.lat,n,slope,SE,p,r2)
  reg.results[nrow(reg.results)+1,1:2] <- results1
  reg.results[nrow(reg.results),3:9] <- results2
}

hist(reg.results$slope, breaks=20)
#plot à la Dornelas
reg.results$bubblesize <- rescale(reg.results$n, c(1,4))
reg.results$bubblecol <- 'gray'
reg.results$bubblecol[reg.results$p < 0.05 & reg.results$slope < 0] <- 'red'
reg.results$bubblecol[reg.results$p < 0.05 & reg.results$slope > 0] <- 'blue'
plot(r2~slope,reg.results,bty='n',pch=16,col=alpha(bubblecol,0.5),cex=bubblesize)
abline(v=0,lty=2)
#not very exciting...

mod <- lmer(mean.size ~ latitude + (latitude-1|species) + (1|species), data = zoo)
# plot(mod)
# summary(mod)
# hist(resid(mod))
# performance(mod)
# icc(mod)


library(randomcoloR)


cols2 <- randomColor(n_distinct(zoo$species))
emptyPlot(xlim = range(zoo$latitude),yaxt='n',xaxt='n',ann=F, ylim=range(zoo$mean.size),bty='l')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab='latitude (degrees)')
title(ylab=expression(body~size~(mm)),line=2.8)
# for(focalsp in levels(mod$model$species)){
#   spdat <- mod$model %>% filter(species == focalsp)
#   spdat <- arrange(spdat, latitude)
#   points(mean.size~latitude,spdat,col=alpha(1,0.5),type='l')
# }
for(i in 1:nlevels(zoo$species)){
  focalsp <- levels(zoo$species)[i]
  spdat <- zoo %>% filter(species == focalsp)
  spdat <- arrange(spdat, latitude)
  #points(mean.size~latitude,spdat,col=alpha(cols2[i],0.5),type='p',pch=16,cex=0.5)
  lmmod <- lm(mean.size~latitude,spdat)
  points(fitted(lmmod)~spdat$latitude,col=alpha(cols2[i],0.5),type='l',lwd=0.8)
}
# predict(mod,)
# plot_smooth(mod, view="latitude", lwd=3, col=cols[1], rm.ranef=T, se=1.96, rug=F, add=T)
# legend('topright',bty='n',legend=bquote(atop(italic('F') == .(testres[1]),italic('p') == .(testres[2]))))


emptyPlot(xlim = range(zoo$latitude),yaxt='n',xaxt='n',ann=F, ylim=range(zoo$scl.size),bty='l')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab='latitude (degrees)')
title(ylab=expression(body~size~(mm)),line=2.8)
# for(focalsp in levels(mod$model$species)){
#   spdat <- mod$model %>% filter(species == focalsp)
#   spdat <- arrange(spdat, latitude)
#   points(scl.size~latitude,spdat,col=alpha(1,0.5),type='l')
# }
for(i in 1:nlevels(zoo$species)){
  focalsp <- levels(zoo$species)[i]
  spdat <- zoo %>% filter(species == focalsp)
  spdat <- arrange(spdat, latitude)
  #points(scl.size~latitude,spdat,col=alpha(cols2[i],0.5),type='p',pch=16,cex=0.5)
  lmmod <- lm(scl.size~latitude,spdat)
  points(fitted(lmmod)~spdat$latitude,col=alpha(cols2[i],0.5),type='l',lwd=0.8)
}


mod <- bam(mean.size ~ s(latitude, k = 8) + s(latitude, species, bs = 'fs', k = 6), data = zoo,nthreads=2)
mod <- bam(scl.size ~ s(latitude, k = 8) + s(latitude, species, bs = 'fs', k = 6), data = zoo,nthreads=2)
summary(mod)
modsum <- summary(mod, re.test=F)
testres <- modsum$s.table[c('s(latitude)'),c('F','p-value')]
testres[1] <- round(testres[1],2)
testres[2] <- round(testres[2],4)
if(testres[2] < 0.0001){testres[2] <- 0.0001}
rsq <- round(modsum$r.sq,2)

ylims <- predict(mod,se.fit=T)
ylims <- range(c(ylims$fit+1.96*ylims$se.fit,ylims$fit-1.96*ylims$se.fit))

emptyPlot(xlim = range(zoo$latitude),yaxt='n',xaxt='n',ann=F, ylim=ylims,bty='l')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab='latitude (degrees)')
title(ylab=expression(body~size~(mm)),line=2.8)
# for(focalsp in levels(mod$model$species)){
#   spdat <- mod$model %>% filter(species == focalsp)
#   xs <- seq(min(spdat$latitude),max(spdat$latitude),by=0.01)
#   conditions <- list(latitude=xs, species=focalsp)
#   ys <- get_predictions(mod, cond = conditions, se=F, print.summary = F, rm.ranef=F) 
#   points(fit~latitude,ys,col=alpha(1,0.5),type='l')
# }
for(i in 1:nlevels(zoo$species)){
  focalsp <- levels(zoo$species)[i]
  spdat <- zoo %>% filter(species == focalsp)
  spdat <- arrange(spdat, latitude)
  lmmod <- gam(mean.size~s(latitude, k = 4),data=spdat)
  points(fitted(lmmod)~spdat$latitude,col=alpha(cols2[i],0.5),type='l',lwd=0.8)
}
plot_smooth(mod, view="latitude", lwd=3, col=1, rm.ranef=T, se=1.96, rug=F, add=T)
legend('topright',bty='n',legend=bquote(atop(italic('F') == .(testres[1]),italic('p') == .(testres[2]))))


#

