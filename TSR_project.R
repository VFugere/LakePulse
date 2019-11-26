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
zoo <- rename(zoo, group = division)
zoo$animal <- 'yes'
zoo$mean.size <- zoo$mean.size*1000
zoo <- rename(zoo, size.um = mean.size)
zoo <- zoo %>% select(Lake_ID:group,animal,size.um)
zoo$group <- tolower(zoo$group)

### phytoplankton

phyto <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/phytoplankton2017/phyto2017_clean.xlsx', sheet='Longform')
phytoT <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/phytoplankton2017/phyto2017_clean.xlsx', sheet='taxonomy_clean')
phytoT$group <- tolower(phytoT$group)

#use clean taxonomy sheet to replace erroneous species name in community matrix
phyto$totalbinomial <- phytoT$clean.name[match(phyto$totalbinomial,phytoT$totalbinomial)]

colnames(phyto) <- tolower(colnames(phyto))

phyto <- phyto %>% rename(Lake_ID = lake_id, size.um = gald.µm, name = totalbinomial) %>%
  select(Lake_ID, name, size.um) %>%
  drop_na

phyto <- left_join(phyto, phytoT, by = c('name' = 'clean.name'))
phyto <- filter(phyto, !is.na(species))
phyto <- phyto %>% group_by(Lake_ID, name) %>% summarize(size.um = mean(size.um)) %>% ungroup
phyto <- phytoT %>% select(group, clean.name) %>% right_join(phyto, by = c('clean.name' = 'name'))
phyto <- rename(phyto, species = clean.name)
phyto$animal <- 'no'
phyto$species <- str_replace(phyto$species, ' ', '_')
phyto <- select(phyto, Lake_ID, species, group, animal, size.um)

### merge both dataframes and format

data <- bind_rows(zoo,phyto) %>% arrange(Lake_ID, size.um)
data <- inner_join(data, lat, by = 'Lake_ID')
data <- left_join(data, rbr, by = 'Lake_ID')

## which variable will be used in analysis?
data$x <- data$latitude
xvarlab <- 'latitude (degrees)'
xvarlab2 <- 'latitudinal range (degrees)'
# data$x <- data$watertemp
# xvarlab <- 'water temperature (°C)'
# xvarlab2 <- 'thermal range (°C)'

data <- data %>% filter(!is.na(x))
data <- data %>% add_count(species)
data <- data %>% filter(n >= 5)

#add scaled variables within species
taxlist <- data %>% distinct(species) %>% pull(species)
data$scl.size <- 0
data$scl.lat <- 0
for(i in 1:length(taxlist)){
  tmp <- filter(data, species == taxlist[i])
  scl.size <- scale(tmp$size.um, scale = T)[,1]
  scl.lat <- scale(tmp$latitude, scale = F)[,1]
  data[data$species == taxlist[i], 'scl.size'] <- scl.size
  data[data$species == taxlist[i], 'scl.lat'] <- scl.lat
}

data$species <- as.factor(data$species)
data$scl.size[!is.finite(data$scl.size)] <- NA

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
                          'r2' = numeric(0), 
                          'mean.size' =numeric(0),
                          stringsAsFactors = F)

for(i in 1:nlevels(data$species)){
  focalsp <- levels(data$species)[i]
  spdat <- data %>% filter(species == focalsp)
  spdat <- spdat %>% drop_na(scl.size)
  if(nrow(spdat) >= 5){
  #lmmod <- lm(size.um~x,spdat)
  lmmod <- lm(scl.size~x,spdat)
  r2 <- summary(lmmod)$r.squared
  slope <- coef(lmmod)[2]
  p <- summary(lmmod)$coefficients[2,4]
  SE <- summary(lmmod)$coefficients[2,2]
  n <- nrow(spdat)
  min.lat <- min(spdat$x)
  max.lat <- max(spdat$x)
  group <- spdat$group[1]
  mean.size <- mean(spdat$size.um)
  results1 <- c(focalsp,group)
  results2 <- c(min.lat,max.lat,n,slope,SE,p,r2,mean.size)
  reg.results[nrow(reg.results)+1,1:2] <- results1
  reg.results[nrow(reg.results),3:10] <- results2
  }
}

#reg.results <- filter(reg.results, abs(slope) < 100)
reg.results$xrange <- reg.results$max.lat-reg.results$min.lat
reg.results <- reg.results %>% filter(xrange > 2)

reg.results$bubblesize <- rescale(reg.results$n, c(1,4))
reg.results$bubblecol <- 'gray'
reg.results$bubblecol[reg.results$p < 0.05 & reg.results$slope < 0] <- 'red'
reg.results$bubblecol[reg.results$p < 0.05 & reg.results$slope > 0] <- 'blue'

hist(reg.results$slope, breaks=50)

#plot à la Dornelas
plot(r2~slope,reg.results,bty='n',pch=16,col=alpha(bubblecol,0.5),cex=bubblesize)
abline(v=0,lty=2)
#not very exciting...

pdf('~/Desktop/results.pdf',width=16,height=5,pointsize=12)
layout(cbind(1,2,3),widths=c(0.4,0.3,0.3))
par(cex=1)
  
emptyPlot(xlim = range(data$x),yaxt='n',xaxt='n',ann=F, ylim=range(data$size.um),bty='l',log='y')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=c(1,10,100,1000))
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab=xvarlab)
title(ylab=expression(mean~body~size~(µm)),line=2.8)
for(i in 1:nlevels(data$species)){
  focalsp <- levels(data$species)[i]
  spdat <- data %>% filter(species == focalsp)
  spdat <- filter(spdat, !is.na(scl.size))
  if(nrow(spdat) >= 5){
    lmmod <- lm(scl.size~x,spdat)
    lncol <- 'gray'
    if(summary(lmmod)$coefficients[2,4] < 0.05){
      if(summary(lmmod)$coefficients[2,1] > 0){lncol <- 'blue'}else{
        lncol <- 'red'}
    }
    lmmod <- lm(size.um~x,spdat)
    points(fitted(lmmod)~spdat$x,col=alpha(lncol,0.5),type='l',lwd=1)
  }
}

plot(slope~xrange,reg.results,bty='n',pch=16,col=alpha(bubblecol,0.5),cex=bubblesize,bty='l',yaxt='n',xaxt='n',ann=F)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab=xvarlab2,line=2.8)
title(ylab=standardized~slope~(beta),line=2.8)
abline(h=0,lty=3)

plot(slope~mean.size,reg.results,bty='n',pch=16,col=alpha(bubblecol,0.5),cex=bubblesize,bty='l',log='x',yaxt='n',xaxt='n',ann=F)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(1,10,100,1000))
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab=standardized~slope~(beta),line=2.8)
title(xlab='mean body size (µm)',line=2.8)
abline(h=0,lty=3)

dev.off()




mod <- lmer(log(size.um) ~ x + (x-1|species) + (1|species), data = data)
plot(mod)
summary(mod)
confint(mod)
hist(resid(mod))
performance(mod)
icc(mod)

library(MCMCglmm)

#MCMCglmm parameters
a <- 1000
priorz <- list(R = list(V = diag(1), nu = 0.002), G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a), G2 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))
prior.s <- list(R=list(V=diag(4),n=0.002),G=list(G1=list(V=diag(1),n=0.002),G2=list(V=diag(1),n=0.002)))
m <- reg.results$SE^2 #variance of coefficients
mod1 <- MCMCglmm(slope ~ 1, data=reg.results, mev=m, prior=prior.s, random=~ system + traitID, rcov=~idh(conditions):units, nitt=103000, verbose=F)

MCMCglmm()

library(randomcoloR)


#cols2 <- randomColor(n_distinct(data$species))

emptyPlot(xlim = range(data$x),yaxt='n',xaxt='n',ann=F, ylim=range(data$size.um),bty='l',log='y')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=c(1,10,100,1000))
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab=xvarlab)
title(ylab=expression(mean~body~size~(µm)),line=2.8)
# for(focalsp in levels(mod$model$species)){
#   spdat <- mod$model %>% filter(species == focalsp)
#   spdat <- arrange(spdat, x)
#   points(size.um~x,spdat,col=alpha(1,0.5),type='l')
# }
for(i in 1:nlevels(data$species)){
  focalsp <- levels(data$species)[i]
  spdat <- data %>% filter(species == focalsp)
  spdat <- arrange(spdat, x)
  #points(size.um~x,spdat,col=alpha(cols2[i],0.5),type='p',pch=16,cex=0.5)
  if(nrow(spdat) >= 5){
    lmmod <- lm(size.um~x,spdat)
    lncol <- 'gray'
    if(summary(lmmod)$coefficients[2,4] < 0.05){
      if(summary(lmmod)$coefficients[2,1] > 0){lncol <- 'blue'}else{
        lncol <- 'red'}
    }
    points(fitted(lmmod)~spdat$x,col=alpha(lncol,0.5),type='l',lwd=1)
  }
}
# predict(mod,)
# plot_smooth(mod, view="x", lwd=3, col=cols[1], rm.ranef=T, se=1.96, rug=F, add=T)
# legend('topright',bty='n',legend=bquote(atop(italic('F') == .(testres[1]),italic('p') == .(testres[2]))))


emptyPlot(xlim = range(data$x),yaxt='n',xaxt='n',ann=F, ylim=range(data$scl.size, na.rm = T),bty='l')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab=xvarlab)
title(ylab=expression(body~size~(mm)),line=2.8)
# for(focalsp in levels(mod$model$species)){
#   spdat <- mod$model %>% filter(species == focalsp)
#   spdat <- arrange(spdat, x)
#   points(scl.size~x,spdat,col=alpha(1,0.5),type='l')
# }
for(i in 1:nlevels(data$species)){
  focalsp <- levels(data$species)[i]
  spdat <- data %>% filter(species == focalsp)
  spdat <- arrange(spdat, x)
  #points(scl.size~x,spdat,col=alpha(cols2[i],0.5),type='p',pch=16,cex=0.5)
  if(sum(is.na(spdat$scl.size))<1){
  lmmod <- lm(scl.size~x,spdat)
  points(fitted(lmmod)~spdat$x,col=alpha(cols2[i],0.5),type='l',lwd=0.8)
  }
}


mod <- bam(size.um ~ s(x, k = 8) + s(x, species, bs = 'fs', k = 6), data = data,nthreads=2)
mod <- bam(scl.size ~ s(x, k = 8) + s(x, species, bs = 'fs', k = 6), data = data,nthreads=2)
summary(mod)
modsum <- summary(mod, re.test=F)
testres <- modsum$s.table[c('s(x)'),c('F','p-value')]
testres[1] <- round(testres[1],2)
testres[2] <- round(testres[2],4)
if(testres[2] < 0.0001){testres[2] <- 0.0001}
rsq <- round(modsum$r.sq,2)

ylims <- predict(mod,se.fit=T)
ylims <- range(c(ylims$fit+1.96*ylims$se.fit,ylims$fit-1.96*ylims$se.fit))
ylims <- c(1,ylims[2])

emptyPlot(xlim = range(data$x),yaxt='n',xaxt='n',ann=F, ylim=ylims,bty='l',log='y')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab=xvarlab)
title(ylab=expression(body~size~(µm)),line=2.8)
# for(focalsp in levels(mod$model$species)){
#   spdat <- mod$model %>% filter(species == focalsp)
#   xs <- seq(min(spdat$x),max(spdat$x),by=0.01)
#   conditions <- list(x=xs, species=focalsp)
#   ys <- get_predictions(mod, cond = conditions, se=F, print.summary = F, rm.ranef=F) 
#   points(fit~x,ys,col=alpha(1,0.5),type='l')
# }
for(i in 1:nlevels(data$species)){
  focalsp <- levels(data$species)[i]
  spdat <- data %>% filter(species == focalsp)
  spdat <- arrange(spdat, x)
  lmmod <- gam(size.um~s(x, k = 4),data=spdat)
  points(fitted(lmmod)~spdat$x,col=alpha(cols2[i],0.5),type='l',lwd=0.8)
}
plot_smooth(mod, view="x", lwd=3, col=1, rm.ranef=T, se=1.96, rug=F, add=T)
legend('topright',bty='n',legend=bquote(atop(italic('F') == .(testres[1]),italic('p') == .(testres[2]))))


#

