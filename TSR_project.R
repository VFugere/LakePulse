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
kestrel <- read.csv2('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2017_kestrel_QC.csv')
colnames(kestrel)[1] <- 'Lake_ID'
rbr <- read.csv2('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2017_RBR_top_bottom_QC.csv')
colnames(rbr)[1] <- 'Lake_ID'
alex <- read.csv('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/other sources/Alex_weather_station_data.csv')
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

#bad.zoo.samples <- c('07-057','17-050')#,'08-205','07-029') # remove these two if anything weird (see email Cindy 30-Sept-2019)

d2017 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/all raw data 2017.xlsx', sheet='raw') %>% select(ID_lakepulse,genus,species,division,`#individuals counted`,`# / L`,D1:D10)
d2018 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/all raw data 2018.xlsx', sheet='raw') %>% select(ID_lakepulse,genus,species,division,`#individuals counted`,`# / L`,D1:D10)
d2019 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/all raw data 2019.xlsx', sheet='raw') %>% select(ID_lakepulse,genus,species,division,`#individuals counted`,`# / L`,D1:D10)

zoo <- bind_rows(d2017,d2018,d2019) %>% rename(Lake_ID = ID_lakepulse, abundance = `#individuals counted`, density = `# / L`)
zoo <- filter(zoo, division %in% c('Cladocera','Copepoda'))
rm(d2017,d2018,d2019)

#calculate mean body size, excluding 0s and NAs
sizes <- select(zoo, D1:D10)
sizes[sizes == 0] <- NA
sizes$mean <- apply(sizes, MARGIN = 1, FUN = mean.nona)
zoo$mean.size <- sizes$mean*1000 #in um
rm(sizes)

## community weighted mean body size

zoo$tot.size <- zoo$density*zoo$mean.size
cmw <- zoo %>% group_by(Lake_ID) %>% summarize(total.zoo.length = sum(tot.size), total.density = sum(density))
cmw$cwms <- cmw$total.zoo.length/cmw$total.density
cmw <- left_join(cmw,basic.data)
plot(cwms~latitude,cmw,pch=16,cex=0.5)
library(party)
dat <- cmw %>% select(cwms, latitude, longitude, area:ecozone, Shore_len:feow)
dat <- dat %>% mutate_if(is.character, ~as.factor(.))
plot(ctree(cwms ~ ., dat))

#removing all taxa that aren't id'd to species
zoo$species <- paste(zoo$genus,zoo$species,sep='_')
distinct(zoo, species) %>% pull(species) -> taxlist
to.rm <- taxlist[str_detect(taxlist, '.spp')]
to.rm <- c(to.rm, taxlist[str_detect(taxlist, 'copepodid')])
to.rm <- c(to.rm, taxlist[str_detect(taxlist, '_NA')])
to.rm <- c(to.rm, c('Pleuroxus_sp.','Latona_sp.','Alona_sp.','Skistodiaptomus_sp.','immature_cladoceran','Macrothrix_sp.','Acanthocyclops_sp.'))
zoo <- filter(zoo, species %!in% to.rm)

#alternative: utiliser mean biomass instead of mean length. This is already calculated in column 'biomass factor'
zoo$mean.size <- zoo$`biomass factor (µg/ind)`

#clean up, remove species with a single site, and add latitude
zoo <- zoo %>% select(Lake_ID,species,division,mean.size) %>% drop_na %>% filter(is.numeric(mean.size))
zoo <- rename(zoo, group = division)
zoo$animal <- 'yes'
zoo <- rename(zoo, ind.size = mean.size)
zoo <- zoo %>% select(Lake_ID:group,animal,ind.size)
zoo$group <- tolower(zoo$group)

### phytoplankton

phyto <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/phytoplankton/phyto2017_clean.xlsx', sheet='Longform')
phytoT <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/phytoplankton/phyto2017_clean.xlsx', sheet='taxonomy_clean')
phytoT$group <- tolower(phytoT$group)

#use clean taxonomy sheet to replace erroneous species name in community matrix
phyto$totalbinomial <- phytoT$clean.name[match(phyto$totalbinomial,phytoT$totalbinomial)]

colnames(phyto) <- tolower(colnames(phyto))

# phyto <- phyto %>% rename(Lake_ID = lake_id, ind.size = gald.µm, name = totalbinomial) %>%
#   select(Lake_ID, name, ind.size) %>%
#   drop_na

# #alternative: use individual cell mass instead of gald length
phyto$mg.per.ml <- phyto$biomass.mgm3/1000000
phyto$bv.um3.per.ml <- phyto$mg.per.ml * 1000000000 #same bv measure as as NLA
phyto$cells.per.ml <- phyto$density.orgl / 1000 # same density measure as NLA
phyto$um3.per.cell <- phyto$bv.um3.per.ml/phyto$cells.per.ml
phyto$g.per.cell <- phyto$um3.per.cell * 0.000000000001
phyto$ug.per.cell <- phyto$g.per.cell * 1000000

phyto <- phyto %>% rename(Lake_ID = lake_id, ind.size = ug.per.cell, name = totalbinomial) %>%
  select(Lake_ID, name, ind.size) %>%
  drop_na

phyto <- left_join(phyto, phytoT, by = c('name' = 'clean.name'))
phyto <- filter(phyto, !is.na(species))
phyto <- phyto %>% group_by(Lake_ID, name) %>% summarize(ind.size = mean(ind.size)) %>% ungroup
phyto <- phytoT %>% select(group, clean.name) %>% right_join(phyto, by = c('clean.name' = 'name'))
phyto <- rename(phyto, species = clean.name)
phyto$animal <- 'no'
phyto$species <- str_replace(phyto$species, ' ', '_')
phyto <- select(phyto, Lake_ID, species, group, animal, ind.size)

###### adding NLA #####

nla_phyto <- read_csv('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/NLA/2012/nla2012_wide_phytoplankton_count_02122014.csv')
colnames(nla_phyto) <- tolower(colnames(nla_phyto))
nla_phyto$species <- Hmisc::capitalize(tolower(nla_phyto$species))
nla_phyto$species <- str_replace(nla_phyto$species, ' ', '_')
nla_phyto <- nla_phyto %>% drop_na(species)

# distinct(nla_phyto, species) %>% pull(species) -> nla_list
# nla_list %in% taxlist
# sum(taxlist %in% nla_list)/142
# #44% of LP species present in NLA dataset

nla_phyto$ug.per.cell <- (nla_phyto$biovolume/nla_phyto$density) * 0.000000000001 * 1000000
# nla_phyto %>% filter(species == 'Ceratium_hirundinella') %>% pull(ug.per.cell)
# phyto %>% filter(str_detect(species, 'Ceratium')) #comparable

nla_phyto <- nla_phyto %>% select(site_id, species, algal_group, ug.per.cell) %>%
  mutate('animal' = 'no') %>% rename(Lake_ID = site_id, group = algal_group, ind.size = ug.per.cell) %>%
  select(Lake_ID, species, group, animal, ind.size)

# nla_phyto07 <- read_csv('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/NLA/2007/nla2007_phytoplankton_softalgaecount_20091023.csv')
# nla_phyto07 %>% distinct(TAXANAME) %>% pull(TAXANAME)
# nla_phyto07 %>% filter(TAXANAME == 'Kirchneriella lunaris') %>% pull(CELL_VOLUME)
# nla_phyto07 %>% filter(TAXANAME == 'Coelastrum microporum') %>% pull(CELL_VOLUME)
# #no measurement of cell size: all cells of a any given species are given the same cell volume

nla_zoo <- read_csv('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/NLA/2012/nla2012_zoopcnt_04032014.csv')
colnames(nla_zoo) <- tolower(colnames(nla_zoo))
nla_zoo$ug.per.ind <- with(nla_zoo, biomass/density)
nla_zoo <- select(nla_zoo, site_id, taxa_id, ug.per.ind)
nla_zoo <- nla_zoo %>% rename(Lake_ID = site_id)

nla_zoolist <- read_csv('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/NLA/2012/nla2012_zooptaxa_wide_10272015.csv')
colnames(nla_zoolist) <- tolower(colnames(nla_zoolist))
nla_zoolist$species <- tolower(nla_zoolist$species)
nla_zoolist <- nla_zoolist %>% drop_na(species)
nla_zoolist$genus <- Hmisc::capitalize(tolower(nla_zoolist$genus))
nla_zoolist$species <- paste(nla_zoolist$genus,nla_zoolist$species,sep='_')
nla_zoolist$subclass <- tolower(nla_zoolist$subclass)
nla_zoolist <- nla_zoolist %>% rename(group = subclass)
nla_zoolist <- nla_zoolist %>% select(taxa_id, species, group)

nla_zoo <- inner_join(nla_zoo, nla_zoolist, by = 'taxa_id')
nla_zoo$animal <- 'yes'
nla_zoo$ind.size <- nla_zoo$ug.per.ind
nla_zoo <- nla_zoo %>% select(Lake_ID, species, group, animal, ind.size)

nla_meta <- read_csv('~/Desktop/NLA/2012/nla2012_wide_siteinfo_08232016.csv') %>%
  rename(Lake_ID = SITE_ID, latitude = LAT_DD83) %>% select(Lake_ID, latitude) %>% as.data.frame

nla_data <- bind_rows(nla_phyto, nla_zoo) %>% arrange(Lake_ID)

nla_data$latitude <- nla_meta$latitude[match(nla_data$Lake_ID,nla_meta$Lake_ID)]
nla_data <- nla_data %>% drop_na(latitude, ind.size)

#### merge both dataframes and format ####

#only zoo
data <- zoo %>% arrange(Lake_ID, ind.size)

data <- bind_rows(zoo,phyto) %>% arrange(Lake_ID, ind.size)
data <- inner_join(data, lat, by = 'Lake_ID')
#data <- left_join(data, rbr, by = 'Lake_ID')

data <- bind_rows(data, nla_data) %>% arrange(species, ind.size)

## which variable will be used in analysis?
data$x <- data$latitude
xvarlab <- 'latitude (degrees)'
xvarlab2 <- 'latitudinal range (degrees)'
# data$x <- data$watertemp
# xvarlab <- 'water temperature (°C)'
# xvarlab2 <- 'thermal range (°C)'

##label for size measurement
sizelab <- expression(mean~body~size~(µm))
#sizelab <- expression(mean~body~mass~(µg))

data <- data %>% filter(!is.na(x))
data <- data %>% add_count(species)

#add scaled variables within species
taxlist <- data %>% distinct(species) %>% pull(species)
data$scl.size <- 0
data$scl.lat <- 0
for(i in 1:length(taxlist)){
  tmp <- filter(data, species == taxlist[i])
  scl.size <- scale(tmp$ind.size, scale = T)[,1]
  scl.lat <- scale(tmp$latitude, scale = F)[,1]
  data[data$species == taxlist[i], 'scl.size'] <- scl.size
  data[data$species == taxlist[i], 'scl.lat'] <- scl.lat
}

data$species <- as.factor(data$species)
data$scl.size[!is.finite(data$scl.size)] <- NA
data <- data %>% filter(ind.size > 0)

#### analysis. approach 1: fit a linear regression to each species with 10+ data points, and check distribution of slopes ####

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
  if(nrow(spdat) >= 10){
  #lmmod <- lm(ind.size~x,spdat)
  lmmod <- lm(scl.size~x,spdat)
  r2 <- summary(lmmod)$r.squared
  slope <- coef(lmmod)[2]
  p <- summary(lmmod)$coefficients[2,4]
  SE <- summary(lmmod)$coefficients[2,2]
  n <- nrow(spdat)
  min.lat <- min(spdat$x)
  max.lat <- max(spdat$x)
  group <- spdat$group[1]
  mean.size <- mean(spdat$ind.size)
  results1 <- c(focalsp,group)
  results2 <- c(min.lat,max.lat,n,slope,SE,p,r2,mean.size)
  reg.results[nrow(reg.results)+1,1:2] <- results1
  reg.results[nrow(reg.results),3:10] <- results2
  }
}

#reg.results <- filter(reg.results, abs(slope) < 100)
reg.results$xrange <- reg.results$max.lat-reg.results$min.lat
#reg.results <- reg.results %>% filter(xrange > 0.8)

reg.results$bubblesize <- rescale(reg.results$r2, c(1,4))
reg.results$bubblecol <- 'grey50'
reg.results$bubblecol[reg.results$p < 0.05 & reg.results$slope < 0] <- 'red'
reg.results$bubblecol[reg.results$p < 0.05 & reg.results$slope > 0] <- 'blue'
reg.results$bubble.pch <- 16
reg.results$bubble.pch[reg.results$bubblecol == 'grey50'] <- 1

hist(reg.results$slope, breaks=50)

#plot à la Dornelas
plot(r2~slope,reg.results,bty='n',pch=16,col=alpha(bubblecol,0.5),cex=bubblesize)
abline(v=0,lty=2)
#not very exciting...

#####

pdf('~/Desktop/regressions.pdf',width=8,height=7,pointsize=14)

labels <- parse(text=paste(10, '^', seq(-8,2,1), sep=''))

emptyPlot(xlim = range(data$x),yaxt='n',xaxt='n',ann=F, ylim=range(data$ind.size),bty='l',log='y')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=10^(-8:2),labels = labels)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab=xvarlab)
title(ylab=sizelab,line=2.8)
for(i in 1:nlevels(data$species)){
  focalsp <- levels(data$species)[i]
  spdat <- data %>% filter(species == focalsp) %>% arrange(x)
  spdat <- filter(spdat, !is.na(scl.size))
  if(nrow(spdat) >= 10){
    lmmod <- lm(scl.size~x,spdat)
    lncol <- 'grey50'
    lnwd <- 1
    if(summary(lmmod)$coefficients[2,4] < 0.05){
      lnwd <- 1
      if(summary(lmmod)$coefficients[2,1] > 0){lncol <- 'blue'}else{
        lncol <- 'red'}
    }
    lmmod <- lm(ind.size~x,spdat)
    points(fitted(lmmod)~spdat$x,col=alpha(lncol,0.5),type='l',lwd=lnwd)
  }
}

dev.off()

##### slopes  #####

pdf('~/Desktop/slopes.pdf',width=14,height=6,pointsize=14)
par(mfrow=c(1,3),cex=1,mar=c(4,2.5,0,0),oma=c(0,2,0,0))

plot(slope~n,reg.results,bty='n',pch=bubble.pch,col=alpha(bubblecol,0.5),cex=bubblesize,bty='l',yaxt='n',xaxt='n',ann=F, log = 'x')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(10,100,1000))
title(xlab=number~of~lakes,line=2.8)
abline(h=0,lty=3)

mtext(expression(Delta~body~mass~per~degree~latitude~(sd)),outer = T,side=2,line=0.5)

# which.max(reg.results$n)
# reg.results[244,]
# reg.results %>% arrange(n) %>% pull(focalsp)

plot(slope~xrange,reg.results,bty='n',pch=bubble.pch,col=alpha(bubblecol,0.5),cex=bubblesize,bty='l',yaxt='n',xaxt='n',ann=F,xlim=(range(reg.results$xrange) + c(-0.5,0.5)))
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab=xvarlab2,line=2.8)
#title(ylab=slope~(beta),line=2.8)
abline(h=0,lty=3)

plot(slope~mean.size,reg.results,bty='n',pch=bubble.pch,col=alpha(bubblecol,0.5),cex=bubblesize,bty='l',log='x',yaxt='n',xaxt='n',ann=F)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=10^(-8:2),labels = labels)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
#title(ylab=slope~(beta),line=2.8)
title(xlab=sizelab,line=2.8)
abline(h=0,lty=3)

dev.off()

##### #####

sum(reg.results$bubblecol == 'blue') / nrow(reg.results)
#20% bergmann cline
sum(reg.results$bubblecol == 'red') / nrow(reg.results)
#5% inverse bergmann
sum(reg.results$bubblecol == 'grey50') / nrow(reg.results)
#75% no trend

###### models to test if there is an overall trend #####

#mixed model

mod <- lmer(log(ind.size) ~ x + (x-1|species) + (1|species), data = data)
plot(mod)
summary(mod)
confint(mod)
hist(resid(mod))
performance(mod)
icc(mod)
coefplot2(mod)

#bayesian meta-analytical model

library(MCMCglmm)
#MCMCglmm parameters
a <- 1000
priorz <- list(R = list(V = diag(1), nu = 0.002), G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))
prior.s <- list(R=list(V=diag(4),n=0.002),G=list(G1=list(V=diag(1),n=0.002),G2=list(V=diag(1),n=0.002)))
m <-reg.results$SE^2 #variance of coefficients
bLMM <- MCMCglmm(slope ~ 1, data=reg.results, mev=m, nitt=103000, verbose=T)
summary(bLMM)

######

# intra vs. interspecific contribution to variation in community-weighted mean body size



##### GAMM #####

mod <- bam(ind.size ~ s(x, k = 8) + s(x, species, bs = 'fs', k = 6), data = data,nthreads=2)
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
for(i in 1:nlevels(data$species)){
  focalsp <- levels(data$species)[i]
  spdat <- data %>% filter(species == focalsp)
  spdat <- arrange(spdat, x)
  gammod <- gam(ind.size~s(x, k = 4),data=spdat)
  points(fitted(gammod)~spdat$x,col=alpha(cols2[i],0.5),type='l',lwd=0.8)
}
plot_smooth(mod, view="x", lwd=3, col=1, rm.ranef=T, se=1.96, rug=F, add=T)
legend('topright',bty='n',legend=bquote(atop(italic('F') == .(testres[1]),italic('p') == .(testres[2]))))


#######

size.diffs <- numeric(0)
for(i in 1:nlevels(data$species)){
  focalsp <- levels(data$species)[i]
  spdat <- data %>% filter(species == focalsp)
  diff <- (max(spdat$ind.size) - min(spdat$ind.size))/min(spdat$ind.size)
  size.diffs[i] <-
}

####lake mean size###
zoo <- 

