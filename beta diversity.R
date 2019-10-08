rm(list=ls())

#libraries
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(adespatial)
library(devtools)
library(sp)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
source_url("https://raw.githubusercontent.com/VFugere/LakePulse/master/mybxp.R")

#cols
cols <- brewer.pal(7, 'Dark2')[1:3] 
cols2 <- brewer.pal(8, 'Dark2')[c(4,5,6,8)] 

#data
load('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/2017data.RData')

basic.data <- filter(basic.data, year == 2017)

plot(log1p(HI)~log(area),basic.data)
plot(HI~area,basic.data)
hist(basic.data$HI,breaks=100)

basic.data$hi_class2 <- 'low'
basic.data$hi_class2[basic.data$HI > 0.05] <- 'moderate'
basic.data$hi_class2[basic.data$HI > 0.4] <- 'high'

trt <- basic.data %>% select(Lake_ID, hi_class2, size_class, ecozone, HI, area, latitude, longitude) %>%
  rename(hi_cont = HI, area_cont = area, lat=latitude,long=longitude) %>%
  rename(area = size_class, HI = hi_class2)

#### community matrix ####

com <- phyto

meta <- com %>% select(Lake_ID)
meta <- left_join(meta, trt, by = 'Lake_ID')
meta$HI_Ez <- paste(meta$ecozone,meta$HI,sep='_')
meta <- mutate_at(meta, vars(HI:HI_Ez), as.factor)

meta$HI <- factor(meta$HI, levels=c('low','moderate','high'))
meta$area <- factor(meta$area, levels=c('small','medium','large'))
meta$ecozone <- as.factor(meta$ecozone)
levels(meta$ecozone) <- c('AH','AM','BS','MP')

com <- com[,2:ncol(com)] %>% as.matrix
com[is.na(com)] <- 0

row.names(com) <- 1:nrow(com)
row.names(meta) <- 1:nrow(meta)

# dm <- vegdist(log1p(com),method = 'bray')
# adonis(dm~meta$HI*meta$ecozone*meta$area)
# adonis(dm~meta$HI+meta$ecozone+meta$area+meta$HI:meta$ecozone)
# #area, ecozone, HI have a signigicant impact on composition.
# #significant HI by ecozone interaction

# #for bacteria, say I only want to use abundant ASVs. Does not change anything.
# to.rm <- as.numeric(which(colSums(com) < 1000))
# com <- com[,-(to.rm)]

##### including all sites together ####

transform.names <- c('no transformation','presence/absence','proportional abundance',
                     'Hellinger transform','Log 1+x','Square root','Wisconsin double')

pdf('~/Desktop/bd_phyto.pdf',width=5,height = 6,pointsize = 8, onefile=T)
par(mfrow=c(3,2),cex=1,mar=c(3,4,1,1),oma=c(0,0,2,0))

for(t in 1:7){
  
  #distance matrix
  com.sub <- com
  meta.sub <- meta
  
  if(t == 2){com.sub[com.sub > 0] <- 1}
  if(t == 3){com.sub <- decostand(com.sub,method='total')}
  if(t == 4){com.sub <- decostand(com.sub,method='hellinger')}
  if(t == 5){com.sub <- log1p(com.sub)}
  if(t == 6){com.sub <- sqrt(com.sub)}
  if(t == 7){com.sub <- wisconsin(com.sub)}
  
  dm <- vegdist(com.sub,method = 'bray')
  
  #HI impact
  cc <- adonis(dm~meta.sub$HI)
  p.cc <- round(cc$aov.tab$`Pr(>F)`,2)[1]
  bd <- anova(betadisper(d=dm, group=meta.sub$HI, type='centroid'))
  p.bd <- round(bd$`Pr(>F)`[1],2)
  plot(betadisper(d=dm, group=meta.sub$HI, type='centroid'),ellipse=F,hull=F,label=F,seg.col=alpha(cols,0.1),segments=T,col=cols,pch=rep(1,3),cex=0.1,main=NULL,sub=NULL,ylab=NULL,xlab=NULL)
  mtext(bquote('Community composition,'~italic(p) == .(p.cc)),side=3,adj=0,cex=1)
  legend('topright',bty='n',legend=levels(meta.sub$HI),pch=16,col=cols)
  distances <- (betadisper(d=dm, group=meta.sub$HI, type='centroid'))$distances
  mybxp(y=distances,x=meta.sub$HI,colvec=cols,ylab='distance to centroid')
  mtext(bquote(italic(beta)~'diversity, '~italic(p) == .(p.bd)),side=3,adj=0,cex=1)
  
  #Area impact
  cc <- adonis(dm~meta.sub$area)
  p.cc <- round(cc$aov.tab$`Pr(>F)`,4)[1]
  bd <- anova(betadisper(d=dm, group=meta.sub$area, type='centroid'))
  p.bd <- round(bd$`Pr(>F)`[1],4)
  plot(betadisper(d=dm, group=meta.sub$area, type='centroid'),ellipse=F,hull=F,label=F,seg.col=alpha(cols,0.1),segments=T,col=cols,pch=rep(1,3),cex=0.1,main=NULL,sub=NULL,ylab=NULL,xlab=NULL)
  mtext(bquote('Community composition,'~italic(p) == .(p.cc)),side=3,adj=0,cex=1)
  legend('topright',bty='n',legend=levels(meta.sub$area),pch=16,col=cols)
  distances <- (betadisper(d=dm, group=meta.sub$area, type='centroid'))$distances
  mybxp(y=distances,x=meta.sub$area,colvec=cols,ylab='distance to centroid')
  mtext(bquote(italic(beta)~'diversity, '~italic(p) == .(p.bd)),side=3,adj=0,cex=1)
  
  #Ecozone
  cc <- adonis(dm~meta.sub$ecozone)
  p.cc <- round(cc$aov.tab$`Pr(>F)`,4)[1]
  bd <- anova(betadisper(d=dm, group=meta.sub$ecozone, type='centroid'))
  p.bd <- round(bd$`Pr(>F)`[1],4)
  plot(betadisper(d=dm, group=meta.sub$ecozone, type='centroid'),ellipse=F,hull=F,label=F,seg.col=alpha(cols2,0.1),segments=T,col=cols2,pch=rep(1,4),cex=0.1,main=NULL,sub=NULL,ylab=NULL,xlab=NULL)
  mtext(bquote('Community composition,'~italic(p) == .(p.cc)),side=3,adj=0,cex=1)
  legend('topright',bty='n',legend=levels(meta.sub$ecozone),pch=16,col=cols2)
  distances <- (betadisper(d=dm, group=meta.sub$ecozone, type='centroid'))$distances
  mybxp(y=distances,x=meta.sub$ecozone,colvec=cols2,ylab='distance to centroid')
  mtext(bquote(italic(beta)~'diversity, '~italic(p) == .(p.bd)),side=3,adj=0,cex=1)
  
  mtext(bquote(bold(.(transform.names[t]))),side=3,outer=T,adj=0.5,cex=1)
  
}

dev.off()

##### beta diversity partitioning #####

pdf('~/Desktop/betadivpart.pdf',width=3,height=4.5,pointsize = 8)

par(mfrow=c(2,1),mar=c(4,2,3,1),oma=c(0,4,0,0),cex=1)

z<-beta.div.comp(zoo.biomass.grouped[,2:ncol(zoo.biomass.grouped)], coef='J', quant=F)
p<-beta.div.comp(phyto[,2:ncol(phyto)], coef='J', quant=F)
b<-beta.div.comp(bacterio[,2:ncol(bacterio)], coef='J', quant=F)
barplot(cbind(b$part[4:5],p$part[4:5],z$part[4:5]),names.arg=c('bacterio','phyto','zoo'),ylim=c(0,1),xlim=c(0,1.75),width=0.5,cex.lab=1,cex.axis=1,cex.names=1,col=c('navyblue','gray82'),ylab=NULL)
title(xlab='species occurence',cex=1.2)
legend(horiz=T,'topleft',inset=c(0.1,-0.2),xpd=T,cex=1, pt.cex=1, legend=c(expression(paste(Delta,' richness',sep='')),'replacement'),pch=22,pt.bg=c('gray82','navyblue'),bty='n')

z<-beta.div.comp(zoo.biomass.grouped[,2:ncol(zoo.biomass.grouped)], coef='J', quant=T)
p<-beta.div.comp(phyto[,2:ncol(phyto)], coef='J', quant=T)
b<-beta.div.comp(bacterio[,2:ncol(bacterio)], coef='J', quant=T)
barplot(cbind(b$part[4:5],p$part[4:5],z$part[4:5]),names.arg=c('bacterio','phyto','zoo'),ylim=c(0,1),xlim=c(0,1.75),width=0.5,cex.lab=1,cex.axis=1,cex.names=1,col=c('navyblue','gray82'),ylab=NULL)
title(xlab='species abundance',cex=1.2)
legend(horiz=T,'topleft',inset=c(0.1,-0.2),xpd=T,cex=1, pt.cex=1, legend=c(expression(paste(Delta,' abundance',sep='')),'replacement'),pch=22,pt.bg=c('gray82','navyblue'),bty='n')

mtext(expression('contribution to'~italic(beta)~diversity~'(proportion)'),side=2,outer=T,line=2,cex=1.2)
dev.off()

# #does not change anything
# z<-beta.div.comp(zoo.biomass.grouped[,2:ncol(zoo.biomass.grouped)], coef='S', quant=F)
# p<-beta.div.comp(phyto[,2:ncol(phyto)], coef='S', quant=F)
# b<-beta.div.comp(bacterio[,2:ncol(bacterio)], coef='S', quant=F)
# barplot(cbind(b$part[4:5],p$part[4:5],z$part[4:5]),names.arg=c('bacterio','phyto','zoo'),ylim=c(0,1),xlim=c(0,2),width=0.5,cex.lab=0.7,cex.axis=0.7,cex.names=0.7,col=c('navyblue','gray82'),ylab='relative contribution to total dissimilarity')
# legend(horiz=T,'topleft',inset=c(0.1,-0.1),xpd=T,cex=0.7, pt.cex=1, legend=c(expression(paste(Delta,' richness',sep='')),'replacement'),pch=22,pt.bg=c('gray82','navyblue'),bty='n')
# 
# #does not change anything
# z<-beta.div.comp(zoo.biomass.grouped[,2:ncol(zoo.biomass.grouped)], coef='S', quant=T)
# p<-beta.div.comp(phyto[,2:ncol(phyto)], coef='S', quant=T)
# b<-beta.div.comp(bacterio[,2:ncol(bacterio)], coef='S', quant=T)
# barplot(cbind(b$part[4:5],p$part[4:5],z$part[4:5]),names.arg=c('bacterio','phyto','zoo'),ylim=c(0,1),xlim=c(0,2),width=0.5,cex.lab=0.7,cex.axis=0.7,cex.names=0.7,col=c('navyblue','gray82'),ylab='relative contribution to total dissimilarity')
# legend(horiz=T,'topleft',inset=c(0.1,-0.1),xpd=T,cex=0.7, pt.cex=1, legend=c(expression(paste(Delta,' abundance',sep='')),'replacement'),pch=22,pt.bg=c('gray82','navyblue'),bty='n')

#split by factor levels

z.trt <- zoo.biomass.grouped %>% select(Lake_ID) %>% left_join(trt, by='Lake_ID')
p.trt <- phyto %>% select(Lake_ID) %>% left_join(trt, by='Lake_ID')
b.trt <- bacterio %>% select(Lake_ID) %>% left_join(trt, by='Lake_ID')

pdf('~/Desktop/betadivpart_levels.pdf',width=7,height=4,pointsize = 10)
par(mfrow=c(2,3),mar=c(4,2,1,1),oma=c(0,4,0,0),cex=1)

sizes <- c('small','medium','large')

for(i in 1:3){
  lev <- sizes[i]
  z<-beta.div.comp(zoo.biomass.grouped[z.trt$area == lev,2:ncol(zoo.biomass.grouped)], coef='J', quant=F)
  p<-beta.div.comp(phyto[p.trt$area == lev,2:ncol(phyto)], coef='J', quant=F)
  b<-beta.div.comp(bacterio[b.trt$area == lev,2:ncol(bacterio)], coef='J', quant=F)
  barplot(cbind(b$part[4:5],p$part[4:5],z$part[4:5]),names.arg=c('bacterio','phyto','zoo'),ylim=c(0,1),xlim=c(0,1.75),width=0.5,cex.lab=1,cex.axis=1,cex.names=0.8,col=c(cols[i],'gray82'),ylab=NULL)
  #if(i==1){legend(horiz=T,'topleft',inset=c(0.1,-0.2),xpd=T,cex=1, pt.cex=1, legend=c(expression(paste(Delta,' richness',sep='')),'replacement'),pch=22,pt.bg=c('gray82','navyblue'),bty='n')}
  title(xlab=lev,cex=1.2)
}


#HI

HI <- c('low','moderate','high')


for(i in 1:3){
  lev <- HI[i]
  z<-beta.div.comp(zoo.biomass.grouped[z.trt$HI == lev,2:ncol(zoo.biomass.grouped)], coef='J', quant=F)
  p<-beta.div.comp(phyto[p.trt$HI == lev,2:ncol(phyto)], coef='J', quant=F)
  b<-beta.div.comp(bacterio[b.trt$HI == lev,2:ncol(bacterio)], coef='J', quant=F)
  barplot(cbind(b$part[4:5],p$part[4:5],z$part[4:5]),names.arg=c('bacterio','phyto','zoo'),ylim=c(0,1),xlim=c(0,1.75),width=0.5,cex.lab=1,cex.axis=1,cex.names=0.8,col=c(cols[i],'gray82'),ylab=NULL)
  #if(i==1){legend(horiz=T,'topleft',inset=c(0.1,-0.15),xpd=T,cex=0.9, pt.cex=1, legend=c(expression(paste(Delta,' richness',sep='')),'replacement'),pch=22,pt.bg=c('gray82','navyblue'),bty='n')}
  title(xlab=lev,cex=1.2)
}

mtext(expression('contribution to'~italic(beta)~diversity~'(proportion)'),side=2,outer=T, cex=1.2, line=2)
dev.off()

###### Turnover along continuous gradient #####

z.trt <- zoo.biomass.grouped %>% select(Lake_ID) %>% left_join(trt, by='Lake_ID')
p.trt <- phyto %>% select(Lake_ID) %>% left_join(trt, by='Lake_ID')
b.trt <- bacterio %>% select(Lake_ID) %>% left_join(trt, by='Lake_ID')

pdf('~/Desktop/distance_decay_Jaccard.pdf', height = 6,width=5,pointsize = 8)
par(mfrow=c(3,2),cex=1,mar=c(4,4,1,1),oma=c(0,0,0,0))

com <- bacterio
meta <- b.trt
com <- com[,2:ncol(com)] %>% as.matrix
dist.hi <- as.numeric(dist(meta$hi_cont))
dist.spatial <- spDists(as.matrix(meta[,c('long','lat')]),longlat = T)
dist.spatial <- as.numeric(dist.spatial[lower.tri(dist.spatial, diag=F)])
sim.cc <- 1-(as.numeric(vegdist(com,method='jaccard',binary=T)))
plot(log1p(sim.cc)~dist.hi,pch=16,cex=0.2,bty='l',col=alpha(1,0.3),ylab='community similarity',xlab='HI difference')
mod1 <- lm(log1p(sim.cc)~dist.hi)
legend('topright',bty='n',legend=bquote(R^2 == .(round(summary(mod1)$r.squared,2))))
abline(mod1,lwd=3,col='navy blue')
legend('topleft',bty='n',legend=expression(italic(bacterio)))
plot(log1p(sim.cc)~dist.spatial,pch=16,cex=0.2,bty='l',col=alpha(1,0.3),ylab='community similarity',xlab='spatial distance (km)')
mod2 <- lm(log1p(sim.cc)~dist.spatial)
legend('topright',bty='n',legend=bquote(R^2 == .(round(summary(mod2)$r.squared,2))))
abline(mod2,lwd=3,col='dark orange')

com <- phyto
meta <- p.trt
com <- com[,2:ncol(com)] %>% as.matrix
dist.hi <- as.numeric(dist(meta$hi_cont))
dist.spatial <- spDists(as.matrix(meta[,c('long','lat')]),longlat = T)
dist.spatial <- as.numeric(dist.spatial[lower.tri(dist.spatial, diag=F)])
sim.cc <- 1-(as.numeric(vegdist(com,method='jaccard',binary=T)))
plot(log1p(sim.cc)~dist.hi,pch=16,cex=0.2,bty='l',col=alpha(1,0.3),ylab='community similarity',xlab='HI difference')
mod1 <- lm(log1p(sim.cc)~dist.hi)
legend('topright',bty='n',legend=bquote(R^2 == .(round(summary(mod1)$r.squared,2))))
abline(mod1,lwd=3,col='navy blue')
legend('topleft',bty='n',legend=expression(italic(phyto)))
plot(log1p(sim.cc)~dist.spatial,pch=16,cex=0.2,bty='l',col=alpha(1,0.3),ylab='community similarity',xlab='spatial distance (km)')
mod2 <- lm(log1p(sim.cc)~dist.spatial)
legend('topright',bty='n',legend=bquote(R^2 == .(round(summary(mod2)$r.squared,2))))
abline(mod2,lwd=3,col='dark orange')

com <- zoo.biomass.grouped
meta <- z.trt
com <- com[,2:ncol(com)] %>% as.matrix
dist.hi <- as.numeric(dist(meta$hi_cont))
dist.spatial <- spDists(as.matrix(meta[,c('long','lat')]),longlat = T)
dist.spatial <- as.numeric(dist.spatial[lower.tri(dist.spatial, diag=F)])
sim.cc <- 1-(as.numeric(vegdist(com,method='jaccard',binary=T)))
plot(log1p(sim.cc)~dist.hi,pch=16,cex=0.2,bty='l',col=alpha(1,0.3),ylab='community similarity',xlab='HI difference')
mod1 <- lm(log1p(sim.cc)~dist.hi)
legend('topright',bty='n',legend=bquote(R^2 == .(round(summary(mod1)$r.squared,2))))
abline(mod1,lwd=3,col='navy blue')
legend('topleft',bty='n',legend=expression(italic(zoo)))
plot(log1p(sim.cc)~dist.spatial,pch=16,cex=0.2,bty='l',col=alpha(1,0.3),ylab='community similarity',xlab='spatial distance (km)')
mod2 <- lm(log1p(sim.cc)~dist.spatial)
legend('topright',bty='n',legend=bquote(R^2 == .(round(summary(mod2)$r.squared,2))))
abline(mod2,lwd=3,col='dark orange')

dev.off()

##

pdf('~/Desktop/distance_decay_Bray_logged.pdf', height = 6,width=5,pointsize = 8)
par(mfrow=c(3,2),cex=1,mar=c(4,4,1,1),oma=c(0,0,0,0))

com <- bacterio
meta <- b.trt
com <- com[,2:ncol(com)] %>% as.matrix
dist.hi <- as.numeric(dist(meta$hi_cont))
dist.spatial <- spDists(as.matrix(meta[,c('long','lat')]),longlat = T)
dist.spatial <- as.numeric(dist.spatial[lower.tri(dist.spatial, diag=F)])
sim.cc <- 1-(as.numeric(vegdist(log1p(com),method='bray',binary=F)))
plot(log1p(sim.cc)~dist.hi,pch=16,cex=0.2,bty='l',col=alpha(1,0.3),ylab='community similarity',xlab='HI difference')
mod1 <- lm(log1p(sim.cc)~dist.hi)
legend('topright',bty='n',legend=bquote(R^2 == .(round(summary(mod1)$r.squared,2))))
abline(mod1,lwd=3,col='navy blue')
legend('topleft',bty='n',legend=expression(italic(bacterio)))
plot(log1p(sim.cc)~dist.spatial,pch=16,cex=0.2,bty='l',col=alpha(1,0.3),ylab='community similarity',xlab='spatial distance (km)')
mod2 <- lm(log1p(sim.cc)~dist.spatial)
legend('topright',bty='n',legend=bquote(R^2 == .(round(summary(mod2)$r.squared,2))))
abline(mod2,lwd=3,col='dark orange')

com <- phyto
meta <- p.trt
com <- com[,2:ncol(com)] %>% as.matrix
dist.hi <- as.numeric(dist(meta$hi_cont))
dist.spatial <- spDists(as.matrix(meta[,c('long','lat')]),longlat = T)
dist.spatial <- as.numeric(dist.spatial[lower.tri(dist.spatial, diag=F)])
sim.cc <- 1-(as.numeric(vegdist(log1p(com),method='bray',binary=F)))
plot(log1p(sim.cc)~dist.hi,pch=16,cex=0.2,bty='l',col=alpha(1,0.3),ylab='community similarity',xlab='HI difference')
mod1 <- lm(log1p(sim.cc)~dist.hi)
legend('topright',bty='n',legend=bquote(R^2 == .(round(summary(mod1)$r.squared,2))))
abline(mod1,lwd=3,col='navy blue')
legend('topleft',bty='n',legend=expression(italic(phyto)))
plot(log1p(sim.cc)~dist.spatial,pch=16,cex=0.2,bty='l',col=alpha(1,0.3),ylab='community similarity',xlab='spatial distance (km)')
mod2 <- lm(log1p(sim.cc)~dist.spatial)
legend('topright',bty='n',legend=bquote(R^2 == .(round(summary(mod2)$r.squared,2))))
abline(mod2,lwd=3,col='dark orange')

com <- zoo.biomass.grouped
meta <- z.trt
com <- com[,2:ncol(com)] %>% as.matrix
dist.hi <- as.numeric(dist(meta$hi_cont))
dist.spatial <- spDists(as.matrix(meta[,c('long','lat')]),longlat = T)
dist.spatial <- as.numeric(dist.spatial[lower.tri(dist.spatial, diag=F)])
sim.cc <- 1-(as.numeric(vegdist(log1p(com),method='bray',binary=F)))
plot(log1p(sim.cc)~dist.hi,pch=16,cex=0.2,bty='l',col=alpha(1,0.3),ylab='community similarity',xlab='HI difference')
mod1 <- lm(log1p(sim.cc)~dist.hi)
legend('topright',bty='n',legend=bquote(R^2 == .(round(summary(mod1)$r.squared,2))))
abline(mod1,lwd=3,col='navy blue')
legend('topleft',bty='n',legend=expression(italic(zoo)))
plot(log1p(sim.cc)~dist.spatial,pch=16,cex=0.2,bty='l',col=alpha(1,0.3),ylab='community similarity',xlab='spatial distance (km)')
mod2 <- lm(log1p(sim.cc)~dist.spatial)
legend('topright',bty='n',legend=bquote(R^2 == .(round(summary(mod2)$r.squared,2))))
abline(mod2,lwd=3,col='dark orange')

dev.off()