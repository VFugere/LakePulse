rm(list=ls())

#libraries
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(adespatial)
library(devtools)
library(sp)
library(ape)
library(party)
library(rworldmap)
library(gllvm)

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

#remove lakes not in qc or ontario

lakes.to.keep <- basic.data %>% filter(province %in% c('QUEBEC','ONTARIO')) %>% pull(Lake_ID)

basic.data <- filter(basic.data, Lake_ID %in% lakes.to.keep)
phyto <- filter(phyto, Lake_ID %in% lakes.to.keep)
bacterio <- filter(bacterio, Lake_ID %in% lakes.to.keep)
zoo <- zoo.biomass.grouped %>% filter(Lake_ID %in% lakes.to.keep)

fish <- as.data.frame(fish)

#### figure ####

#pdf('~/Desktop/betadiv1.pdf',width=13,height = 7,pointsize = 12)
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

#dev.off()

#### LVM for GRIL 2020 poster ####

fish.lakes <- fish$Lake_ID

env <- basic.data %>% select(Lake_ID, province, latitude, longitude, area, Elevation, depth_m, HI, Shore_dev) %>%
  filter(Lake_ID %in% fish.lakes) %>%
  rename(elevation = Elevation, depth = depth_m, SD = Shore_dev)

zoo.com <- zoo[,2:ncol(zoo)]
zoo.pa <- zoo.com
zoo.pa[zoo.pa > 0] <- 1
zoo.size <- zoo.pa
bad.sp <- names(zoo.com)[names(zoo.com) %!in% zooT$taxon]
bad.sp <- c(bad.sp,zooT$taxon[is.na(zooT$length)])
zoo.size <- zoo.size[,(names(zoo.size) %!in% bad.sp)]
for(i in 1:ncol(zoo.size)){
  col.sub <- zoo.size[,i]
  sp <- names(zoo.size)[i]
  sp.len <- zooT$length[zooT$taxon == sp]
  zoo.size[,i] <- zoo.size[,i] * sp.len
}

#zoo.ordi<-metaMDS(comm = zoo.pa,distance = 'jaccard', k=2, trymax = 20)
dm <- vegdist(log1p(zoo.com),method = 'bray')
ordi <- betadisper(d=dm, group=rep(1,nrow(zoo.com)), type='centroid')
zoo.cc <- scores(ordi)$sites

zoo.df <- data.frame(zoo$Lake_ID)
zoo.df$bm <- apply(zoo.com, 1, sum)
zoo.df$max.size <- apply(zoo.size, 1, max)
zoo.df <- cbind(zoo.df, zoo.cc)

zoo.df <- filter(zoo.df, zoo.Lake_ID %in% fish.lakes)
fish <- filter(fish, Lake_ID %in% zoo.df$zoo.Lake_ID)
env <- filter(env, Lake_ID %in% zoo.df$zoo.Lake_ID)

zoo.df <- zoo.df %>% rename(Lake_ID = zoo.Lake_ID) %>% arrange(Lake_ID)
fish <- arrange(fish, Lake_ID)
env <- arrange(env, Lake_ID)

env2 <- select(env, -Lake_ID)
env.zoo <- left_join(env,zoo.df,'Lake_ID') %>% select(-Lake_ID)

fish.com <- select(fish, -Lake_ID)
fish.com <- fish.com[,apply(fish.com, 2, sum) != 0]

## models at last

library(corrplot)
library(gclus)

fish.com -> backup

fish.com <- backup
common.tax <- apply(fish.com, 2, sum)
common.tax <- common.tax[common.tax > 5]
common.tax <- names(common.tax)
fish.com <- fish.com %>% select(common.tax)

fit_ord <- gllvm(fish.com, family='binomial', num.lv=2)
fit_ord
plot(fit_ord)
ordiplot(fit_ord, biplot=T, ind.spp = 15, symbols=T, pch=1)

cr <- getResidualCor(fit_ord)
pdf('~/Desktop/corrs_ord.pdf',width = 12,height = 12,pointsize = 18)
corrplot(cr[order.single(cr), order.single(cr)],diag=F,type='lower',method='square',tl.cex=1,tl.srt=5,tl.col = 1)
dev.off()

env.zoo$province <- as.numeric(env.zoo$province)
env.zoo <- mutate_all(env.zoo, scale)

fit_env <- gllvm(fish.com, env.zoo, family='binomial', num.lv = 2, formula = ~ latitude)
fit_env
fit_zoo <- gllvm(fish.com, env.zoo, family='binomial', num.lv = 2, formula = ~ latitude+longitude+area+depth+HI+max.size)
fit_zoo
summary(fit_ord)

cr <- getResidualCor(fit_env)
pdf('~/Desktop/corrs.pdf',width = 12,height = 12,pointsize = 8)
corrplot(cr[order.single(cr), order.single(cr)],diag=F,type='lower',method='square',tl.cex=1,tl.srt=45,tl.col = 1)
dev.off()

fit_env <- gllvm(fish.com, env.zoo, family='binomial', num.lv = 2, formula = ~ latitude+longitude+depth+latitude:longitude)
fit_env
fit_env2 <- gllvm(fish.com, env.zoo, family='binomial', num.lv = 2, formula = ~ latitude+longitude+depth+province)
fit_env2

fit_pr <- gllvm(fish.com, env.zoo, family='binomial', num.lv = 2, formula = ~ province)
fit_pr

pdf('~/Desktop/coefs.pdf',width = 3,height = 8,pointsize = 12)
coefplot(fit_env,cex.ylab=0.7,mar=c(4,9,2,1))
dev.off()


ordi <- metaMDS(fish.com,distance='jaccard',k=2)

g<-ordi$points[,1:2]

zvec <- scales::rescale(env$depth, to=c(0.5,2.5))

pdf('~/Desktop/ordiplot.pdf',width=8,height=8,pointsize=20)
plot(g[,2] ~ g[,1], type = "p",pch=16,col=alpha('black',0.3),cex=zvec,yaxt='n',xaxt='n',ann=F)
title(xlab='NMDS dimension 1',cex.lab=1,line = 2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='NMDS dimension 2',cex.lab=1,line = 2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
common.tax <- apply(fish.com, 2, sum)
common.tax <- common.tax[common.tax > 11]
common.tax <- names(common.tax)
labels <- c('AMRU','AMNE','CACO','ESLU','LEGI','MIDO','MISA','PEFL','SANA','SAVI')
label.subset <- ordi$species[rownames(ordi$species) %in% common.tax,]
text(label.subset, make.italic(labels), cex = 1, col = scales::alpha('blue',0.8))
legend('bottomleft',bty='n',legend='Stress = 0.11')
legend('topright',bty='n',legend=make.italic('Point size indicates lake depth'))
dev.off()

mapplot2 <- function(x,y,z,name){
  xrange <- range(x)+c(-0.5,0.5)
  yrange <- range(y)+c(-1,1)
  plot(map, xlim = xrange, ylim = yrange,col='light gray',border=0,asp=1.2,axes=F,cex.lab=0.5)
  points(x=x,y=y,pch=1,col=1,cex=1.2)
  points(x=x,y=y,pch=16,col=c('white','darkblue')[z+1],cex=1.1)
  legend('topleft',legend=make.italic(name),bty='n',text.col=1)
}

for(i in 1:ncol(fish.com)){
  sub <- fish.com[,i]
  mapplot2(env$longitude,env$latitude,sub,name=names(fish.com)[i])
}

cols <- c(0,RColorBrewer::brewer.pal(11,'RdYlBu')[c(1,6,11)])

fish.com[,c('Sander_vitreus','Salvelinus_namaycush')]

colvec<-c(1,3,4,2,3,2,4,2,2,2,2,2,1,2,2,4,2,2,4,4,4,4,1,4,4,2,2,4,4,4)

xrange <- range(env$longitude)+c(-0.5,0.5)
yrange <- range(env$latitude)+c(-1,1)
pdf('~/Desktop/fishmap.pdf',width=8,height = 5,pointsize = 18)
plot(map, xlim = xrange, ylim = yrange,col='light gray',border=0,asp=1.5,axes=F,cex.lab=0.5)
points(x=env$longitude,y=env$latitude,pch=1,col=1,cex=2)
points(x=env$longitude,y=env$latitude,pch=16,col=cols[colvec],cex=1.9)
#legend('topleft',bty='n',text.col=1,pch=16,ncol=2,legend=c('SAVI','SANA','both','neither'),col=cols[c(2,4,3,1)])
dev.off()
