rm(list=ls())

#libraries
library(tidyverse)
library(RColorBrewer)
library(rworldmap)
library(readxl)
library(scales)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
devtools::source_url("https://raw.githubusercontent.com/VFugere/LakePulse/master/custom_plots.R")

#cols
cols <- brewer.pal(3, 'Dark2')
cols2 <- brewer.pal(8, 'Dark2')[c(4,5,6,8)] 

#data
load('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/basic_data.RData')
# zoo <- zoo.biomass.grouped 
# basic.data <- filter(basic.data, year == 2017)
# basic.data$ecozone <- as.factor(basic.data$ecozone)
# levels(basic.data$ecozone) <- c('AH','AM','BS','MP')

#map
map <- getMap(resolution = "low")

# #LP sites with 3 overlapping datasets
# merged <- inner_join(basic.data,zoo)
# merged <- inner_join(basic.data, bacterio) %>% 
#   inner_join(phyto) %>% inner_join(zoo)
# name <- 'all (147/217 sites)'
# col <- 'goldenrod4'
# 
# x <- basic.data$longitude
# y <- basic.data$latitude
# xrange <- range(x)+c(-2,2)
# yrange <- range(y)+c(-1,1)
# plot(map, xlim = xrange, ylim = yrange,col='light gray',border=0,asp=1.2,axes=F,cex.lab=0.5)
# points(x=x,y=y,pch=1,col=1,cex=1.2)
# points(x=merged$longitude,y=merged$latitude,pch=16,col=alpha(col,0.5),cex=1.1)
# legend('topright',legend=make.italic(name),bty='n',text.col=1)
# 
# x <- basic.data$longitude
# y <- basic.data$latitude
# xrange <- range(x)+c(-2,2)
# yrange <- range(y)+c(-1,1)
# plot(map, xlim = xrange, ylim = yrange,col='light gray',border=0,asp=1.2,axes=F,cex.lab=0.5)
# points(x=x,y=y,pch=1,col=1,cex=1.2)
# points(x=basic.data$longitude,y=basic.data$latitude,pch=16,col=alpha(cols2[basic.data$ecozone],0.5),cex=1.1)
# legend('topright',legend=make.italic('2017 Lake Pulse campaign'),bty='n',text.col=1)
# legend('bottomright',bty='n',legend=levels(basic.data$ecozone),pch=15,col=cols2)
# 
# pdf('~/Desktop/lakes.pdf',width=8,height = 4,pointsize = 12)
# par(mfrow=c(1,2))
# plot(depth_m~area,basic.data,log='x',xlab='lake area',ylab='max depth',pch=16,col=alpha(cols2[basic.data$ecozone],0.5),bty ='l')
# plot(HI~area,basic.data,log='x',xlab='lake area',ylab='human impact index',pch=16,col=alpha(cols2[basic.data$ecozone],0.5),bty ='l')
# dev.off()

## fish

fish <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/fish_sites.xlsx')

basic.coords <- select(basic.data, longitude, latitude, Lake_ID)

fish <- left_join(fish, basic.coords, 'Lake_ID')

name <- 'fish data available'
colz <- '#183EA8'

x <- fish$longitude
y <- fish$latitude
xrange <- range(x)+c(-2,2)
yrange <- range(y)+c(-1,1)
leg.x <- xrange[2]+(xrange[1]-xrange[2])/2

fish.sub <- filter(fish, fish > 0)

plot(map, xlim = xrange, ylim = yrange,col='gray90',border=0,asp=1.5,axes=F,cex.lab=0.5)
points(x=x,y=y,pch=1,col=1,bg=0,cex=1.1)
points(x=x,y=y,pch=16,col=0,cex=0.95)
points(x=fish.sub$longitude,y=fish.sub$latitude,pch=16,col=colz,cex=0.95)
#legend(x=leg.x,y=yrange[2],legend=make.italic(name),bty='n',text.col=col)

colfunc <- colorRampPalette(RColorBrewer::brewer.pal(11,'RdYlBu'))
mapcols <- colfunc(100)[100:1]

zvec <- scales::rescale(fish$HI, to=c(0,100))
zvec <- round(zvec,0)
cexvec <- fish$cex
  
plot(map, xlim = xrange, ylim = yrange,col=0,border=1,asp=1.5,axes=F,cex.lab=0.5)
points(x=x,y=y,pch=1,col=1,cex=cexvec)
points(x=x,y=y,pch=16,col=alpha(mapcols[zvec],0.8),cex=cexvec)
polygon(x=c(-40,-40,-99,-99),y=c(75,62,62,75),col=0,border=0)

leg.x.rg <- range(xrange)[1]-(range(xrange)[1]-range(xrange)[2])*0.75
leg.x.rg <- seq(leg.x.rg,range(xrange)[2],length.out = 100)
leg.x.rg <- leg.x.rg - 5
ypos <- range(yrange)[2] - 3
points(rep(ypos,100)~leg.x.rg, col=mapcols, pch=16)
text(x=leg.x.rg[c(1,100)],y=ypos-0.2,cex=1,label=c('low','high'),pos=1)
text(x=leg.x.rg[50],y=ypos,cex=1,label='land use intensity',pos=3)

xpoints <- leg.x.rg[c(30,50,70)] - 23
points(rep(ypos,3)~xpoints, pch=1, col = 1, cex = c(0.6,1,1.4))
text(x=xpoints[2],y=ypos,cex=1,label='lake size',pos=3)
