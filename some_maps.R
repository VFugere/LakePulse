rm(list=ls())

#libraries
library(tidyverse)
library(RColorBrewer)
library(rworldmap)
library(readxl)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
devtools::source_url("https://raw.githubusercontent.com/VFugere/LakePulse/master/custom_plots.R")

#cols
cols <- brewer.pal(3, 'Dark2')
cols2 <- brewer.pal(8, 'Dark2')[c(4,5,6,8)] 

#data
load('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/2017data.RData')
zoo <- zoo.biomass.grouped 
basic.data <- filter(basic.data, year == 2017)
basic.data$ecozone <- as.factor(basic.data$ecozone)
levels(basic.data$ecozone) <- c('AH','AM','BS','MP')

#map
map <- getMap(resolution = "low")

#LP sites with 3 overlapping datasets
merged <- inner_join(basic.data,zoo)
merged <- inner_join(basic.data, bacterio) %>% 
  inner_join(phyto) %>% inner_join(zoo)
name <- 'all (147/217 sites)'
col <- 'goldenrod4'

x <- basic.data$longitude
y <- basic.data$latitude
xrange <- range(x)+c(-2,2)
yrange <- range(y)+c(-1,1)
plot(map, xlim = xrange, ylim = yrange,col='light gray',border=0,asp=1.2,axes=F,cex.lab=0.5)
points(x=x,y=y,pch=1,col=1,cex=1.2)
points(x=merged$longitude,y=merged$latitude,pch=16,col=alpha(col,0.5),cex=1.1)
legend('topright',legend=make.italic(name),bty='n',text.col=1)

x <- basic.data$longitude
y <- basic.data$latitude
xrange <- range(x)+c(-2,2)
yrange <- range(y)+c(-1,1)
plot(map, xlim = xrange, ylim = yrange,col='light gray',border=0,asp=1.2,axes=F,cex.lab=0.5)
points(x=x,y=y,pch=1,col=1,cex=1.2)
points(x=basic.data$longitude,y=basic.data$latitude,pch=16,col=alpha(cols2[basic.data$ecozone],0.5),cex=1.1)
legend('topright',legend=make.italic('2017 Lake Pulse campaign'),bty='n',text.col=1)
legend('bottomright',bty='n',legend=levels(basic.data$ecozone),pch=15,col=cols2)

pdf('~/Desktop/lakes.pdf',width=8,height = 4,pointsize = 12)
par(mfrow=c(1,2))
plot(depth_m~area,basic.data,log='x',xlab='lake area',ylab='max depth',pch=16,col=alpha(cols2[basic.data$ecozone],0.5),bty ='l')
plot(HI~area,basic.data,log='x',xlab='lake area',ylab='human impact index',pch=16,col=alpha(cols2[basic.data$ecozone],0.5),bty ='l')
dev.off()

## fish

fish <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish/Lake_ID_fish.xlsx')

merged <- inner_join(basic.data,fish, by = 'Lake_ID')
name <- 'fish occurence (56/217 sites)'
col <- 'coral3'

x <- basic.data$longitude
y <- basic.data$latitude
xrange <- range(x)+c(-2,2)
yrange <- range(y)+c(-1,1)
plot(map, xlim = xrange, ylim = yrange,col='light gray',border=0,asp=1.2,axes=F,cex.lab=0.5)
points(x=x,y=y,pch=1,col=1,cex=1.2)
points(x=merged$longitude,y=merged$latitude,pch=16,col=alpha(col,0.5),cex=1.1)
legend('topright',legend=make.italic(name),bty='n',text.col=1)
