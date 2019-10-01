rm(list=ls())

#libraries
library(tidyverse)
library(vegan)
library(RColorBrewer)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))

#cols
cols <- brewer.pal(3, 'Dark2')

#data
load('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/2017data.RData')

trt <- env.data %>% select(Lake_ID, HI_class, Area_class, Ecozone.x) %>% rename(ecozone = Ecozone.x, area = Area_class, HI = HI_class)

#### zooplankton ####

com <- bacterio

meta <- com %>% select(Lake_ID)
meta <- left_join(meta, trt, by = 'Lake_ID')
meta$HI_Ez <- paste(meta$ecozone,meta$HI,sep='_')
meta <- mutate_at(meta, vars(HI:HI_Ez), as.factor)

com <- com[,2:ncol(com)] %>% as.matrix
com[is.na(com)] <- 0

# dm <- vegdist(log1p(com),method = 'bray')
# adonis(dm~meta$HI*meta$ecozone*meta$area)
# adonis(dm~meta$HI+meta$ecozone+meta$area+meta$HI:meta$ecozone)
# #area, ecozone, HI have a signigicant impact on composition.
# #significant HI by ecozone interaction

##### including all sites together ####

par(mfrow=c(2,2))

#HI impact
lines <- which(meta$area != 'large')
com.sub <- com[lines,]
meta.sub <- meta[lines,]
# com.sub[com.sub > 0] <- 1
# com.sub <- decostand(com.sub,method='total')
# com.sub <- decostand(com.sub,method='hellinger')
# com.sub <- log1p(com.sub)
# com.sub <- sqrt(com.sub)
com.sub <- wisconsin(com.sub)
dm <- vegdist(com.sub,method = 'bray')
cc <- adonis(dm~meta.sub$HI)
p.cc <- round(cc$aov.tab$`Pr(>F)`,2)[1]
bd <- anova(betadisper(d=dm, group=meta.sub$HI, type='centroid'))
p.bd <- round(bd$`Pr(>F)`[1],2)
plot(betadisper(d=dm, group=meta.sub$HI, type='centroid'),hull=T,label=F,segments=F,col=cols,pch=rep(1,3),cex=0.2)
legend('topleft',bty='n',legend=c(p.cc,p.bd))
legend('topright',bty='n',legend=levels(meta.sub$HI),pch=16,col=cols)
distances <- (betadisper(d=dm, group=meta.sub$HI, type='centroid'))$distances
boxplot(distances~meta.sub$HI)

#Area impact
lines <- which(meta$HI != 'high')
com.sub <- com[lines,]
meta.sub <- meta[lines,]
# com.sub[com.sub > 0] <- 1
# com.sub <- decostand(com.sub,method='total')
# com.sub <- decostand(com.sub,method='hellinger')
# com.sub <- log1p(com.sub)
# com.sub <- sqrt(com.sub)
com.sub <- wisconsin(com.sub)
dm <- vegdist(com.sub,method = 'bray')
cc <- adonis(dm~meta.sub$area)
p.cc <- round(cc$aov.tab$`Pr(>F)`,2)[1]
bd <- anova(betadisper(d=dm, group=meta.sub$area, type='centroid'))
p.bd <- round(bd$`Pr(>F)`[1],2)
plot(betadisper(d=dm, group=meta.sub$area, type='centroid'),hull=T,label=F,segments=F,col=cols,pch=rep(1,3),cex=0.2)
legend('topleft',bty='n',legend=c(p.cc,p.bd))
legend('topright',bty='n',legend=levels(meta.sub$area),pch=16,col=cols)
distances <- (betadisper(d=dm, group=meta.sub$area, type='centroid'))$distances
boxplot(distances~meta.sub$area)

######

### separating ecozones

HI_classes <- levels(meta$HI)
ecozones <- levels(meta$ecozone)
areas <- levels(meta$area)

par(mfrow=c(2,2))

#effect of HI

for(ecozone in ecozones){
  
  lines <- which(meta$ecozone == ecozone & meta$area != 'large')
  com.sub <- com[lines,]
  meta.sub <- meta[lines,]
  
  # com.sub[com.sub > 0] <- 1
  # com.sub <- decostand(com.sub,method='total')
  # com.sub <- decostand(com.sub,method='hellinger')
  # com.sub <- log1p(com.sub)
  # com.sub <- sqrt(com.sub)
  com.sub <- wisconsin(com.sub)
   
  dm <- vegdist(com.sub,method = 'bray')
  
  #testing whether HI differ in centroid location and mean dissimilarity
  cc <- adonis(dm~meta.sub$HI)
  p.cc <- round(cc$aov.tab$`Pr(>F)`,2)[1]
  
  bd <- anova(betadisper(d=dm, group=meta.sub$HI, type='centroid'))
  p.bd <- round(bd$`Pr(>F)`[1],2)
  
  plot(betadisper(d=dm, group=meta.sub$HI, type='centroid'),hull=T,label=F,segments=F,col=cols,pch=rep(1,3),cex=0.2)
  legend('topleft',bty='n',legend=c(p.cc,p.bd))
  legend('topright',bty='n',legend=levels(meta.sub$HI),pch=16,col=cols)
  
}

#effect of area
for(ecozone in ecozones){
  
  lines <- which(meta$ecozone == ecozone & meta$HI != 'high')
  com.sub <- com[lines,]
  meta.sub <- meta[lines,]
  
  # com.sub[com.sub > 0] <- 1
  # com.sub <- decostand(com.sub,method='total')
  # com.sub <- decostand(com.sub,method='hellinger')
  # com.sub <- log1p(com.sub)
  # com.sub <- sqrt(com.sub)
  com.sub <- wisconsin(com.sub)
   
  dm <- vegdist(com.sub,method = 'bray')
  
  #testing whether area differ in centroid location and mean dissimilarity
  cc <- adonis(dm~meta.sub$area)
  p.cc <- round(cc$aov.tab$`Pr(>F)`,2)[1]
  
  bd <- anova(betadisper(d=dm, group=meta.sub$area, type='centroid'))
  p.bd <- round(bd$`Pr(>F)`[1],2)
  
  plot(betadisper(d=dm, group=meta.sub$area, type='centroid'),hull=T,label=F,segments=F,col=cols,pch=rep(1,3),cex=0.2)
  legend('topleft',bty='n',legend=c(p.cc,p.bd))
  legend('topright',bty='n',legend=levels(meta.sub$area),pch=16,col=cols)
  
}
