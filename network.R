rm(list=ls())

#libraries
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(scales)
library(netassoc)
library(party)
library(igraph)
library(cooccur)
library(bipartite)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
devtools::source_url("https://raw.githubusercontent.com/VFugere/LakePulse/master/custom_plots.R")

#cols
cols <- brewer.pal(3, 'Dark2')
cols2 <- brewer.pal(8, 'Dark2')[c(4,5,6,8)] 

#data
load('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/2017data.RData')
basic.data <- filter(basic.data, year == 2017)
basic.data$ecozone <- as.factor(basic.data$ecozone)
levels(basic.data$ecozone) <- c('AH','AM','BS','MP')

#land use data
lulc <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/environmental/LakePulse2017_LULC_QC.csv', sep = ';', stringsAsFactors = F)
lulc <- select(lulc, -comment,-flag)
basic.data <- left_join(basic.data,lulc, by = c('Lake_ID' = 'lakepulse_id'))
rm(lulc)

## formatting and adding some taxonomic info

phyto.long <- gather(phyto, taxon, biov ,-Lake_ID)
phyto.long$class <- phytoT$KINDGOM[match(phyto.long$taxon, phytoT$totalbinomial)]
phytoT$GENUS[is.na(phytoT$GENUS)] <- phytoT$GENUS.UNKNOWN[is.na(phytoT$GENUS)]
phyto.long$genus <- phytoT$GENUS[match(phyto.long$taxon, phytoT$totalbinomial)]

phyto.genus <- phyto.long %>% filter(class != 'CYANOBACTERIA') %>% 
  group_by(Lake_ID, genus) %>%
  summarize(biov = sum(biov)) %>%
  ungroup %>%
  spread(genus, biov)

zoo <- zoo.biomass.grouped 
zooT <- readxl::read_xlsx('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/traits/zoo_trait_summary.xlsx')

zoo.long <- gather(zoo.biomass.grouped, taxon, biom ,-Lake_ID)
zoo.long$genus <- sub('\\s.*', '', zoo.long$taxon)
zoo.long$class <- zooT$class[match(zoo.long$taxon,zooT$taxon)]
zoo.long$TG <- zooT$TG[match(zoo.long$taxon,zooT$taxon)]
zoo.long$TGS <- zooT$TGS[match(zoo.long$taxon,zooT$taxon)]

zoo.genus <- zoo.long %>% filter(TG == 'herb') %>%
  group_by(Lake_ID, genus) %>%
  summarize(biom = sum(biom)) %>%
  ungroup %>%
  spread(genus, biom)

zoo.genus <- select(zoo.genus, -calanoid, -cyclopoid)

data <- inner_join(basic.data,phyto.genus) %>%
  inner_join(zoo.genus)

#excluding few large lakes
data <- filter(data, depth_m < 39)

com <- data %>% select(Acanthoceras:Simocephalus) %>% as.matrix
row.names(com) <- as.character(data$Lake_ID)
com[com>0] <- 1

## splitting lakes

#table(Hmisc::cut2(data$HI,g=3))

com.low <- com[which(data$HI < 0.0836),]
com.low  <- com.low[,which(colSums(com.low) != 0)]

com.mod <- com[which(data$HI >=0.0836 & data$HI < 0.2568),]
com.mod <- com.mod[,which(colSums(com.mod) != 0)]

com.high <- com[which(data$HI >= 0.2568),]
com.high <- com.high[,which(colSums(com.high) != 0)]

##all possible bipartite interactions

colnames(phyto.genus)[2:ncol(phyto.genus)] -> p.genera
colnames(zoo.genus)[2:ncol(zoo.genus)] -> z.genera
temp <- expand.grid(p.genera,z.genera)
all.int <- paste(temp$Var1,temp$Var2,sep='_')

#test <- make_netassoc_network(t(com.nat))

## association matrrices formatting

temp <- cooccur(t(com.low), spp_names = T, thresh=F)
cm <- temp$results
cm <- mutate_if(cm, is.factor, as.character)
cm$pair <- paste(cm$sp1_name,cm$sp2_name,sep='_')
cm.ah <- filter(cm, pair %in% all.int) #algae-herbivore inter
inter1 <- data.frame('higher' = cm.ah$sp2_name, 'lower' = cm.ah$sp1_name, 'webID' = 'low impact', freq = cm.ah$obs_cooccur, stringsAsFactors = F)

temp <- cooccur(t(com.mod), spp_names = T, thresh=F)
cm <- temp$results
cm <- mutate_if(cm, is.factor, as.character)
cm$pair <- paste(cm$sp1_name,cm$sp2_name,sep='_')
cm.ah <- filter(cm, pair %in% all.int) #algae-herbivore inter
inter2 <- data.frame('higher' = cm.ah$sp2_name, 'lower' = cm.ah$sp1_name, 'webID' = 'moderate impact', freq = cm.ah$obs_cooccur, stringsAsFactors = F)

temp <- cooccur(t(com.high), spp_names = T, thresh=F)
cm <- temp$results
cm <- mutate_if(cm, is.factor, as.character)
cm$pair <- paste(cm$sp1_name,cm$sp2_name,sep='_')
cm.ah <- filter(cm, pair %in% all.int) #algae-herbivore inter
inter3 <- data.frame('higher' = cm.ah$sp2_name, 'lower' = cm.ah$sp1_name, 'webID' = 'high impact', freq = cm.ah$obs_cooccur, stringsAsFactors = F)

inter <- bind_rows(inter1,inter2,inter3)
web <- frame2webs(inter)

#### plotting ####

cols <- c('chocolate4', 'dark green', 'darkgoldenrod3')
pdf('~/Desktop/bipartite.pdf',width = 13,height=18,pointsize = 12)
par(mfrow=c(3,1),cex=1)

# low impact web
phyto.web <- data.frame('genus' = row.names(web[[2]]),stringsAsFactors = F)
phyto.web$class <- phytoT$KINDGOM[match(phyto.web$genus, phytoT$GENUS)]
phyto.web$class.s <- 'other'
phyto.web$class.s[phyto.web$class == 'DIATOMS'] <- 'brown'
phyto.web$class.s[phyto.web$class == 'CHLOROPHYCEAE'] <- 'green'
phyto.web <- arrange(phyto.web, class.s)
nb.browns <- sum(phyto.web$class.s == 'brown')
nb.greens <- sum(phyto.web$class.s == 'green')
nb.others <- sum(phyto.web$class.s == 'other')
seq.web <- list('seq.lower' = phyto.web$genus, 'seq.higher' = colnames(web[[2]]))
web[[2]] <- sortweb(web[[2]], sort.order="seq", sequence=seq.web)
colvec <- c(rep(cols[1],nb.browns),rep(cols[2],nb.greens),rep(cols[3],nb.others))
plotweb(web[[2]], y.width.low=0.03, y.width.high=0.001, ybig=1.5,text.rot = 90, method='normal',col.low=colvec,bor.col.low = colvec,bor.col.interaction = alpha(1,0),col.interaction=rep(colvec,each=ncol(web[[2]])))

# moderate impact web
phyto.web <- data.frame('genus' = row.names(web[[3]]),stringsAsFactors = F)
phyto.web$class <- phytoT$KINDGOM[match(phyto.web$genus, phytoT$GENUS)]
phyto.web$class.s <- 'other'
phyto.web$class.s[phyto.web$class == 'DIATOMS'] <- 'brown'
phyto.web$class.s[phyto.web$class == 'CHLOROPHYCEAE'] <- 'green'
phyto.web <- arrange(phyto.web, class.s)
nb.browns <- sum(phyto.web$class.s == 'brown')
nb.greens <- sum(phyto.web$class.s == 'green')
nb.others <- sum(phyto.web$class.s == 'other')
seq.web <- list('seq.lower' = phyto.web$genus, 'seq.higher' = colnames(web[[3]]))
web[[3]] <- sortweb(web[[3]], sort.order="seq", sequence=seq.web)
colvec <- c(rep(cols[1],nb.browns),rep(cols[2],nb.greens),rep(cols[3],nb.others))
plotweb(web[[3]], y.width.low=0.03, y.width.high=0.001, ybig=1.5,text.rot = 90, method='normal',col.low=colvec,bor.col.low = colvec,bor.col.interaction = alpha(1,0),col.interaction=rep(colvec,each=ncol(web[[3]])))

# high impact web
phyto.web <- data.frame('genus' = row.names(web[[1]]),stringsAsFactors = F)
phyto.web$class <- phytoT$KINDGOM[match(phyto.web$genus, phytoT$GENUS)]
phyto.web$class.s <- 'other'
phyto.web$class.s[phyto.web$class == 'DIATOMS'] <- 'brown'
phyto.web$class.s[phyto.web$class == 'CHLOROPHYCEAE'] <- 'green'
phyto.web <- arrange(phyto.web, class.s)
nb.browns <- sum(phyto.web$class.s == 'brown')
nb.greens <- sum(phyto.web$class.s == 'green')
nb.others <- sum(phyto.web$class.s == 'other')
seq.web <- list('seq.lower' = phyto.web$genus, 'seq.higher' = colnames(web[[1]]))
web[[1]] <- sortweb(web[[1]], sort.order="seq", sequence=seq.web)
colvec <- c(rep(cols[1],nb.browns),rep(cols[2],nb.greens),rep(cols[3],nb.others))
plotweb(web[[1]], y.width.low=0.03, y.width.high=0.001, ybig=1.5,text.rot = 90,method='normal',col.low=colvec,bor.col.low = colvec,bor.col.interaction = alpha(1,0),col.interaction=rep(colvec,each=ncol(web[[1]])))

dev.off()

#### extracting a few metrics from the web ####

networklevel(web[[2]], 'connectance')
networklevel(web[[2]], 'weighted connectance')
networklevel(web[[2]], 'web asymmetry')
networklevel(web[[2]], 'links per species')
networklevel(web[[2]], 'number of compartments')
networklevel(web[[2]], 'nestedness')
networklevel(web[[2]], 'linkage density')
networklevel(web[[2]], 'Fisher alpha')
networklevel(web[[2]], 'interaction evenness')

networklevel(web[[2]], 'binary')
networklevel(web[[2]], 'topology')

networklevel(web[[2]], 'info')
networklevel(web[[3]], 'info')
networklevel(web[[1]], 'info')

networklevel(web[[2]], 'binary')
networklevel(web[[3]], 'binary')
networklevel(web[[1]], 'binary')
