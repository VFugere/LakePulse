# Analysis script to model max trophic level of fish communities in Lake Pulse sites,
# using either only abiotic parameters or using also zooplankton. Does including zoopplankton
# increase predictive power?

rm(list=ls())

#libraries
library(tidyverse)
library(readxl)
library(vegan)
library(devtools)
library(party)
library(rworldmap)

#functions
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
source_url("https://raw.githubusercontent.com/VFugere/LakePulse/master/custom_plots.R")
mean.nona <- function(x){mean(x, na.rm=T)}
sum.nona <- function(x){sum(x, na.rm=T)}

### load data ###

map <- getMap(resolution='low')

## lake pulse abiotic data

load('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/allenvdata.RData')

#need to update this and select relevant variables
env <- left_join(basic.data, chla) %>% left_join(rbr) %>% left_join(secchi) %>% left_join(strati) %>% left_join(TN) %>% left_join(climate)
land.use <- left_join(basic.data, lulc.watershed) %>% left_join(lulc.buffer) %>% select(Lake_ID, HII:prop.urban.buffer)

## fish

load(file='~/Google Drive/Recherche/Lake Pulse Postdoc/R/FisHab_Git/formatted_open_data/species_codes.RData')
load(file='~/Google Drive/Recherche/Lake Pulse Postdoc/R/FisHab_Git/formatted_open_data/fishbase.RData')
fish <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/LP_fish_communitymatrix.xlsx')
fish.long <- fish %>% pivot_longer(FS005:FS364, names_to = 'species', values_to = 'presence')

#### calculating max TL ####

fish.long$species.name <- species_codes$clean.species.name[match(fish.long$species,species_codes$fish_species_ID)]
#only keeping stuff id'd to species
fish.long <- drop_na(fish.long)
#merging species codes that are the same species by grouping by name instead of FS code
fish.long <- fish.long %>% group_by(id_lakepulse,species.name) %>% summarize(pres = max(presence)) %>% ungroup
#to wide
fish.wide <- fish.long %>% spread(species.name,pres)
fish.wide <- fish.wide %>% replace(is.na(.), 0)
#adding trophic level
for(co in 2:ncol(fish.wide)){
  sp <- colnames(fish.wide)[co]
  sp <- str_replace(sp, '_', ' ')
  tl <- estimate.df %>% filter(Species == sp) %>% pull(Troph)
  fish.wide[,co] <- fish.wide[,co]*tl
}

fish.wide$max.trophic.level <- apply(fish.wide[,2:ncol(fish.wide)], 1, max)

# send to Cindy
# cindy <- select(fish.wide, id_lakepulse, max.trophic.level) %>% filter(max.trophic.level > 0)
# writexl::write_xlsx(cindy, '~/Desktop/Cindy.xlsx')

#### adding zooplankton metrics ####

d2019 <- read.csv('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/ALLfinal_grouping2019.csv', stringsAsFactors = F)
colnames(d2019) <- str_replace(colnames(d2019), '\\.\\.', ' ')
colnames(d2019) <- str_replace(colnames(d2019), 'spp\\.', 'sp\\.')
colnames(d2019) <- str_replace(colnames(d2019), 'sp\\.', 'sp')
colnames(d2019) <- str_replace(colnames(d2019), '\\.', ' ')
d2017 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/ALLfinal_grouping2017.csv', stringsAsFactors = F) %>%
  mutate_at(vars(Ergasilus.spp.:Simocephalus..spp.), as.numeric)
colnames(d2017) <- str_replace(colnames(d2017), '\\.\\.', ' ')
colnames(d2017) <- str_replace(colnames(d2017), 'spp\\.', 'sp\\.')
colnames(d2017) <- str_replace(colnames(d2017), 'sp\\.', 'sp')
colnames(d2017) <- str_replace(colnames(d2017), '\\.', ' ')
d2018 <- read.csv2('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/ALLfinal_grouping2018.csv', stringsAsFactors = F) %>%
  mutate_at(vars(Ergasilus.spp.:Simocephalus..spp.), as.numeric)
colnames(d2018) <- str_replace(colnames(d2018), '\\.\\.', ' ')
colnames(d2018) <- str_replace(colnames(d2018), 'spp\\.', 'sp\\.')
colnames(d2018) <- str_replace(colnames(d2018), 'sp\\.', 'sp')
colnames(d2018) <- str_replace(colnames(d2018), '\\.', ' ')
zoo <- bind_rows(d2017,d2018,d2019)
zoo[is.na(zoo)] <- 0
zoo <- zoo[,c(1,order(names(zoo[,2:ncol(zoo)]))+1)]
rm(d2017,d2018,d2019)

species_list <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/zoo_species_list.xlsx')
cladocerans <- filter(species_list, division_1 == 'Cladocera') %>% pull(species)
daphniids <- filter(species_list, division_2 == 'Daphniidae') %>% pull(species)
chydorids <- filter(species_list, division_2 == 'Chydoridae') %>% pull(species)
sidids <- filter(species_list, division_2 == 'Sididae') %>% pull(species)
copepods <- filter(species_list, division_1 == 'Copepoda') %>% pull(species)
cyclopoids <- filter(species_list, division_2 == 'Cyclopoid') %>% pull(species)
calanoids <- filter(species_list, division_2 == 'Calanoid') %>% pull(species)

zoobiomass <- data.frame('Lake_ID' = zoo$Lake_ID,'biomass' = rowSums(zoo[,2:108]))

#saving RUE for later
zoochla <- inner_join(zoobiomass,chla) %>% mutate(zoo.rue = biomass/chla) %>% select(Lake_ID, zoo.rue)

zoocompo <- data.frame('cladocerans' = rowSums(zoo[,cladocerans]),
                       'daphniids' = rowSums(zoo[,daphniids]),
                       'chydorids' = rowSums(zoo[,chydorids]),
                       'sidids' = rowSums(zoo[,sidids]),
                       'copepods' = rowSums(zoo[,copepods]),
                       'cyclopoids' = rowSums(zoo[,cyclopoids]),
                       'calanoids' = rowSums(zoo[,calanoids]))

zoorelcompo <- zoocompo / zoobiomass$biomass
colnames(zoorelcompo) <- paste0('rel_',colnames(zoorelcompo))

#adding diversity
div <- data.frame('Lake_ID' = zoo$Lake_ID)
div$richness <- specnumber(zoo[,2:108])
div$alpha <- exp(diversity(zoo[,2:108]))

zoobiomass <- bind_cols(zoobiomass, zoocompo, zoorelcompo) %>% left_join(div)
rm(zoocompo,zoorelcompo,div)

# com <- zoo[,2:108]
# com <- decostand(com, 'hellinger')
# ordi <- rda(com ~ 1)
# g <- scores(ordi)$sites[,1:2]
# g <- as.data.frame(g)
# 
# zoobiomass <- bind_cols(zoobiomass,g)

d2017 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/all raw data 2017.xlsx', sheet='raw') %>% rename(`biomass factor` = `biomass factor (µg/ind)`) %>% select(ID_lakepulse,genus,species,division,`#individuals counted`,`# / L`,`biomass factor`,D1:D10)
d2018 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/all raw data 2018.xlsx', sheet='raw') %>% select(ID_lakepulse,genus,species,division,`#individuals counted`,`# / L`,`biomass factor`,D1:D10)
d2019 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/all raw data 2019.xlsx', sheet='raw') %>% select(ID_lakepulse,genus,species,division,`#individuals counted`,`# / L`,`biomass factor`,D1:D10)

zoo <- bind_rows(d2017,d2018,d2019) %>% rename(Lake_ID = ID_lakepulse, abundance = `#individuals counted`, density = `# / L`, mean.mass = `biomass factor`)

#chaoborus presence
filter(zoo, genus == 'Chaoborus') %>% select(Lake_ID, density)
alllakes <- select(basic.data, Lake_ID)
chaoborus <- filter(zoo, genus == 'Chaoborus') %>% select(Lake_ID, density) %>% right_join(alllakes) %>% as.data.frame
chaoborus[is.na(chaoborus)] <- 0
chaoborus <- chaoborus %>% group_by(Lake_ID) %>% summarize(chaoborus = sum(density))

#remove
zoo <- filter(zoo, division %in% c('Cladocera','Copepoda'))
rm(d2017,d2018,d2019)

#calculate mean body size, excluding 0s and NAs
sizes <- select(zoo, D1:D10)
sizes[sizes == 0] <- NA
sizes$mean <- apply(sizes, MARGIN = 1, FUN = mean.nona)
zoo$mean.size <- sizes$mean*1000 #in um
rm(sizes)
#zoo <- select(zoo, -abundance, -(D1:D10))
zoo$taxon.biomass <- zoo$density*zoo$mean.mass

## community weighted mean body size

zoo$tot.size <- zoo$density*zoo$mean.size
cwm <- zoo %>% group_by(Lake_ID) %>% summarize(total.zoo.length = sum(tot.size), total.density = sum(density))
cwm$cwms <- cwm$total.zoo.length/cwm$total.density
cwm <-cwm %>% select(Lake_ID,cwms)
# test <- filter(zoo, Lake_ID == '04-572') #should get 504
# weighted.mean(test$mean.size,test$density) #good
cwm.clad <- zoo %>% filter(division == 'Cladocera') %>% group_by(Lake_ID) %>% summarize(total.zoo.length = sum(tot.size), total.density = sum(density))
cwm.clad$cwms <- cwm.clad$total.zoo.length/cwm.clad$total.density
cwm.clad <-cwm.clad %>% select(Lake_ID,cwms) %>% rename(cwms.clad = cwms)
cwm.cop <- zoo %>% filter(division == 'Copepoda') %>% group_by(Lake_ID) %>% summarize(total.zoo.length = sum(tot.size), total.density = sum(density))
cwm.cop$cwms <- cwm.cop$total.zoo.length/cwm.cop$total.density
cwm.cop <-cwm.cop %>% select(Lake_ID,cwms) %>% rename(cwms.cop = cwms)
cwm <- left_join(cwm,cwm.clad) %>% left_join(cwm.cop)

## chla to TN

chlatn <- left_join(chla, TN) %>% mutate(phyto.rue = chla/TN) %>% dplyr::select(Lake_ID, phyto.rue)

zoobiomass <- left_join(zoobiomass,cwm) %>% left_join(chaoborus) %>% left_join(zoochla) %>% left_join(chlatn)

zoo <- zoobiomass %>% arrange(Lake_ID)
rm(zoobiomass,cwm,cwm.clad,cwm.cop,chaoborus,alllakes,chlatn,zoochla)

# # does body size project look promising?
# 
# test <- left_join(select(env, Lake_ID, latitude, altitude, mean_temp, rbr_temp), select(zoo, Lake_ID, cwms:cwms.cop)) %>% left_join(land.use)
# corrgram::corrgram(test)
# library(mgcv)
# library(itsadug)
# m <- gam(cwms ~ ti(latitude) + ti(altitude) + ti(latitude, altitude), data=subset(test, !is.na(cwms)))
# summary(m)
# plot(m)

#### machine learning ####

library(randomForest)
library(party)
library(caret)
library(gbm)
library(skimr)

TL <- fish.wide %>% select(id_lakepulse, max.trophic.level) %>% rename(Lake_ID = id_lakepulse, fishTL = max.trophic.level) %>% arrange(Lake_ID)
TL <- filter(TL, Lake_ID %in% basic.data$Lake_ID)

skim(filter(env, Lake_ID %in% TL$Lake_ID))
subenv <- env %>% select(-feow, -cont.watershed, -(DO_sat:pH), -(SPC:hypolimnion_depth))

skim(filter(land.use, Lake_ID %in% TL$Lake_ID))
skim(filter(zoo, Lake_ID %in% TL$Lake_ID))

dat <- filter(TL, fishTL > 3) %>% left_join(subenv) %>% left_join(land.use) %>% left_join(zoo) %>% drop_na %>%
  mutate_if(is.character, as.factor) %>% select(-Lake_ID)

dat.e <- select(dat, fishTL:total_precip)
dat.l <- select(dat, fishTL, HII:prop.urban.buffer)
dat.b <- select(dat, fishTL, biomass:phyto.rue)
dat.el <- select(dat, fishTL:prop.urban.buffer)
dat.lb <- select(dat, fishTL, HII:phyto.rue)
dat.eb <- select(dat, fishTL:total_precip, biomass:phyto.rue)
dat.elb <- dat

results <- data.frame('Env' = numeric(0), 'LU' = numeric(0), 'Bio' = numeric(0),
                      'Env_LU' = numeric(0), 'LU_Bio' = numeric(0), 'Env_Bio' = numeric(0),
                      'All_pred' = numeric(0))

for(i in 1:100){
  
  smpl <- sample(nrow(dat), round(nrow(dat)*.75))
  
  #e
  training.set <- dat.e[smpl,]
  test.set <- dat.e[-smpl,]
  mod <-cforest(fishTL ~ ., training.set, controls = cforest_unbiased(mtry=6))
  test.results <- predict(mod, newdata = test.set)
  results[i,1] <- postResample(pred = test.results, obs = test.set$fishTL)[2]

  #lu
  training.set <- dat.l[smpl,]
  test.set <- dat.l[-smpl,]
  mod <-cforest(fishTL ~ ., training.set, controls = cforest_unbiased(mtry=6))
  test.results <- predict(mod, newdata = test.set)
  results[i,2] <- postResample(pred = test.results, obs = test.set$fishTL)[2]
  
  #bio
  training.set <- dat.b[smpl,]
  test.set <- dat.b[-smpl,]
  mod <-cforest(fishTL ~ ., training.set, controls = cforest_unbiased(mtry=6))
  test.results <- predict(mod, newdata = test.set)
  results[i,3] <- postResample(pred = test.results, obs = test.set$fishTL)[2]
  
  #el
  training.set <- dat.el[smpl,]
  test.set <- dat.el[-smpl,]
  mod <-cforest(fishTL ~ ., training.set, controls = cforest_unbiased(mtry=6))
  test.results <- predict(mod, newdata = test.set)
  results[i,4] <- postResample(pred = test.results, obs = test.set$fishTL)[2]
  
  #lb
  training.set <- dat.lb[smpl,]
  test.set <- dat.lb[-smpl,]
  mod <-cforest(fishTL ~ ., training.set, controls = cforest_unbiased(mtry=6))
  test.results <- predict(mod, newdata = test.set)
  results[i,5] <- postResample(pred = test.results, obs = test.set$fishTL)[2]
  
  #eb
  training.set <- dat.eb[smpl,]
  test.set <- dat.eb[-smpl,]
  mod <-cforest(fishTL ~ ., training.set, controls = cforest_unbiased(mtry=6))
  test.results <- predict(mod, newdata = test.set)
  results[i,6] <- postResample(pred = test.results, obs = test.set$fishTL)[2]
  
  #elb
  training.set <- dat.elb[smpl,]
  test.set <- dat.elb[-smpl,]
  mod <-cforest(fishTL ~ ., training.set, controls = cforest_unbiased(mtry=6))
  test.results <- predict(mod, newdata = test.set)
  results[i,7] <- postResample(pred = test.results, obs = test.set$fishTL)[2]
  
}

#boxplot

colors <- wesanderson::wes_palette('FantasticFox1',3)

colnames(results) <- 1:7
plotdat <- gather(results) %>% mutate(key = as.numeric(key))
boxplot(value~key, plotdat, lty=1, staplewex=0, whisklwd=1,boxwex=0.6, outline=F,boxlwd=1, medlwd=1, col=alpha(colors[2],0.3), ylab=expression(model~R^2), xlab=NULL, xaxt = 'n')
points(x=jitter(as.numeric(plotdat$key),0.5),y=plotdat$value,pch=1,cex=0.5,col=colors[1])
axis(1,at=1:7, labels = c('E','L','B','E+L','L+B','E+B','E+L+B'),lwd=0,lwd.ticks = 1)
title(xlab='Predictors included in model')

#barplot of variable importance using last model

vars <- varimp(mod, conditional = T)
vars[vars < 0] <- 0
vars <- vars * 100
vars <- sort(vars)
names(vars) <- c('chlorophyll a','zoo biomass','cladoceran biomass','Sididae biomass','copepod biomass',
                 'cyclopoid biomass','calanoid biomass','% Daphnidae','% Sididae','% cyclopoid',
                 'mean copepod body size','zoo/phyto biomass','phyto biomass/TN','Chydoridae biomass','% Chydoridae',
                 'Daphnidae biomass','stratified','max air temp','human impact index','% agriculture near shore',
                 '% agriculture watershed','% mines near shore','mean zoo body size','% urban watershed','% urban near shore',
                 'total nitrogen','mean cladoceran body size','% natural near shore','% cladocerans','residence time',
                 '% copepod','min air temp','% pasture near shore','% forestry near shore','% mines watershed',
                 '% natural watershed','Chaoborus density','mean air temp','zoo richness','% forestry watershed',
                 'altitude','conductivity','shoreline development','% pasture watershed','water temp',
                 '% calanoids','zoo alpha diversity','annual precipations','salinity','latitude',
                 'max lake depth','province','longitude','ecozone','lake area')

pdf('~/Desktop/var_influence.pdf',width = 5,height = 8,pointsize = 8)
par(mar=c(4,12,1,1))
barplot(vars,col=colors[3],border=0,horiz=T,las=1)
title(xlab='relative influence (%)')
dev.off()
par(mar=c(4,4,1,1))

##### GAMMs #####

library(mgcv)
library(itsadug)

dat$log.area <- log10(dat$area)

par(mfrow=c(3,3),mar=c(4,4,1,1),oma=c(1,1,1,1))

for(prov in levels(dat$province)){
  subdat <- filter(dat, province == prov)
  if(nrow(subdat)>10){
    mod <- gam(fishTL ~ s(log.area, k = floor(nrow(subdat)/2)), data = subdat)
    plot_smooth(mod, view = 'log.area',col=colors[1],rug=F,print.summary=F,xlab = 'log10 lake area',ylab = 'MaxTL', hide.label = T,bty='o')
    legend('topleft',legend=prov,bty='n')
    pval <- round(summary(mod)$s.table[1,4],4)
    if(pval < 0.0001){pval <- 0.0001}
    legend('bottomright',legend=bquote(italic(p) < .(pval)),bty='n')
  }
}

par(mfrow=c(1,1),mar=c(4,4,2,2),oma=c(1,1,1,1))

cols <- viridis::viridis(50)
mod <- gam(fishTL ~ te(latitude,longitude,bs='gp') + ti(log.area, k=8) + ti(alpha, k=8) + ti(log.area,alpha,k=8) + s(province, bs='re'), data=dat, correlation = corSpher(form = ~ latitude + longitude))
summary(mod)
fvisgam(mod, view=c('log.area','alpha'),dec=1,color=cols,hide.label=T,rm.ranef = T,xlab='log10 lake area',ylab='zooplankton diversity',main = 'Max trophic level')

fvisgam(mod, view=c('log.area','alpha'),color=cols,hide.label=T,rm.ranef = T,plot.type='persp',theta=45,main=NULL,zlab = 'max trophic level',xlab='lake area',ylab='zooplankton diversity')

tree <- ctree(fishTL ~ log.area + alpha, dat)
plot(tree)

# map

theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")
library(viridis)

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world) + 
  geom_sf(color = "black", fill = "lightgreen")

xlims <- range(dat$longitude)+c(-5,5)
ylims <- range(dat$latitude)+c(-8,8)

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = xlims, ylim = ylims) +
  #annotation_scale(location = "bl", width_hint = 0.35) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  geom_point(data=dat,aes(x=longitude,y=latitude,col=fishTL)) + scale_color_viridis()

dat$resid <- resid(mod)

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = xlims, ylim = ylims) +
  #annotation_scale(location = "bl", width_hint = 0.35) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  geom_point(data=dat,aes(x=longitude,y=latitude,col=resid)) + scale_color_viridis()

gam.check(mod)
summary(mod)

simpler.mod <- gam(fishTL ~ te(latitude,longitude,bs='gp') + s(log.area, k=8) + s(province, bs='re'), data=dat, correlation = corSpher(form = ~ latitude + longitude))
summary(simpler.mod)
gam.check(simpler.mod)
plot(simpler.mod)

plot_smooth(simpler.mod, view = 'log.area',col=colors[1],rug=F,print.summary=F,xlab = 'log10 lake area',ylab = 'MaxTL', hide.label = T,bty='o')
pval <- round(summary(simpler.mod)$s.table[2,4],4)
if(pval < 0.0001){pval <- 0.0001}
legend('top',legend='p < 0.0001',bty='n')

par(mfrow=c(1,2))
plot_smooth(mod, view = 'log.area',col=colors[1],rug=F,print.summary=F,xlab = 'log10 lake area',ylab = 'MaxTL', hide.label = T,bty='o')
legend('top',legend='p < 0.0001',bty='n')
plot_smooth(mod, view = 'alpha',col=colors[3],rug=F,print.summary=F,xlab = 'zooplankton diversity',ylab = 'MaxTL', hide.label = T,bty='o')
legend('top',legend='p = 0.17',bty='n')


########
  
## simplest of random forest with randomForest

forest <- randomForest(fishTL ~ ., training.set, localImp = TRUE)
importance(forest, scale=T)
forest

## conditional random forest with party

cforest1 <-cforest(fishTL ~ ., training.set, controls = cforest_unbiased(mtry=6)) # mtry optimizied with a loop
train.results <- predict(cforest1, OOB=T)
test.results <- predict(cforest1, newdata = test.set)
postResample(pred = train.results, obs = training.set$fishTL)
postResample(pred = test.results, obs = test.set$fishTL)
plot(test.results~test.set$fishTL)
vars <- varimp(cforest1, conditional = T)
sort(vars)

## GBM

gbm1 <- gbm(fishTL ~ ., data=training.set, n.trees = 500, interaction.depth = 3)
#gbm.perf(gbm1)
plot(predict(gbm1, newdata = test.set, n.trees = 500)~test.set$fishTL)
cor(predict(gbm1, newdata = test.set, n.trees = 500),test.set$fishTL)^2

# #trying to optimize with caret
# 
# fitControl <- trainControl(## 10-fold CV
#   method = "repeatedcv",
#   number = 10,
#   ## repeated ten times
#   repeats = 10)
# 
# gbmGrid <-  expand.grid(interaction.depth = 1:9, 
#                        n.trees = 500, 
#                        shrinkage = c(0.1,0.5,0.9),
#                        n.minobsinnode = 20)
# 
# 
# gbm2 <- train(fishTL ~ ., data=training.set, 
#                  method = "gbm", 
#                  trControl = fitControl, 
#                  verbose = FALSE, 
#                  tuneGrid = gbmGrid)
# 
# gbm2
# plot(gbm2)
# 
# plot(predict(gbm2, newdata = test.set, n.trees = 500)~test.set$fishTL) #exact same thing

# ## Bayesian Additive Regression Trees (BART)
# 
# library(BayesTree)
# 
# bart1 <- bart(x.train = as.data.frame(training.set[,-1]), y.train = training.set$fishTL, x.test=as.data.frame(test.set[,-1]))
# bart1
# plot(bart1) #awful
# 
# ## BART with caret & bartMachine
# 
# library(bartMachine)
# 
#
#### refitting RF, cforest and GBM with zoo predictors ####

training.set <- dat.zoo[sample,]
test.set <- dat.zoo[-sample,]

forest2<-randomForest(fishTL ~ ., training.set, localImp = TRUE)
vars <- importance(forest2, scale=T)
test.results <- predict(forest2, newdata = test.set)
cor(test.results,test.set$fishTL)^2

cforest2 <-cforest(fishTL ~ ., training.set, controls = cforest_unbiased(mtry=6)) # mtry optimizied with a loop
train.results <- predict(cforest2, OOB=T)
test.results <- predict(cforest2, newdata = test.set)
postResample(pred = train.results, obs = training.set$fishTL)
postResample(pred = test.results, obs = test.set$fishTL)
#cor(test.results,test.set$fishTL)^2 #same indeed
plot(test.results~test.set$fishTL)
vars <- varimp(cforest2, conditional = T)
sort(vars)

gbm2 <- gbm(fishTL ~ ., data=training.set, n.trees = 500, interaction.depth = 3)
plot(predict(gbm2, newdata = test.set, n.trees = 500)~test.set$fishTL)
cor(predict(gbm2, newdata = test.set, n.trees = 500),test.set$fishTL)^2

##### traditional modelling ####

corrgram::corrgram(dat.zoo)

dat.zoo$chao.pa <- dat.zoo$chaoborus
dat.zoo$chao.pa[dat.zoo$chao.pa > 1] <- 1

## varpart

#space <- scale(dat.zoo[,c('latitude','altitude')])
morpho <- scale(dat.zoo[,c('area','depth_m','Shore_dev')])
impact <- scale(dat.zoo[,c('HII','prop.ag.buffer','TN')])
zooplank <- scale(dat.zoo[,31:34])

mod <- varpart(Y=dat.zoo$fishTL,morpho,land.use,zooplank)
mod
plot(mod)

## gam

library(mgcv)

long <- basic.data %>% filter(Lake_ID %in% sites) %>% pull(longitude)

dat.zoo$long <- long

m1 <- gam(fishTL ~ s(area,k=5) + s(depth_m,k=5) + te(latitude,long,k=4), correlation = corSpher(form = ~ latitude + long), data = dat.zoo)
summary(m1)
gam.check(m1)
plot(m1)

m2 <- gam(fishTL ~ s(area) + s(depth_m) + te(latitude,long) + s(cwms), correlation = corSpher(form = ~ latitude + long), data = dat.zoo)
summary(m2)
