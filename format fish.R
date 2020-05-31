### LOAD FISH DATASETS (PRSENCE/ABSENCE) FOR LAKE PULSE SITES
### FORMAT INTO A COMMUNITY MATRIX AFTER CLEANING UP TAXONOMIC DIFFERENCES AMONG PROVINCES

# code by Vincent Fugère (2020)

rm(list=ls())

library(tidyverse)
library(readxl)
library(writexl)
library(vegan)

comma2dot <- function(x){as.numeric(str_replace(x,',','.'))}
return1stval <- function(x){x[1]}
remove_empty_cols <- function(x) {return(!(all(is.na(x)) | all(x == "")))}
'%!in%' <- function(x,y)!('%in%'(x,y))

# The script below formats Qc & On fish matrices. This outputs the file 'clean_community_matrix' that I sent to Annick, and 
# which she incorporated in her excel book with all provincial sheets

# ## quebec
# 
# codes <- read_xlsx('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/taxonomic_codes.xlsx', sheet='QC')
# Qc <- read_xlsx('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Quebec_MFFP.xlsx', sheet=1)
# Qc <- select(Qc, Lake_ID, AMNE:SESP) %>% select_if(remove_empty_cols)
# Qc <- Qc %>% replace(is.na(.), 0)
# Qc <- as.data.frame(Qc)
# #nb.species <- ncol(Qc)-1
# Qc[,2:47][Qc[,2:47]>=1] <- 1 #converts to p/a
# 
# #remove bad names, change codes to binomial
# species <- colnames(Qc)
# species <- species[-1]
# species.idx <- codes$binomial[match(species,codes$code)] %>% is.na
# good.species <- species[!species.idx]
# bad.species <- species[species.idx]
# Qc <- select(Qc, -bad.species)
# colnames(Qc)[2:44] <- codes$binomial[match(good.species,codes$code)]
# 
# ## Ontario
# 
# codes <- read_xlsx('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/taxonomic_codes.xlsx', sheet='ON')
# On <- read_xlsx('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Ontario_BsM.xlsx', sheet=1)
# On <- select(On, LP.ID, `longnose gar`:`spoonhead sculpin`) %>% select_if(remove_empty_cols)
# On <- On %>% replace(is.na(.), 0)
# colnames(On)[1] <- 'Lake_ID'
# colnames(On)[44] <- 'mottled sculpin'
# On <- as.data.frame(On)
# 
# #remove bad names, change codes to binomial
# species <- colnames(On)
# species <- species[-1]
# species.idx <- codes$binomial[match(species,codes$english)] %>% is.na
# good.species <- species[!species.idx]
# bad.species <- species[species.idx]
# On <- select(On, -bad.species)
# colnames(On)[2:44] <- codes$binomial[match(good.species,codes$english)]
# 
# # verifying that names match
# 
# allfish <- bind_rows(Qc,On)
# colnames(allfish)[2:63] %>% sort #looks good except for cyprinidae.sp
# allfish <- select(allfish, -Cyprinidae_sp.)
# allfish <- allfish %>% replace(is.na(.), 0)
# 
# com <- allfish[,2:62]
# com <- com[,order(colnames(com))]
# 
# allfish <- cbind(allfish$Lake_ID,com)
# 
# write_xlsx(allfish, '/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/clean_community_matrix.xlsx')

#load other basic data

load(file='~/Google Drive/Recherche/Lake Pulse Postdoc/R/FisHab_Git/formatted_open_data/species_codes.RData')
load(file='~/Google Drive/Recherche/Lake Pulse Postdoc/R/FisHab_Git/formatted_open_data/fishbase.RData')

#load spreadsheets compiled by Annick, where she cleaned up the taxonomy for all provinces.

AL <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Annick_workbook.xlsx', sheet='AL') %>%
  select(-data_source) %>% gather(fish_species_ID,presence, FS351:FS264) %>% arrange(id_lakepulse)

BC1 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Annick_workbook.xlsx', sheet='BC-FishObs') %>%
  select(-data_source) %>% mutate_at(vars(-id_lakepulse), as.numeric) %>% gather(fish_species_ID,presence, FS005:FS360) %>% arrange(id_lakepulse)

BC2 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Annick_workbook.xlsx', sheet='BC-FishPres') %>%
  select(-data_source) %>% mutate_at(vars(-id_lakepulse), as.numeric) %>% gather(fish_species_ID,presence, FS176:FS061) %>% arrange(id_lakepulse)

NB <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Annick_workbook.xlsx', sheet='NB') %>%
  select(-data_source, -Smelt) %>% gather(fish_species_ID,presence, FS009:FS264) %>% arrange(id_lakepulse)

NS1 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Annick_workbook.xlsx', sheet='NS-FishCap') %>%
  select(-data_source) %>% mutate_at(vars(-id_lakepulse), as.numeric) %>% gather(fish_species_ID,presence, FS009:FS264) %>% arrange(id_lakepulse)

NS2 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Annick_workbook.xlsx', sheet='NS-Stock') %>%
  select(-data_source) %>% mutate_at(vars(-id_lakepulse), as.numeric) %>% gather(fish_species_ID,presence, FS327:FS252) %>% arrange(id_lakepulse)

QC1 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Annick_workbook.xlsx', sheet='QC-PE') %>%
  select(-data_source) %>% gather(fish_species_ID,presence, FS014:FS363) %>% arrange(id_lakepulse)

ON1 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Annick_workbook.xlsx', sheet='ON-MNRF') %>%
  select(-data_source) %>% mutate_at(vars(-id_lakepulse), as.numeric) %>% gather(fish_species_ID,presence, FS010:FS342) %>% arrange(id_lakepulse)

ON2 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Annick_workbook.xlsx', sheet='ON-MPO') %>%
  select(-data_source, -Centrarchidae) %>% mutate_at(vars(-id_lakepulse), as.numeric) %>% gather(fish_species_ID,presence, FS014:FS289) %>% arrange(id_lakepulse)

PEI <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Annick_workbook.xlsx', sheet='PEI') %>%
  select(-data_source) %>% mutate_at(vars(-id_lakepulse), as.numeric) %>% gather(fish_species_ID,presence, FS121:FS198) %>% arrange(id_lakepulse)

SK1 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Annick_workbook.xlsx', sheet='SK') %>%
  select(-data_source) %>% mutate_at(vars(-id_lakepulse), as.numeric) %>% gather(fish_species_ID,presence, FS264:FS327) %>% arrange(id_lakepulse)

SK2 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Annick_workbook.xlsx', sheet='HABISask_Fish') %>%
  select(-data_source) %>% mutate_at(vars(-id_lakepulse), as.numeric) %>% gather(fish_species_ID,presence, FS146:FS264) %>% arrange(id_lakepulse)

SK3 <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Annick_workbook.xlsx', sheet='HABISask_Stock') %>%
  select(-data_source) %>% mutate_at(vars(-id_lakepulse), as.numeric) %>% gather(fish_species_ID,presence, FS024:FS264) %>% arrange(id_lakepulse)

dat <- bind_rows(AL, BC1, BC2, NB, NS1, NS2, QC1, ON1, ON2, PEI, SK1, SK2, SK3)

dat <- dat %>% group_by(id_lakepulse,fish_species_ID) %>% summarize(pres = max(presence)) %>% ungroup
sum(dat$fish_species_ID %!in% species_codes$fish_species_ID) #all codes in there. phew!

LP.fish.long <- dat
LP.fish.wide <- dat %>% spread(fish_species_ID,pres)
LP.fish.wide <- LP.fish.wide %>% replace(is.na(.), 0)

writexl::write_xlsx(LP.fish.wide, '~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/LP_fish_communitymatrix.xlsx')
