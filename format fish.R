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

## quebec

codes <- read_xlsx('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/taxonomic_codes.xlsx', sheet='QC')
Qc <- read_xlsx('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Quebec_MFFP.xlsx', sheet=1)
Qc <- select(Qc, Lake_ID, AMNE:SESP) %>% select_if(remove_empty_cols)
Qc <- Qc %>% replace(is.na(.), 0)
Qc <- as.data.frame(Qc)
#nb.species <- ncol(Qc)-1
Qc[,2:47][Qc[,2:47]>=1] <- 1 #converts to p/a

#remove bad names, change codes to binomial
species <- colnames(Qc)
species <- species[-1]
species.idx <- codes$binomial[match(species,codes$code)] %>% is.na
good.species <- species[!species.idx]
bad.species <- species[species.idx]
Qc <- select(Qc, -bad.species)
colnames(Qc)[2:44] <- codes$binomial[match(good.species,codes$code)]

## Ontario

codes <- read_xlsx('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/taxonomic_codes.xlsx', sheet='ON')
On <- read_xlsx('/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/Ontario_BsM.xlsx', sheet=1)
On <- select(On, LP.ID, `longnose gar`:`spoonhead sculpin`) %>% select_if(remove_empty_cols)
On <- On %>% replace(is.na(.), 0)
colnames(On)[1] <- 'Lake_ID'
colnames(On)[44] <- 'mottled sculpin'
On <- as.data.frame(On)

#remove bad names, change codes to binomial
species <- colnames(On)
species <- species[-1]
species.idx <- codes$binomial[match(species,codes$english)] %>% is.na
good.species <- species[!species.idx]
bad.species <- species[species.idx]
On <- select(On, -bad.species)
colnames(On)[2:44] <- codes$binomial[match(good.species,codes$english)]

# verifying that names match

allfish <- bind_rows(Qc,On)
colnames(allfish)[2:63] %>% sort #looks good except for cyprinidae.sp
allfish <- select(allfish, -Cyprinidae_sp.)
allfish <- allfish %>% replace(is.na(.), 0)

com <- allfish[,2:62]
com <- com[,order(colnames(com))]

allfish <- cbind(allfish$Lake_ID,com)

write_xlsx(allfish, '/Users/vincentfugere/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/fish_output/clean_community_matrix.xlsx')
