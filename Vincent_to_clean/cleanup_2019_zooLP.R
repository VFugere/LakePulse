#for 2019, do what Cindy did for other years, splitting biomass of un-Id'd species across species, and grouping copepodids with adults

dat <- read_xlsx('~/Google Drive/Recherche/Lake Pulse Postdoc/data/LP/zooplankton/all raw data 2019.xlsx', sheet = 'clean biomass')

# first group all copepodids id'd to species level with adults
taxo <- dat %>% filter(division == 'Copepoda') %>%
  distinct(`species name`) %>%
  filter(`species name` %!in% c("calanoid copepodid", "cyclopoid copepodid")) %>%
  rename(name = `species name`) %>%
  arrange(name)
taxo$out <- taxo$name
taxo$out <- str_remove(taxo$out, pattern = ' copepodid')
taxo <- as.data.frame(taxo)
taxo$out[taxo$name == 'Heterocope  septentrionalis'] <- 'Heterocope septentrionalis'

dat[dat$`species name` %in% taxo$name,'species name'] <- taxo$out[match(dat$`species name`[dat$`species name` %in% taxo$name],taxo$name)]
colnames(dat)[2:4] <- c('name','division','biomass')
dat <- dat %>% group_by(ID_lakepulse,name,division) %>% summarize(biomass = sum(biomass)) %>% ungroup

#making a taxonomy data frame to distinguish calanoids from cyclopoids
taxo <- dat %>% filter(division == 'Copepoda') %>% distinct(name)
taxo$order <- 'calanoid'
taxo <- arrange(taxo, name)
taxo$order[12] <- 'poecilostomatoid'
#finding cyclopoids
taxo$order[c(1,5,6,7,8,13,14,15,28,29,30,31,33,36)] <- 'cyclopoid'

dat$genus <- str_remove(dat$name, '\\ .*')
dat$cyc_order <- taxo$order[match(dat$name,taxo$name)]

out.dat <- data.frame()
lakes <- unique(dat$ID_lakepulse)

for(i in 1:length(lakes)){
  
  sub <- dat %>% filter(ID_lakepulse == lakes[i])
  
  #split copepodids based on species, as long as one species was id'd to genus or species
  sub$copepopid <- 'no'
  sub$copepopid[str_detect(sub$name, 'copepodid')] <- 'yes'
  sub.adults <- filter(sub, copepopid == 'no', !is.na(cyc_order))
  sub.copepodid <- filter(sub, copepopid == 'yes')
  if(nrow(sub.adults)>0){
    for(c in unique(sub.copepodid$cyc_order)){
      order <- c
      babybiomass <- sub.copepodid %>% filter(cyc_order == c) %>% pull(biomass)
      prop.dat <- sub.adults %>% filter(cyc_order == c)
      if(nrow(prop.dat) > 0){
        prop.dat$props <- prop.dat$biomass/sum(prop.dat$biomass)
        prop.dat$bm.to.add <- prop.dat$props*babybiomass
        prop.dat$final.bm <- prop.dat$biomass+prop.dat$bm.to.add
        sub$biomass[match(prop.dat$name,sub$name)] <- prop.dat$final.bm
        sub <- filter(sub, name != paste(order,'copepodid'))
      }
    }
  }
  
  #add data to output data frame
  out.dat <- bind_rows(out.dat,sub)
  
}

#adding proportion of each genus for cladoceran split
out.dat$genus <- gsub(" .*$", "", out.dat$name) #A space (), then any character (.) any number of times (*) until the end of the string ($)

#finishing split in excel, too many exceptions
writexl::write_xlsx(out.dat, '~/Desktop/2019.xlsx')

#reformat to wide
dat <- read_xlsx('~/Desktop/2019.xlsx') %>% rename(Lake_ID = 'ID_lakepulse') %>%
  select(Lake_ID,name,biomass) %>% arrange(name) %>% group_by(Lake_ID) %>% 
  pivot_wider(names_from = name, values_from = biomass) %>% arrange(Lake_ID)

dat <- as.data.frame(dat)
dat[is.na(dat)] <- 0

write.csv(dat, '~/Desktop/ALLfinal_grouping2019.csv')
