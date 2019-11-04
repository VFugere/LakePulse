library(tidyverse)
library(skimr)
library(readxl)

dat <- read.csv2('~/Desktop/LakePulse2017_kestrel_QC.csv')
skim(dat)

dat <- read.csv2('~/Desktop/LakePulse2017_RBR_top_bottom_QC.csv')
skim(dat)

dat <- read_xlsx('~/Desktop/hydroLakes.xlsx')
skim(dat)

HydroLakes_Canada <- filter(dat, Country == 'Canada')

