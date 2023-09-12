library(data.table)
library(scales)
library(ggplot2)
library(ggpubr)
library(sqldf)
library(dplyr)
library(tidyverse)
library(stringr)
library(plyr)

# workplace setup 
setwd("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/aug23LGR_daily")
filenames <- list.files(pattern='f00',full.names=T)

# read in LGR datafiles to list and go from .txt to .csv 
data_s <- lapply(filenames,read.csv,skip=1)
data_s <- rbindlist(data_s)
data_s <- subset(data_s, X.CH4._ppm!="NA")
as.character(data_s$Time)

date_time <- strptime(data_s$Time,format='%d/%m/%Y %H:%M:%S')

# formatting dates and times ----
# add year,month,day,JD,hour,min,sec columns to dataframe
Year <- as.numeric(format(date_time,'%Y'))
Month <- as.numeric(format(date_time,'%m'))
Day <- as.numeric(format(date_time,'%d'))
fDOY <- as.numeric(julian(date_time,'2023-01-01'))  #Change for year
Hour <- as.numeric(format(date_time,'%k'))
Min <- as.numeric(format(date_time,'%M'))
Sec <- as.numeric(format(date_time,'%S'))
data_clean_day <- cbind(date_time,Year,Month,Day,fDOY,Hour,Min,Sec,data_s[,-1])

# save LGR data as dt
data_clean_day <- data.table(data_clean_day)

meta_TMI_day <- read.csv('/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/TMI_flux_metadata_aug_daily.csv')

as.Date(meta_TMI_day$measurement_day_actual, '%m/%d/%Y')

## Pull out date and time data
date_time_start <- strptime(paste(meta_TMI_day$measurement_day_actual,meta_TMI_day$flux_start),'%m/%d/%Y %H:%M:%S')

## Add 30 seconds to set buffered start time
date_time_start <- date_time_start+30

## Add 90 sec to determine end time
date_time_end <- date_time_start+90

## Add fDOY columns for start and endtime to dataframe
fDOY_start <- as.numeric(julian(date_time_start,'2023-01-01'))  #Change for year
fDOY_end <- as.numeric(julian(date_time_end,'2023-01-01'))  #Change for year
meta_TMI_day <- cbind(fDOY_start,fDOY_end,meta_TMI_day)

# calculate average respiration temperature - add column
meta_TMI_day$avg_respT <- rowMeans(meta_TMI_day[ , c(9:11)], na.rm = TRUE)

# change sample to factor level variables 
meta_TMI_day$sample <- as.factor(meta_TMI_day$sample) 

# rename co2 and ch4 columns because naming format including '.' messes up the sql code 
data_clean_day <- data_clean_day %>%
  dplyr::rename(CH4 = X.CH4._ppm,
                CO2 = X.CO2._ppm,
                CH4_dry = X.CH4.d_ppm,
                CO2_dry = X.CO2.d_ppm)

# complete merge for respiration data and metadata (aug)
data_merge_day <- sqldf("select data_clean_day.fDOY, data_clean_day.CO2_dry, data_clean_day.CH4_dry, meta_TMI_day.sample, meta_TMI_day.avg_respT, meta_TMI_day.flux_source, meta_TMI_day.position, meta_TMI_day.tag, meta_TMI_day.directon, meta_TMI_day.chamber, meta_TMI_day.chamber_SA_cm2, meta_TMI_day.fDOY_start 
from data_clean_day LEFT JOIN meta_TMI_day ON (data_clean_day.fDOY BETWEEN meta_TMI_day.fDOY_start AND meta_TMI_day.fDOY_end)")

# use 'complete.cases'
data_merge_day <- data_merge_day[complete.cases(data_merge_day), ]

# check to make sure all of the samples are there (9 mounds total)
summary(unique(data_merge_day$sample))

# use known volume of chamber (in m3) to convert ppm of gas to liters of gas
# small, h = 5.7 cm, d = 7.8 cm
sm <- pi*((7.8/2)^2)*5.7

# medium, h = 10.8 cm, d = 10.3 cm
med <- pi*((10.3/2)^2)*10.8

# large, h = 16.1 cm, d = 15.3 cm
lg <- pi*((15.3/2)^2)*16.1

# mini, v = 119 ml
mini <- 119

volume_cm3 <- c(sm, med, lg, mini)
chamber <- c("sm", "med", "lg", "mini")
chamber_vols <- data.frame(volume_cm3, chamber)

# chamber volumes are in cm3 - convert to m3 for the next step 
chamber_vols$volume_m3 <- chamber_vols$volume_cm3*0.000001

data_merge_day <- merge(data_merge_day, chamber_vols, by = c("chamber")) 
data_merge_day <- data_merge_day %>% relocate(chamber, .after = chamber_SA_cm2)

# also convert SA chamber to be in m2 (currently in cm2)
data_merge_day$chamber_SA_m2 <- data_merge_day$chamber_SA_cm2 * 0.0001
data_merge_day <- data_merge_day %>% relocate(chamber_SA_m2, .after = chamber_SA_cm2)

# CH4
data_merge_day$CH4_dry_L=
  # parts CH4 per million parts air * volume of air in chamber (m3) * 1000 L per m3
  (data_merge_day$CH4_dry/1000000) * data_merge_day$volume_m3 * 1000 

# CO2
data_merge_day$CO2_dry_L=
  # parts CO2 per million parts air * volume of air in chamber (m3) * 1000 L per m3
  (data_merge_day$CO2_dry/1000000) * data_merge_day$volume_m3 * 1000 

# Use ideal gas law to calculate umol of CH4 or mmol of CO2
# CH4
data_merge_day$CH4_dry_umol=
  # (atm pressure * L CH4) / (R in L*atm/?K*mol * ?K temp) * 10^6 umol/mol
  ((1*data_merge_day$CH4_dry_L)/(0.08206*(data_merge_day$avg_respT+273)))*10^6

# CO2
data_merge_day$CO2_dry_mmol=
  # (atm pressure * L CO2) / (R in L*atm/?K*mol * ?K temp) * 10^3 mmol/mol
  ((1*data_merge_day$CO2_dry_L)/(0.08206*(data_merge_day$avg_respT+273)))*10^3

## Identify start of fluxes
Time_start <- meta_TMI_day[,c("fDOY_start","sample")]
flux.times=unique(Time_start$fDOY_start)

# Create new dataframes to hold final fluxes 
fluxes.CH4_aug=data.frame(matrix(NA,nrow=length(flux.times)))
fluxes.CO2_aug=data.frame(matrix(NA,nrow=length(flux.times)))

## Add named columns
# CH4
fluxes.CH4_aug$Time_start=flux.times
fluxes.CH4_aug$flux.CH4=0
fluxes.CH4_aug$R2.CH4=0
fluxes.CH4_aug$p.CH4=0

# CO2
fluxes.CO2_aug$Time_start=flux.times
fluxes.CO2_aug$flux.CO2=0
fluxes.CO2_aug$R2.CO2=0
fluxes.CO2_aug$p.CO2=0

## Remove initial empty column
fluxes.CO2_aug=fluxes.CO2_aug[,-1]
fluxes.CH4_aug=fluxes.CH4_aug[,-1]

## For each start time (aug 2023)
for (i in flux.times) {
  ## CH4 ##
  # Subset data for one chamber measurement
  temp1=subset(data_merge_day,fDOY_start==i)
  # Set corresponding row of output table
  j=which(flux.times==i)
  # Determine if start time has a CH4 flux
  if (nrow(temp1)>0) {
    # If so:  
    # Calulate flux in umol/day using linear regression
    mod=with(temp1,lm(CH4_dry_umol~fDOY))
    # Save flux rate and R2 and p-value of slope in corresponding row of dataframe
    # flux rate
    fluxes.CH4_aug$flux.CH4[j]=coef(mod)[2]/temp1$chamber_SA_m2 # chamber surface area
    # R2 of slope
    fluxes.CH4_aug$R2.CH4[j]=summary(mod)$r.squared
    # p-value of slope
    fluxes.CH4_aug$p.CH4[j]=summary(mod)$coefficients[2,4]
    # If not:
    # Fill rows of table with NA    
  } else {
    fluxes.CH4_aug$flux.CH4[j]=NA
    fluxes.CH4_aug$R2.CH4[j]=NA
    fluxes.CH4_aug$p.CH4[j]=NA
  }
  ## CO2 ##
  # Subset data for one chamber measurement
  temp2=subset(data_merge_day,fDOY_start==i)
  # Calulate flux in mmol/day using linear regression
  mod=with(temp2,lm(CO2_dry_mmol~fDOY))
  # Save flux rate and R2 and p-value of slope in corresponding row of dataframe
  # flux rate
  fluxes.CO2_aug$flux.CO2[j]=coef(mod)[2]/temp2$chamber_SA_m2 # chamber surface area
  # R2 of slope
  fluxes.CO2_aug$R2.CO2[j]=summary(mod)$r.squared
  # p-value of slope
  fluxes.CO2_aug$p.CO2[j]=summary(mod)$coefficients[2,4]
}

# merge metadata back in to flux files 
# CO2
meta_TMI_day <- meta_TMI_day %>%
  dplyr::rename(Time_start = fDOY_start)

CO2_fluxfinal_aug23_day <- merge(fluxes.CO2_aug, meta_TMI_day, by="Time_start")
CO2_fluxfinal_aug23_day$campaign <- "aug23"

# CH4
CH4_fluxfinal_aug23_day <- merge(fluxes.CH4_aug, meta_TMI_day, by="Time_start")
CH4_fluxfinal_aug23_day$campaign <- "aug23"

# merge in species data to flux data 
species <- read.csv("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/finalfluxes/daily_species.csv")
species <- species[c("X", "ID_cleaned")]

# merge sp and CH4flux_final, CO2_fluxfinal
colnames(species)[1] = "sample"

CH4_fluxfinal <- merge(CH4_fluxfinal_aug23_day, species, by = c("sample")) 
CO2_fluxfinal <- merge(CO2_fluxfinal_aug23_day, species, by = c("sample")) 

# new column for methane:co2 -first need to convert CO2 from mmol to umol (in order to match CH4 in units of umol)
# (conversion factor 1000umol = 1mmol)
CO2_fluxfinal$flux.CO2_umol <- CO2_fluxfinal$flux.CO2*1000
CO2_fluxfinal <- CO2_fluxfinal %>% relocate(flux.CO2_umol, .after = flux.CO2)

# take out just CO2 flux to paste in working CH4_fluxfinal df
a <- CO2_fluxfinal[c("Time_start", "flux.CO2_umol")]
CH4_fluxfinal <- merge(CH4_fluxfinal, a, by = c("Time_start")) 
CH4_fluxfinal <- CH4_fluxfinal %>% relocate(flux.CO2_umol, .after = flux.CH4)

CH4_fluxfinal$ch4_co2 <-  CH4_fluxfinal$flux.CH4/CH4_fluxfinal$flux.CO2_umol
CH4_fluxfinal <- CH4_fluxfinal %>% relocate(ch4_co2, .after = flux.CH4)

# create categories for time of day 


ggplot(data = CH4_fluxfinal, aes(x = Time_start, y = flux.CH4, color = sample)) + 
  geom_point()+ 
  theme_classic()+
  ylab("CH4 flux umol/d/m2")+
  facet_wrap(~ID_cleaned, scales = "free")+
  ylim(0,30000)


