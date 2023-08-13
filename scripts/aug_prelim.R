#TMI flux data load and prep (aug 2023 data)

# workplace setup 
setwd("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/aug23LGR")
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
data_clean2 <- cbind(date_time,Year,Month,Day,fDOY,Hour,Min,Sec,data_s[,-1])
# save LGR data as dt
data_clean2 <- data.table(data_clean2)

# read in file with metadata
meta_TMI2 <- read.csv('/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/TMI_flux_metadata_aug.csv')

as.Date(meta_TMI2$measurement_day_actual, '%m/%d/%Y')

## Pull out date and time data
date_time_start <- strptime(paste(meta_TMI2$measurement_day_actual,meta_TMI2$flux_start),'%m/%d/%Y %H:%M:%S')

## Add 30 seconds to set buffered start time
date_time_start <- date_time_start+30

## Add 90 sec to determine end time
date_time_end <- date_time_start+90

## Add fDOY columns for start and endtime to dataframe
fDOY_start <- as.numeric(julian(date_time_start,'2023-01-01'))  #Change for year
fDOY_end <- as.numeric(julian(date_time_end,'2023-01-01'))  #Change for year
meta_TMI2 <- cbind(fDOY_start,fDOY_end,meta_TMI2)

# calculate average respiration temperature - add column
meta_TMI2$avg_respT <- rowMeans(meta_TMI2[ , c(9:11)], na.rm = TRUE)

# change sample to factor level variables 
meta_TMI2$sample <- as.factor(meta_TMI2$sample) 

# rename co2 and ch4 columns because naming format including '.' messes up the sql code 
data_clean2 <- data_clean2 %>%
  dplyr::rename(CH4 = X.CH4._ppm,
         CO2 = X.CO2._ppm,
         CH4_dry = X.CH4.d_ppm,
         CO2_dry = X.CO2.d_ppm)

## 2. Merge LGR files with metadata 

# complete merge for respiration data and metadata
data_merge2 <- sqldf("select data_clean2.fDOY, data_clean2.CO2_dry, data_clean2.CH4_dry, meta_TMI2.sample, meta_TMI2.avg_respT, meta_TMI2.flux_source, meta_TMI2.position, meta_TMI2.tag, meta_TMI2.directon, meta_TMI2.chamber, meta_TMI2.chamber_SA_cm2, meta_TMI2.fDOY_start 
from data_clean2 LEFT JOIN meta_TMI2 ON (data_clean2.fDOY BETWEEN meta_TMI2.fDOY_start AND meta_TMI2.fDOY_end)")

# use 'complete.cases'
data_merge2 <- data_merge2[complete.cases(data_merge2), ]

# check to make sure all of the samples are there 
summary(unique(data_merge2$sample))

## Visually check clipped bits of data 

# co2 (xlim=c(0.4247222, 0.4332870) corresponds to only MD1 first file)
par(mfrow=c(2,1),mar=c(4,4,1,1))
with(data_clean2,plot(fDOY,CO2_dry,ylim=c(400,1300), xlim=c(218.4970, 218.5038)))
with(data_merge2,plot(fDOY,CO2_dry,ylim=c(400,1300), xlim=c(218.4970, 218.5038)))

# ch4 
par(mfrow=c(2,1),mar=c(4,4,1,1))
with(data_clean2,plot(fDOY,CH4_dry,ylim=c(1, 10), xlim=c(218.4970, 218.5038)))
with(data_merge2,plot(fDOY,CH4_dry,ylim=c(1,10), xlim=c(218.4970, 218.5038)))

# use known volume of chamber (in m3) to convert ppm of gas to liters of gas

# small
# h = 5.7 cm
# d = 7.8 cm
sm <- pi*((7.8/2)^2)*5.7

# medium
# h = 10.8 cm
# d = 10.3 cm
med <- pi*((10.3/2)^2)*10.8

# large
# h = 16.1 cm
# d = 15.3 cm
lg <- pi*((15.3/2)^2)*16.1

# mini
# v = 119 ml
mini <- 119

volume_cm3 <- c(sm, med, lg, mini)
chamber <- c("sm", "med", "lg", "mini")
chamber_vols <- data.frame(volume_cm3, chamber)

# chamber volumes are in cm3 - convert to m3 for the next step 
chamber_vols$volume_m3 <- chamber_vols$volume_cm3*0.000001

data_merge2 <- merge(data_merge2, chamber_vols, by = c("chamber")) 
data_merge2 <- data_merge2 %>% relocate(chamber, .after = chamber_SA_cm2)

# also convert SA chamber to be in m2 (currently in cm2)
data_merge2$chamber_SA_m2 <- data_merge2$chamber_SA_cm2 * 0.0001
data_merge2 <- data_merge2 %>% relocate(chamber_SA_m2, .after = chamber_SA_cm2)

# CH4
data_merge2$CH4_dry_L=
  # parts CH4 per million parts air * volume of air in chamber (m3) * 1000 L per m3
  (data_merge2$CH4_dry/1000000) * data_merge2$volume_m3 * 1000 

# CO2
data_merge2$CO2_dry_L=
  # parts CO2 per million parts air * volume of air in chamber (m3) * 1000 L per m3
  (data_merge2$CO2_dry/1000000) * data_merge2$volume_m3 * 1000 

# Use ideal gas law to calculate umol of CH4 or mmol of CO2
# CH4
data_merge2$CH4_dry_umol=
  # (atm pressure * L CH4) / (R in L*atm/?K*mol * ?K temp) * 10^6 umol/mol
  ((1*data_merge2$CH4_dry_L)/(0.08206*(data_merge2$avg_respT+273)))*10^6

# CO2
data_merge2$CO2_dry_mmol=
  # (atm pressure * L CO2) / (R in L*atm/?K*mol * ?K temp) * 10^3 mmol/mol
  ((1*data_merge2$CO2_dry_L)/(0.08206*(data_merge2$avg_respT+273)))*10^3

## Calculate chamber fluxes
## Identify start of fluxes
Time_start2 <- meta_TMI2[,c("fDOY_start","sample")]
flux.times2=unique(Time_start2$fDOY_start)

# Create new dataframes to hold final fluxes 
fluxes.CH4_aug=data.frame(matrix(NA,nrow=length(flux.times2)))
fluxes.CO2_aug=data.frame(matrix(NA,nrow=length(flux.times2)))

## Add named columns
# CH4
fluxes.CH4_aug$Time_start=flux.times2
fluxes.CH4_aug$flux.CH4=0
fluxes.CH4_aug$R2.CH4=0
fluxes.CH4_aug$p.CH4=0

# CO2
fluxes.CO2_aug$Time_start=flux.times2
fluxes.CO2_aug$flux.CO2=0
fluxes.CO2_aug$R2.CO2=0
fluxes.CO2_aug$p.CO2=0

## Remove initial empty column
fluxes.CO2_aug=fluxes.CO2_aug[,-1]
fluxes.CH4_aug=fluxes.CH4_aug[,-1]

## For each start time (aug 23)
for (i in flux.times2) {
  ## CH4 ##
  # Subset data for one chamber measurement
  temp1=subset(data_merge2,fDOY_start==i)
  # Set corresponding row of output table
  j=which(flux.times2==i)
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
  temp2=subset(data_merge2,fDOY_start==i)
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

## Finalize fluxes, merge back with metadata and extract

# merge metadata back in to flux files 
# CO2
meta_TMI <- meta_TMI %>%
  rename(Time_start = fDOY_start)
meta_TMI2 <- meta_TMI2 %>%
  rename(Time_start = fDOY_start)

CO2_fluxfinal_may22 <- merge(fluxes.CO2_may, meta_TMI, by="Time_start")
CO2_fluxfinal_may22$campaign <- "may22"
CO2_fluxfinal_nov22 <- merge(fluxes.CO2_nov, meta_TMI2, by="Time_start")
CO2_fluxfinal_nov22$campaign <- "nov22"

CO2_fluxfinal_all <- rbind.fill(CO2_fluxfinal_may22, CO2_fluxfinal_nov22)

# CH4
CH4_fluxfinal_may22 <- merge(fluxes.CH4_may, meta_TMI, by="Time_start")
CH4_fluxfinal_may22$campaign <- "may22"
CH4_fluxfinal_nov22 <- merge(fluxes.CH4_nov, meta_TMI2, by="Time_start")
CH4_fluxfinal_nov22$campaign <- "nov22"

CH4_fluxfinal_all <- rbind.fill(CH4_fluxfinal_may22, CH4_fluxfinal_nov22)

## Export fluxes as .csv file
# CO2 (units are mmol/day)
# write.csv(CO2_fluxfinal_all,"/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/finalfluxes/CO2_fluxfinal_all.csv", row.names = FALSE)

# CH4 (units are umol/day)
# write.csv(CH4_fluxfinal_all,"/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/finalfluxes/CH4_fluxfinal_all.csv", row.names = FALSE)
```