# Load packages and set working directory
library(tidyverse)
library(lubridate)
library(data.table)

# read in IRGA files
setwd("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/naked_termite/IRGA")
filenames <- list.files(pattern='f000',full.names=T)

# read in LGR datafiles to list and go from .txt to .csv 
dat <- lapply(filenames,read.csv,skip=1)
dat <- rbindlist(dat)
dat <- subset(dat, X.CH4._ppm!="NA")
as.character(dat$Time)

date_time <- strptime(dat$Time,format='%d/%m/%Y %H:%M:%S')

# add year,month,day,JD,hour,min,sec columns to dataframe
Year <- as.numeric(format(date_time,'%Y'))
Month <- as.numeric(format(date_time,'%m'))
Day <- as.numeric(format(date_time,'%d'))
fDOY <- as.numeric(julian(date_time,'2023-01-01'))  #Change for year
Hour <- as.numeric(format(date_time,'%k'))
Min <- as.numeric(format(date_time,'%M'))
Sec <- as.numeric(format(date_time,'%S'))
data_clean <- cbind(date_time,Year,Month,Day,fDOY,Hour,Min,Sec,dat)

# save LGR data as dt
data_clean <- data.table(data_clean)

# Preliminary visualization - there should be two measurement times for mid and end of Aug, field and lab campaigns
co2 <- ggplot(data_clean, aes(date_time, X.CO2.d_ppm)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("CO2 flux")

ch4 <- ggplot(data_clean, aes(date_time, X.CH4.d_ppm)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("CH4 flux")

# read in chamber files 
setwd("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/naked_termite/arduino_output")
chamb <- list.files(path = "/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/naked_termite/arduino_output", pattern = "*.csv") %>% 
  map_df(~fread(.))

# format dates and times
chamber_out <- chamb %>%
  unite(year, month, day, hour, minute, second, col = Time) %>%
  mutate(Time = floor_date(ymd_hms(Time), "second")) %>%
  distinct(Time, .keep_all = TRUE)

# merge two files to only include relevant data
data_clean <- force_tz(data_clean, "UTC")
test <- inner_join(data_clean, chamber_out, by = c("date_time" = "Time"))

# check to make sure all of the samples are there (13 total)
summary(unique(test$sample))

### WORKING HERE
# Convert ppm to moles for CH4 and CO2
# chamber volume (in m3) to convert ppm of gas to liters of gas
chamber_v <- 1

# chamber volumes are in cm3 - convert to m3 for the next step 
chamber_vols$volume_m3 <- chamber_vols$volume_cm3*0.000001

data_merge2 <- merge(data_merge2, chamber_vols, by = c("chamber")) 
data_merge2 <- data_merge2 %>% relocate(chamber, .after = chamber_SA_cm2)

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

# Calculate chamber fluxes
## Identify start of fluxes
Time_start2 <- meta_TMI2[,c("fDOY_start","sample")]
flux.times2=unique(Time_start2$fDOY_start)

# Create new dataframes to hold final fluxes 
fluxes.CH4_nov=data.frame(matrix(NA,nrow=length(flux.times2)))
fluxes.CO2_nov=data.frame(matrix(NA,nrow=length(flux.times2)))

## Add named columns
# CH4
fluxes.CH4_nov$Time_start=flux.times2
fluxes.CH4_nov$flux.CH4=0
fluxes.CH4_nov$R2.CH4=0
fluxes.CH4_nov$p.CH4=0

# CO2
fluxes.CO2_nov$Time_start=flux.times2
fluxes.CO2_nov$flux.CO2=0
fluxes.CO2_nov$R2.CO2=0
fluxes.CO2_nov$p.CO2=0

## Remove initial empty column
fluxes.CO2_nov=fluxes.CO2_nov[,-1]
fluxes.CH4_nov=fluxes.CH4_nov[,-1]

## For each start time (nov 2022)
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
    fluxes.CH4_nov$flux.CH4[j]=coef(mod)[2]/temp1$chamber_SA_m2 # chamber surface area
    # R2 of slope
    fluxes.CH4_nov$R2.CH4[j]=summary(mod)$r.squared
    # p-value of slope
    fluxes.CH4_nov$p.CH4[j]=summary(mod)$coefficients[2,4]
    # If not:
    # Fill rows of table with NA    
  } else {
    fluxes.CH4_nov$flux.CH4[j]=NA
    fluxes.CH4_nov$R2.CH4[j]=NA
    fluxes.CH4_nov$p.CH4[j]=NA
  }
  ## CO2 ##
  # Subset data for one chamber measurement
  temp2=subset(data_merge2,fDOY_start==i)
  # Calulate flux in mmol/day using linear regression
  mod=with(temp2,lm(CO2_dry_mmol~fDOY))
  # Save flux rate and R2 and p-value of slope in corresponding row of dataframe
  # flux rate
  fluxes.CO2_nov$flux.CO2[j]=coef(mod)[2]/temp2$chamber_SA_m2 # chamber surface area
  # R2 of slope
  fluxes.CO2_nov$R2.CO2[j]=summary(mod)$r.squared
  # p-value of slope
  fluxes.CO2_nov$p.CO2[j]=summary(mod)$coefficients[2,4]
}


