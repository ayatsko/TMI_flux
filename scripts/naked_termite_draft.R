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

# rename variables for easier use 
data_merge <- test %>%
  dplyr::rename(CH4 = X.CH4._ppm,
                CO2 = X.CO2._ppm,
                CH4_dry = X.CH4.d_ppm,
                CO2_dry = X.CO2.d_ppm)

# look at where data measurements are - try to figure out which correspond to samples
data_merge$location <- ifelse(grepl("L", data_merge$sample), "lab", "field")
field <- data_merge[data_merge$location == "field",]
lab <- data_merge[data_merge$location == "lab",]

# field measurements: kinda crap 
a <- ggplot(field, aes(date_time, CO2_dry)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("CO2 flux")
ggplotly(a)

b <- ggplot(field, aes(date_time, CH4_dry)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("CH4 flux")
ggplotly(b)

# lab measurements: 
c <- ggplot(lab, aes(date_time, CO2_dry)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("CO2 flux")
ggplotly(c)

d <- ggplot(lab, aes(date_time, CH4_dry)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("CH4 flux")
ggplotly(d)

# clean beginning and end of measurements - clip off negative measurements and then maybe last few seconds? 
# data_merge <- data_merge[data_merge$Etime >= 0, ]
# data_merge <- data_merge[ data_merge$Etime <= 150, ]

# unit conversions for gas flux calculation

# collar height = 5 or 9 cm, offset = 1 cm 
# collar area = 317.8 cm^2, chamber volume = 4076.1 cm^3 (top of chamber space, not in the PVC collar)
# Vc (volume of chamber) = (((collar-1)*317.8) + 4076.1)/1000 L
chamber_v_cm3 <-  (((5-1)*317.8) + 4076.1)

# chamber volumes are in cm3 - convert to m3 for the next step 
chamber_v_m3 <- chamber_v_cm3*0.000001
data_merge$chamber_v_m3 <- chamber_v_m3

# CH4
data_merge$CH4_dry_L=
  # parts CH4 per million parts air * volume of air in chamber (m3) * 1000 L per m3
  (data_merge$CH4_dry/1000000) * data_merge$chamber_v_m3 * 1000 

# CO2
data_merge$CO2_dry_L=
  # parts CO2 per million parts air * volume of air in chamber (m3) * 1000 L per m3
  (data_merge$CO2_dry/1000000) * data_merge$chamber_v_m3 * 1000

# Use ideal gas law to calculate umol of CH4 or mmol of CO2
# CH4
data_merge$CH4_dry_umol=
  # (atm pressure * L CH4) / (R in L*atm/?K*mol * ?K temp) * 10^6 umol/mol
  ((1*data_merge$CH4_dry_L)/(0.08206*(data_merge$Tcham+273)))*10^6

# CO2
data_merge$CO2_dry_mmol=
  # (atm pressure * L CO2) / (R in L*atm/?K*mol * ?K temp) * 10^3 mmol/mol
  ((1*data_merge$CO2_dry_L)/(0.08206*(data_merge$Tcham+273)))*10^3

# Calculate chamber fluxes
## Identify start of fluxes
md <- read.csv("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/naked_termite/naked_termite_metadata.csv")
md$Time_start <- as.POSIXct(paste(md$measurement_day_actual, md$flux_start), format="%Y-%m-%d %H:%M:%S")
md$Time_start <- force_tz(md$Time_start, "UTC")

Time_start <- md[,c("Time_start","sample")]
sample=unique(Time_start$sample)

# Create new dataframes to hold final fluxes 
fluxes.CH4_nt=data.frame(matrix(NA,nrow=length(sample)))
fluxes.CO2_nt=data.frame(matrix(NA,nrow=length(sample)))

## Add named columns
# CH4
fluxes.CH4_nt$sample=sample
fluxes.CH4_nt$flux.CH4=0
fluxes.CH4_nt$R2.CH4=0
fluxes.CH4_nt$p.CH4=0

# CO2
fluxes.CO2_nt$sample=sample
fluxes.CO2_nt$flux.CO2=0
fluxes.CO2_nt$R2.CO2=0
fluxes.CO2_nt$p.CO2=0

## Remove initial empty column
fluxes.CO2_nt=fluxes.CO2_nt[,-1]
fluxes.CH4_nt=fluxes.CH4_nt[,-1]

## For each start time (nov 2022)
for (i in sample) {
  ## CH4 ##
  # Subset data for one chamber measurement
  temp1=subset(data_merge,sample==i)
  # Set corresponding row of output table
  j=which(sample==i)
  # Determine if start time has a CH4 flux
  if (nrow(temp1)>0) {
    # If so:  
    # Calulate flux in umol/day using linear regression
    mod=with(temp1,lm(CH4_dry_umol~fDOY))
    # Save flux rate and R2 and p-value of slope in corresponding row of dataframe
    # flux rate
    fluxes.CH4_nt$flux.CH4[j]=coef(mod)[2]
    # R2 of slope
    fluxes.CH4_nt$R2.CH4[j]=summary(mod)$r.squared
    # p-value of slope
    fluxes.CH4_nt$p.CH4[j]=summary(mod)$coefficients[2,4]
    # If not:
    # Fill rows of table with NA    
  } else {
    fluxes.CH4_nt$flux.CH4[j]=NA
    fluxes.CH4_nt$R2.CH4[j]=NA
    fluxes.CH4_nt$p.CH4[j]=NA
  }
  ## CO2 ##
  # Subset data for one chamber measurement
  temp2=subset(data_merge,sample==i)
  # Calulate flux in mmol/day using linear regression
  mod=with(temp2,lm(CO2_dry_mmol~fDOY))
  # Save flux rate and R2 and p-value of slope in corresponding row of dataframe
  # flux rate
  fluxes.CO2_nt$flux.CO2[j]=coef(mod)[2]
  # R2 of slope
  fluxes.CO2_nt$R2.CO2[j]=summary(mod)$r.squared
  # p-value of slope
  fluxes.CO2_nt$p.CO2[j]=summary(mod)$coefficients[2,4]
}

# these look pretty bad - inspect 
test <- data_merge[data_merge$sample == "ami1_F", ]
plot(test$fDOY, test$CH4_dry_umol)

