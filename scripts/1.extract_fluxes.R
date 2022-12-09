# workspace prep
library(data.table)
library(scales)
library(ggplot2)
library(ggpubr)
library(sqldf)
library(dplyr)
library(tidyverse)
library(stringr)

# TMIflux data load and prep ----

# get names of all LGR files with GHG concentration, a.k.a the data with the form ('...f####.txt')
# read in LGR datafiles to list and go from .txt to .csv 
# add in file name (f####) as an identifier column to merge with metadata

dataTMI <- list.files(path ="/Users/abbeyyatsko/Desktop/repos/TMI_flux/data", pattern='f00', full.names= T) %>%
  setNames(nm = .) %>% 
  lapply(read.csv,skip=1) %>%
  bind_rows(.id = "file") %>%
  mutate(file = str_extract(file, "f00\\d+"))

# there is a pgp error message at the end of some LGR files. need to remove these lines 
# remove any entry where there is NA for both x.CH4._ppm and x.CH4._ppm_sd
data_clean <- subset(dataTMI, X.CH4._ppm!="NA")

# read in metadata to try and change the year for each f#### file so that the fDOY makes sense ----
meta_TMI <- read.csv("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/TMI_flux_metadata - Sheet1.csv")
file <- unique(meta_TMI$file)
file_years <- as.data.frame(file)
file_years$year_alt <- 2002:2068

# merge in alternative year (year_alt)
data_clean <- merge(data_clean, file_years, by = c("file")) 

# move year_alt 
data_clean <- data_clean %>% relocate(year_alt, .before = X.CH4._ppm)

# replace original year (always 2002) in date_time with the alternative year (needed for fDOY)

# TEST
# d <- data_clean[c("Time", "year_alt")]
# d$Time <- gsub("2002", "", d$Time)
# d$time_new <- with(d, paste(substr(Time, 1, 8),
#                           year_alt,
#                           substr(Time, 9, nchar(Time)),
#                           sep = ""))

# change in actual df (data_clean)

# remove 2002 year from Time column so that unique year can be pasted in
data_clean$Time <- gsub("2002", "", data_clean$Time)

# create altered time column by pasting in alternative year to Time column in the place of the old 2002
data_clean$time_alt <- with(data_clean, paste(substr(Time, 1, 8),
                           year_alt,
                           substr(Time, 9, nchar(Time)),
                           sep = ""))
# move up time_alt to compare with original Time
data_clean <- data_clean %>% relocate(time_alt, .before = year_alt)

# pull out date and time data - use time_alt
date_time <- strptime(data_clean[,3],format='%m/%d/%Y %H:%M:%S')

# formatting dates and times ----
# add year,month,day,JD,hour,min,sec columns to dataframe
Year <- as.numeric(format(date_time,'%Y'))
Month <- as.numeric(format(date_time,'%m'))
Day <- as.numeric(format(date_time,'%d'))
fDOY <- as.numeric(julian(date_time,'2002-01-01'))  #Change for year
Hour <- as.numeric(format(date_time,'%k'))
Min <- as.numeric(format(date_time,'%M'))
Sec <- as.numeric(format(date_time,'%S'))
data_clean <- cbind(date_time,Year,Month,Day,fDOY,Hour,Min,Sec,data_clean[,-3])
# save LGR data as data.table
data_clean <- data.table(data_clean)

# preliminary viz ----
# for co2
co2 <- ggplot(data_clean, aes(date_time, X.CO2.d_ppm)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("CO2 flux")

# for ch4 
ch4 <- ggplot(data_clean, aes(date_time, X.CH4.d_ppm)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("CH4 flux")

# looks like there are the right number of lines that represent the 67 flux files 

# merge with metadata ----

# merge in year_alt to metadata
meta_TMI <- merge(meta_TMI, file_years, by = c("file")) 

# add in column for measurement_day_alt
meta_TMI$measurement_day_alt <- "01/01/"

# paste in year_alt to the end of measurement_day_alt
meta_TMI$measurement_day_alt <- with(meta_TMI, paste(substr(measurement_day_alt, 1, 8),
                                              year_alt,
                                              substr(measurement_day_alt, 9, nchar(measurement_day_alt)),
                                              sep = ""))
# reposition measurement_day_alt in df 
meta_TMI <- meta_TMI %>% relocate(measurement_day_alt, .before = measurement_day_actual)
meta_TMI <- meta_TMI %>% relocate(year_alt, .before = measurement_day_actual)

# format measuremetn_day_alt
as.Date(meta_TMI$measurement_day_alt, '%m/%d/%Y')

## Pull out date and time data
date_time_start <- strptime(paste(meta_TMI$measurement_day_alt,meta_TMI$flux_start),'%m/%d/%Y %H:%M:%S')

## Add 30 seconds to set buffered start time
date_time_start <- date_time_start+30

## Add 1.5 min (90 sec) to determine end time
date_time_end <- date_time_start+90

## Add fDOY columns for start and endtime to dataframe
fDOY_start <- as.numeric(julian(date_time_start,'2002-01-01'))  #Change for year
fDOY_end <- as.numeric(julian(date_time_end,'2002-01-01'))  #Change for year
meta_TMI <- cbind(fDOY_start,fDOY_end,meta_TMI)

# calculate average respiration temperature - add column
meta_TMI$avg_respT <- rowMeans(meta_TMI[ , c(12:15)], na.rm = TRUE)

# change sample to factor level variables 
meta_TMI$sample <- as.factor(meta_TMI$sample) 

# rename co2 and ch4 columns because naming format including '.' messes up the sql code 
data_clean <- data_clean %>%
  rename(CH4 = X.CH4._ppm,
         CO2 = X.CO2._ppm, 
         CH4_dry = X.CH4.d_ppm,
         CO2_dry = X.CO2.d_ppm)

# merge LGR files to plot with metadata ---- 
# goal is to clip out samples and remove measurements for when the chamber was empty 

# complete merge for respiration data and metadata
data_merge <- sqldf("select data_clean.fDOY, data_clean.CO2_dry, data_clean.CH4_dry, meta_TMI.sample, meta_TMI.avg_respT, meta_TMI.flux_source, meta_TMI.position, meta_TMI.tag, meta_TMI.directon, meta_TMI.chamber_SA_cm2, meta_TMI.soil_moisture, meta_TMI.fDOY_start 
from data_clean LEFT JOIN meta_TMI ON (data_clean.fDOY BETWEEN meta_TMI.fDOY_start AND meta_TMI.fDOY_end)")

# filter out rows that  have no value for 'sample' - this means that no mound was being recorded
data_merge <- data_merge[complete.cases(data_merge[,c("sample")]),]

# check to make sure all of the samples are there 
summary(unique(data_merge$sample))







# now visually check clipped bits of data 
# co2 (xlim=c(212.4705, 212.5215) corresponds to only the first day of flux measurement)
par(mfrow=c(2,1),mar=c(4,4,1,1))
with(data,plot(fDOY,CO2_dry,ylim=c(400,1300), xlim=c(212.4705, 212.5215)))
with(data_merge,plot(fDOY,CO2_dry,ylim=c(400,1300), xlim=c(212.4705, 212.5215)))

# ch4 
par(mfrow=c(2,1),mar=c(4,4,1,1))
with(data,plot(fDOY,CH4_dry,ylim=c(1.95,2.4), xlim=c(212.4705, 212.5215)))
with(data_merge,plot(fDOY,CH4_dry,ylim=c(1.95,2.4), xlim=c(212.4705, 212.5215)))
```





