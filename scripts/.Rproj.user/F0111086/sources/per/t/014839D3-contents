# WORKSPACE PREP
library(data.table)
library(scales)
library(ggplot2)
library(ggpubr)
library(sqldf)
library(dplyr)
library(tidyverse)
library(stringr)

# TMIflux data load and prep 

# get names of all LGR files with GHG concentration, a.k.a the data with the form ('...f####.txt')
# read in LGR datafiles to list and go from .txt to .csv 
# add in file name (f####) as an identifier column to merge with metadata
dataTMI <- list.files(path ="/Users/abbeyyatsko/Desktop/TMIflux_LGR", pattern='f00', full.names= T) %>%
  setNames(nm = .) %>% 
  lapply(read.csv,skip=1) %>%
  bind_rows(.id = "file") %>%
  mutate(file = str_extract(file, "f00\\d+"))

# there is a pgp error message at the end of some LGR files. need to remove these lines 
# remove any entry where there is NA for both x.CH4._ppm and x.CH4._ppm_sd
data_clean <- subset(dataTMI, X.CH4._ppm!="NA")

# read in metadata to try and change the year for each f#### file so that the fDOY makes sense
meta_TMI <- read.csv("/Users/abbeyyatsko/Desktop/TMIflux_LGR/TMI_flux_metadata - Sheet1.csv")
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

# formatting dates and times 
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

# preliminary viz
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

# view co2 and ch4 plots together
ggarrange(co2, ch4)

