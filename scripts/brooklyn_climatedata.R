# libraries
library(lubridate)
library(ggplot2)
library(scales)
library(nasapower)
library(dplyr)
library(plotly)
library(readr)

# TEMPERATURE

# using data from brooklyn HQ 
# t <- read.csv("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/supp_data/wtf_hq_awc_ten_min_2023_08_09_17_52_37.csv")
data <- read.table("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/weather/HQ_met_export.dat", header=TRUE, skip=1, sep=",")
data <- data[-c(1:2),]
data$AirTC_Avg <- as.numeric(data$AirTC_Avg)

# extract monthly averages 
# may 2022
m22 <- data[grepl("2022-05", data$TimeStamp), ]
m22_avg <- mean(m22$AirTC_Avg, na.rm = TRUE)

# november 2022
n22 <- data[grepl("2022-11", data$TimeStamp), ]
n22_avg <- mean(n22$AirTC_Avg, na.rm = TRUE)

# august 2023
a23 <- data[grepl("2023-08", data$TimeStamp), ]
a23_avg <- mean(a23$AirTC_Avg, na.rm = TRUE)

# february 2024 - missing because the cyclone destroyed the station: fill with NASA POWER data

# download data from NASA POWER
sites_lonlat <- data_frame(site = c("STCK")) %>%
  mutate( longitude = case_when(site == "STCK" ~ 145.2406),
          latitude = case_when(site == "STCK" ~ -16.61158))

POWER_dat <- data.frame()
for(x in 1:1){
  P1 <- get_power(community = "sb",
                  lonlat = c(sites_lonlat[[x,2]],sites_lonlat[[x,3]]),
                  pars = c("T2M","RH2M","PRECTOTCORR","PS","WS2M","ALLSKY_SFC_SW_DWN","CLRSKY_SFC_SW_DWN","SZA","ALLSKY_SFC_LW_DWN"),
                  dates = c("2022-05-01","2024-03-01"),
                  temporal_api = "hourly",
                  time_standard = "LST") %>%
    mutate(site = sites_lonlat[[x,1]])
  POWER_dat <- rbind(POWER_dat,P1)
}

# extract monthly averages 
# may 2022
m22 <- POWER_dat[grepl("2022-05", POWER_dat$YYYYMMDD), ]
m22_avg <- mean(m22$T2M, na.rm = TRUE)

# november 2022
n22 <- POWER_dat[grepl("2022-11", POWER_dat$YYYYMMDD), ]
n22_avg <- mean(n22$T2M, na.rm = TRUE)

# august 2023
a23 <- POWER_dat[grepl("2023-08", POWER_dat$YYYYMMDD), ]
a23_avg <- mean(a23$T2M, na.rm = TRUE)

# february 2024
f24 <- POWER_dat[grepl("2024-02", POWER_dat$YYYYMMDD), ]
f24_avg <- mean(f24$T2M, na.rm = TRUE)

# PRECIPITATION 
# read in 2022-2024 precip data 
p22 <- read.csv("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/weather/MC_precip_22.csv")
p23 <- read.csv("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/weather/MC_precip_23.csv")
p24 <- read.csv("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/weather/MC_precip_24.csv")

# merge into one sheet 
p <- rbind(p22, p23, p24)

# for each campaign, create 3 dfs 
# May 2022
# 16 April - 16 May 
may_1mo <- p[(p$Year == 2022 & p$Month == 4 & p$Day >= 16) | (p$Year == 2022 & p$Month == 5 & p$Day <= 16), ]
may_1mo_sum <- sum(may_1mo$Rainfall.amount..millimetres.)
# 16 March - 16 May 
may_2mo <- p[(p$Year == 2022 & p$Month == 3 & p$Day >= 16) | (p$Year == 2022 & p$Month == 4) | (p$Year == 2022 & p$Month == 5 & p$Day <= 16), ]
may_2mo_sum <- sum(may_2mo$Rainfall.amount..millimetres.)
# 16 Feb - 16 May
may_3mo <- p[(p$Year == 2022 & p$Month == 2 & p$Day >= 16) | (p$Year == 2022 & p$Month == 3) | (p$Year == 2022 & p$Month == 4) | (p$Year == 2022 & p$Month == 5 & p$Day <= 16), ]
may_3mo_sum <- sum(may_3mo$Rainfall.amount..millimetres.)

# Nov 2022
# 14 Oct - 14 Nov 
nov_1mo <- p[(p$Year == 2022 & p$Month == 10 & p$Day >= 14) | (p$Year == 2022 & p$Month == 11 & p$Day <= 14), ]
nov_1mo_sum <- sum(nov_1mo$Rainfall.amount..millimetres.)
# 14 Sept - 14 Nov
nov_2mo <- p[(p$Year == 2022 & p$Month == 9 & p$Day >= 14) | (p$Year == 2022 & p$Month == 10) | (p$Year == 2022 & p$Month == 11 & p$Day <= 14), ]
nov_2mo_sum <-sum(nov_2mo$Rainfall.amount..millimetres.)
# 14 Aug - 14 Nov 
nov_3mo <- p[(p$Year == 2022 & p$Month == 8 & p$Day >= 14) | (p$Year == 2022 & p$Month == 9) |(p$Year == 2022 & p$Month == 10) |(p$Year == 2022 & p$Month == 11 & p$Day <= 14), ]
nov_3mo_sum <- sum(nov_3mo$Rainfall.amount..millimetres.)

# Aug 2023
# 7 July - 7 Aug
aug_1mo <- p[(p$Year == 2023 & p$Month == 7 & p$Day >= 7) | (p$Year == 2023 & p$Month == 8 & p$Day <= 7), ]
aug_1mo_sum <- sum(aug_1mo$Rainfall.amount..millimetres.)
# 7 June - 7 Aug
aug_2mo <- p[(p$Year == 2023 & p$Month == 6 & p$Day >= 7) | (p$Year == 2023 & p$Month == 7) | (p$Year == 2023 & p$Month == 8 & p$Day <= 7), ]
aug_2mo_sum <- sum(aug_2mo$Rainfall.amount..millimetres.)
# 7 May - 7 Aug
aug_3mo <- p[(p$Year == 2023 & p$Month == 5 & p$Day >= 7) | (p$Year == 2023 & p$Month == 6) |(p$Year == 2023 & p$Month == 7) |(p$Year == 2023 & p$Month == 8 & p$Day <= 7), ]
aug_3mo_sum <- sum(aug_3mo$Rainfall.amount..millimetres.)

# Feb 2024
# 26 Jan - 26 Feb 
feb_1mo <- p[(p$Year == 2024 & p$Month == 1 & p$Day >= 26) | (p$Year == 2024 & p$Month == 2 & p$Day <= 26), ]
feb_1mo_sum <- sum(feb_1mo$Rainfall.amount..millimetres.)
# 26 Dec - 26 Feb
feb_2mo <- p[(p$Year == 2023& p$Month == 12 & p$Day >= 26) | (p$Year == 2024 & p$Month == 1) | (p$Year == 2024 & p$Month == 2 & p$Day <= 26), ]
feb_2mo_sum <- sum(feb_2mo$Rainfall.amount..millimetres.)
# 26 Nov - 26 Feb 
feb_3mo <- p[(p$Year == 2023& p$Month == 11 & p$Day >= 26) | (p$Year == 2023 & p$Month == 12) | (p$Year == 2024 & p$Month == 1) | (p$Year == 2024 & p$Month == 2 & p$Day <= 26), ]
feb_3mo_sum <- sum(feb_3mo$Rainfall.amount..millimetres.)

# create precip dataframe 
precip <- data.frame(may_1mo_sum, may_2mo_sum, may_3mo_sum,
                     aug_1mo_sum, aug_2mo_sum, aug_3mo_sum, 
                     nov_1mo_sum, nov_2mo_sum, nov_3mo_sum,
                     feb_1mo_sum, feb_2mo_sum, feb_3mo_sum)

# write.csv(precip,"/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/supp_data/precip.csv", row.names = FALSE)


#################
# EXTRA 
#################
# format datetime
# POWER_dat$Datetime<-with(POWER_dat, ymd_h(paste(YEAR, MM, DD, HR, sep= ' ')))
# 
# POWER_dat <- POWER_dat[POWER_dat$T2M != -999.00, ]
# 
# # graph temp over year 
# p <- ggplot(POWER_dat, aes(x = Datetime, y = T2M))+
#   geom_point(aes(text=MM), colour="red", alpha=1/2)+
#   scale_x_datetime(labels = date_format("%b-%Y"))
# ggplotly(p)