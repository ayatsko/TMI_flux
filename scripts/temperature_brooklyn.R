# libraries
library(lubridate)
library(ggplot2)
library(scales)
library(nasapower)
library(dplyr)
library(plotly)

# download data from NASA POWER
library(nasapower)
sites_lonlat <- data_frame(site = c("STCK")) %>%
  mutate( longitude = case_when(site == "STCK" ~ 145.2406),
          latitude = case_when(site == "STCK" ~ -16.61158))

POWER_dat <- data.frame()
for(x in 1:1){
  P1 <- get_power(community = "sb",
                  lonlat = c(sites_lonlat[[x,2]],sites_lonlat[[x,3]]),
                  pars = c("T2M","RH2M","PRECTOTCORR","PS","WS2M","ALLSKY_SFC_SW_DWN","CLRSKY_SFC_SW_DWN","SZA","ALLSKY_SFC_LW_DWN"),
                  dates = c("2018-01-01","2023-03-1"),
                  temporal_api = "hourly",
                  time_standard = "LST") %>%
    mutate(site = sites_lonlat[[x,1]])
  POWER_dat <- rbind(POWER_dat,P1)
}

# format datetime
POWER_dat$Datetime<-with(POWER_dat, ymd_h(paste(YEAR, MM, DD, HR, sep= ' ')))

POWER_dat <- POWER_dat[POWER_dat$T2M != -999.00, ]

# graph temp over year 
p <- ggplot(POWER_dat, aes(x = Datetime, y = T2M))+
  geom_point(aes(text=MM), colour="red", alpha=1/2)+
  scale_x_datetime(labels = date_format("%b-%Y"))
ggplotly(p)




