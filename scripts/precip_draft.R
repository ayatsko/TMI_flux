# read in 2022-2024 precip data 
p22 <- read.csv("/Users/abbeyyatsko/Desktop/precip_MC/MC_precip_22.csv")
p23 <- read.csv("/Users/abbeyyatsko/Desktop/precip_MC/MC_precip_23.csv")
p24 <- read.csv("/Users/abbeyyatsko/Desktop/precip_MC/MC_precip_24.csv")

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




