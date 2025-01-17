# species ranges 

# clean datasheets from gbif.org 

c <- read.delim("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/supp_data/range_c.acinaciformis.txt")
n <- read.delim("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/supp_data/range_n.magnus.txt")
a <- read.delim("/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/supp_data/range_a.laurensis.txt")

# coptotermes
c_range <- c[,c("decimalLatitude", "decimalLongitude", "scientificName", "countryCode", "basisOfRecord")]

# remove any rows with na for lat or long
c_range <- c_range[!is.na(c_range$decimalLatitude),]

# remove longitude = 0 
c_range <- c_range[c_range$decimalLongitude != 0,]

# australia only 
c_range <- c_range[c_range$countryCode == "AU",]

# keep rows with scientificName = "Coptotermes acinaciformis (Froggatt, 1898)
c_range <- c_range[c_range$scientificName %in% c("Coptotermes acinaciformis (Froggatt, 1898)", 
                                                 "Coptotermes acinaciformis acinaciformis", 
                                                 "Coptotermes acinaciformis raffrayi (Wasmann, 1900)"), ]

# write.csv(c_range, "/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/supp_data/c_acinaciformis_range_cleaned.csv", row.names = FALSE)

# nasutitermes
n_range <- n[,c("decimalLatitude", "decimalLongitude", "scientificName", "countryCode", "basisOfRecord")]

# remove any rows with na for lat or long
n_range <- n_range[!is.na(n_range$decimalLatitude),]

# australia only 
n_range <- n_range[n_range$countryCode == "AU",]

# write.csv(n_range, "/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/supp_data/n_magnus_range_cleaned.csv", row.names = FALSE)

# amitermes
a_range <- a[,c("decimalLatitude", "decimalLongitude", "scientificName", "countryCode", "basisOfRecord")]

# remove any rows with na for lat or long
a_range <- a_range[!is.na(a_range$decimalLatitude),]

# australia only 
a_range <- a_range[a_range$countryCode == "AU",]

# write.csv(a_range, "/Users/abbeyyatsko/Desktop/repos/TMI_flux/data/supp_data/a_laurensis_range_cleaned.csv", row.names = FALSE)







# gbif range package
library(gbif.range)
library(terra)
library(rnaturalearth)
library(raster)
library(gridExtra)

# plot occurrences
countries = vect(ne_countries(type = "countries",returnclass = "sf"))
plot(countries,col = "#bcbddc")
points(c_range[,c("decimalLongitude","decimalLatitude")],pch=20,col="#99340470",cex=1.5)

plot(countries,col = "#bcbddc")
points(a_range[,c("decimalLongitude","decimalLatitude")],pch=20,col="darkgreen",cex=1.5)

plot(countries,col = "#bcbddc")
points(n_range[,c("decimalLongitude","decimalLatitude")],pch=20,col="darkblue",cex=1.5)

# range estimate 
eco_terra = read_bioreg(bioreg_name = "eco_terra", save_dir = NULL)

range_c = get_range(occ_coord = c_range,
                        bioreg = eco_terra,
                        bioreg_name = "ECO_NAME")

range_a = get_range(occ_coord = a_range,
                    bioreg = eco_terra,
                    bioreg_name = "ECO_NAME")

range_n = get_range(occ_coord = n_range,
                    bioreg = eco_terra,
                    bioreg_name = "ECO_NAME")

plot(countries,col = "#bcbddc", main = "C. acinaciformis distribution", xlim = c(112, 155), ylim = c(-45, -10))
plot(range_c[[2]],col = "darkred",add = TRUE,axes = FALSE,legend = FALSE)

plot(countries,col = "#bcbddc", main = "A. laurensis distribution", xlim = c(112, 155), ylim = c(-45, -10))
plot(range_a[[2]],col = "darkgreen",add = TRUE,axes = FALSE,legend = FALSE)

plot(countries,col = "#bcbddc", main = "N. magnus distribution", , xlim = c(112, 155), ylim = c(-45, -10))
plot(range_n[[2]],col = "darkblue",add = TRUE,axes = FALSE,legend = FALSE)


# area of range (in m2)
area_range_c <- terra::expanse(range_c[[2]])
area_range_a <- terra::expanse(range_a[[2]])
area_range_n <- terra::expanse(range_n[[2]])
