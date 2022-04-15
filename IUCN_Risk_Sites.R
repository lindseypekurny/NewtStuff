#### Creating the IUCN + Risk + Site Map
### Pull in layers created in QGIS. Mask county map to IUCN map. Figure out Risk values for each county; export to a csv and take a gander.

setwd("C:/Users/14842/Documents/Projects")
load("C:/Users/14842/Documents/Projects/7.7.21.RData") #occ models
load("C:/Users/14842/Documents/Projects/postrasterize10Nov2020.RData") # spatial stuff from BAM

library(raster)
library(rgdal)

county <- raster::shapefile("C:/Users/14842/Documents/Projects/GIS/County/county_84.shp")
plot(county)
IUCN <- raster::shapefile("C:/Users/14842/Documents/Projects/GIS/IUCN/iucn84.shp")
plot(IUCN)

crop_county <- raster::crop(county, IUCN)
plot(crop_county)

county_risk <- raster::shapefile("C:/Users/14842/Documents/Projects/GIS/Risk/iucn_risk_clip.shp")
plot(county_risk)

df_crop_county <- as.data.frame(crop_county)
write.csv(df_crop_county, file = "C:/Users/14842/Documents/Projects/GIS/County/Crop_County.csv")
