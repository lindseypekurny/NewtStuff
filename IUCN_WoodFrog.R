library(raster)

county <- raster::shapefile("C:/Users/14842/Documents/Projects/GIS/County/county_84.shp")
plot(county)
crs(county)

IUCN <- raster::shapefile("C:/Users/14842/Documents/Projects/GIS/woodfrog/data_0.shp")
plot(IUCN, add = TRUE, col = "blue", main = "IUCN Range Map for Wood Frog")
crs(IUCN)
