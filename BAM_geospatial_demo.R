# Code to get condec shapefile ready for predicting occupancy
# Written by Lindsey Pekurny and Brittany Mosher

#########################
##Table of Contents    ##
#########################
#1: Load stuff - libraries, HII file, and condec file 
#2: Rasterize the condec shapefile
#3: Align two rasters
#4: Convert to occupancy


rm(list = ls())
# 
#                                                   #########################
#                                                   ##Section 1: Load stuff##
#                                                   #########################


rm(list = ls()) # start fresh

# load workspace called "post_rasterize_10Nov2020.RData"
load("C:/Users/14842/Documents/Projects/postrasterize10Nov2020.RData")
load("C:/Users/14842/Documents/Projects/7.16.21.RData")

# # set working dir
setwd("C:/Users/14842/Documents/Projects")
# 
library(maptools)
library(raster)
library(rgdal)
library(sp)
library(sf)
library(ggplot2)

 
## read in HII raster
#HIIcrop <- raster("HIIcrop2.tif") # took underscore out - did not fix error
# plot(HII_crop)
# 
# # read in deciduous/coniferous shapefile
# condec <- shapefile("condec.shp") #takes a long time
# save.image("files_read_in.RData") #saved workspace here so I won't have to reload these again


#                                                   #########################
#                                                   ##Section 2: Rasterize ##
#                                                   #########################

# summary(condec@data) #the attribute we are interested in is called 'Deciduous'
# 
# # Inspect and check loaded shapefile
# class(condec) # look at class
# head(condec) # look at data
# #str(condec) # look at structure
# crs(condec) # look at coordinate reference system - note can be saved to object & applied to other datasets
# extent(condec) # look at extent
# 
# #compare crs - these match
# crs(HII_crop)
# crs(condec)
# 
# # Rasterize
# # Check out the class of the column of interest
# class(condec$Deciduous) # returns "character"
# 
# # convert to factor
# condec$Deciduous <- as.factor(condec$Deciduous)
# nam <- unique(condec$Deciduous) # assign factor level names 
# 
# # Assign conif/decid to 1s and 2s
nam_df <- data.frame(ID = 0:(length(nam) - 1), nam = nam)
nam_df[1,1] <- 2
nam_df[2,1] <- 1
# 0 is deciduous, 1 is conferous

# 
# # Place IDs that are 1s and 2s and correspond to Coniferous Deciduous
condec$ID <- nam_df$ID[match(condec$Deciduous,nam_df$nam)] 
# 
# # Define temporary RasterLayer object which will be filled
temp.raster <- raster()
# 
# # Define raster extent to be same as original shapefile
extent(temp.raster) <- extent(condec)
# 
head(condec)
# # rasterize
condec_ras <- rasterize(x = condec, y = temp.raster, field = "ID") # rasterize by the ID field which contains 1s and 2s
save.image("post_rasterize_10Nov2020.RData") # save workspace so we don't need to re-rasterize again

# load workspace called "post_rasterize_10Nov2020.RData"
load("C:/Users/14842/Documents/Projects/postrasterize10Nov2020.RData")

#                                                   ###############################
#                                                   ##Section 3: Align 2 rasters ##
#                                                   ###############################
### From lindsey's code
HIIcrop2 <- raster("C:/Users/14842/Documents/Projects/HIIcrop2.tif")
plot(HIIcrop2)
# check out the two rasters -- what is the same and what is different?
condec_ras
HIIcrop <- HIIcrop2

# we need to crop to the same extent and resolution
res(HIIcrop)
res(condec_ras)
# they are different...HII has a smaller resolution

# StackOverflow suggests using aggregate or resample to do this
# https://gis.stackexchange.com/questions/255150/using-resample-vs-aggregate-extend-in-r-to-have-rasters-of-matching-resolutio/255155
#Condec_res <- resample(condec_ras, HII_crop, method="bilinear")
Condec_res <- resample(condec_ras, HIIcrop, method = "ngb") #I first used method bilinear, but got continuous variables out instead of 1s and 2s

#Did it work? Looks like it (we made a new condec raster is called Condec_res)
res(Condec_res) == res(HIIcrop)
extent(Condec_res) == extent(HIIcrop)
origin(Condec_res) == origin(HIIcrop)

library(RColorBrewer)
library(viridis)
viridis.pal <- viridis_pal("magma")
pal <- colorRampPalette(c("yellow","black"))

plot(Condec_res, col = pal(2))
plot(HIIcrop, col = pal(40)) ### read in tif file Lindsey   created in her code 

#                                                   ###############################
#                                                   ##Section 4: Occupancy       ##
#                                                   ###############################

# converting to occupancy. First we need some betas -- I will make them up below, but we will sub in the betas from the best model 
beta_0 <- 0.316276 # intercept term -> summary of best model, 3 values
beta_1 <- -0.458693 # effect of HII on occ
beta_2 <- 0 # effect of Conif/Decid on occ ## put in code

###--------PLOT TOP MODEL ONLY----
str(trapdateHII)
trapdateHII$beta$psi

# generate logit values using the two raster layers and the betas above
logit.raster1 <- trapdateHII$beta$psi["psi.Int","est"] +  
  trapdateHII$beta$psi["psi.psi.HII","est"] * HIIcrop2
# zip up the logits again, right now are un-crunched aka on logit scale
psi_ras_1 <- exp(logit.raster) / (1 + exp(logit.raster))
plot(psi_ras_1)

writeRaster(psi_ras_1,
            filename = "C:/Users/14842/Documents/Projects/GIS/psi_ras1.tif",
            options = c('TFW = YES'),
            overwrite = TRUE)
psi_ras_1 <- raster("C:/Users/14842/Documents/Projects/GIS/c_psiras1.tif")
plot(psi_ras_1,
     col = viridis(1e3),
     main = "Predicted Occ Model 1")

##### HIGH STANDARD ERROR FOR MODEL 1 --------------------------

h.logit.raster1 <- (trapdateHII$beta$psi["psi.Int","est"] +
                    trapdateHII$beta$psi["psi.Int","se"]) + 
                    (trapdateHII$beta$psi["psi.psi.HII","est"] +
                      trapdateHII$beta$psi["psi.psi.HII","se"]) * HIIcrop2
h.psi_ras_1 <- exp(h.logit.raster1) / (1 + exp(h.logit.raster1))
plot(h.psi_ras_1, col = viridis(1e3), main = "High Pred")

writeRaster(h.psi_ras_1, filename = "C:/Users/14842/Documents/Projects/GIS/h.psi_ras_1.tif", options = c('TFW = YES'), overwrite = TRUE)
h.psi_ras_1 <- raster("C:/Users/14842/Documents/Projects/GIS/h.psi_ras_1.tif")

#### LOW STANDARD ERROR FOR MODEL 1 --------------------------

l.logit.raster1 <- (trapdateHII$beta$psi["psi.Int","est"] -
                      trapdateHII$beta$psi["psi.Int","se"]) + #Intercept
  (trapdateHII$beta$psi["psi.psi.HII","est"] -
     trapdateHII$beta$psi["psi.psi.HII","se"]) * HIIcrop2
l.psi_ras_1 <- exp(l.logit.raster1) / (1 + exp(l.logit.raster1))
plot(l.psi_ras_1, col = viridis(1e3), main = "Low Pred")

writeRaster(l.psi_ras_1, filename = "C:/Users/14842/Documents/Projects/GIS/l.psi_ras_1.tif", options = c('TFW = YES'), overwrite = TRUE)
l.psi_ras_1 <- raster("C:/Users/14842/Documents/Projects/GIS/l.psi_ras_1.tif")

par(mfrow = c(1,2))
plot(l.psi_ras_1, col = viridis(1e3), main = "Low Pred")
plot(h.psi_ras_1, col = viridis(1e3), main = "High Pred")

### Model 2
logit.raster2 <- dateHII$beta$psi["psi.Int","est"] +  
  dateHII$beta$psi["psi.psi.HII","est"] * HIIcrop2
psi_ras_2 <- exp(logit.raster2) / (1 + exp(logit.raster2))
plot(psi_ras_2, col = viridis(1e3), main = "Model 2")

### Model 3
logit.raster3 <- total$beta$psi["psi.Int","est"] +  
  total$beta$psi["psi.psi.HII","est"] * HIIcrop2 +
  total$beta$psi["psi.psi.dc.col2", "est"] * Condec_res
psi_ras_3 <- exp(logit.raster3) / (1 + exp(logit.raster3))
plot(psi_ras_3, col = viridis(1e3), main = "Model 3, Total")

### Model 4
logit.raster4 <- dateHIIdc$beta$psi["psi.Int","est"] +  
  dateHIIdc$beta$psi["psi.psi.HII","est"] * HIIcrop2 +
  dateHIIdc$beta$psi["psi.psi.dc.col2","est"] * Condec_res
psi_ras_4 <- exp(logit.raster3) / (1 + exp(logit.raster3))
plot(psi_ras_4, col = viridis(1e3), main = "Model 4")

# Now must convert to the occupancy scale; can use expit function
library(clusterPower) # expit is in library clusterPower
occ.raster <- expit(logit.raster)
plot(occ.raster, col = viridis(1e5), main = "Predicted Occupancy")


### Plot psi_ras_2, psi_ras_3, psi_ras_4 same way as above. Give weights and multiply - make sure everything equals 1

test <- psi_ras_1 * results$table["11","wgt"] + 
        psi_ras_2 * results$table["8","wgt"] +
        psi_ras_3 * results$table["12","wgt"] +
        psi_ras_4 * results$table["9","wgt"]
writeRaster(test, filename = "C:/Users/14842/Documents/Projects/GIS/ModelAvg.tif", options = c('TFW = YES'))
### clipped in QGIS
test <- raster("C:/Users/14842/Documents/Projects/GIS/clip_modavg.tif")
plot(test, col = viridis(1e3), main = "Model Averaged")

### Make sure you make comments about what didn't work for me.

save.image("post_rasterize_10Nov2020.RData")


###------------ IUCN MAP --------------------------------------
# Brining in IUCN shapefile. Will clip in R to STATEBOUND, a shapefile I made in QGIS that I clipped raster layers to. Clipping vectors in not working in QGIS atm.
library(raster)
library(rgdal)
library(sf)
library(viridis)

IUCN <- raster::shapefile("C:/Users/14842/Documents/Projects/GIS/IUCN/iucn/data_0.shp")
states <- raster::shapefile("C:/Users/14842/Documents/Projects/GIS/stateboundaries/STATEBOUND.shp")
plot(IUCN)
plot(states)

crs(IUCN)
crs(states)
crs(sites_sf)

Crop_IUCN <- crop(IUCN, states)
extent(Crop_IUCN)
extent(test)
plot(Crop_IUCN)

par(mfrow = c(2,2))

### ------ Histogram for HII Frequency 
d <- density(freq$hii)

par(mfrow = c(2,2))
plot(d, main = "HII Value Frequency",
     ylab = "Density",
     xlab = "HII Value",
      fill = "red")
      rug(freq$hii, 
          ticksize = 0.03,
          side = 1,
          lwd = 0.5,
          col = par("fg"),
          quiet = getOption("warn") < 0)

plot(psi_ras_1, 
     col = viridis(1e3), 
     main = "Predicted Occ Model 1", 
     axes = FALSE, box = FALSE)
plot(test, col = viridis(1e3), 
     main = "Model Averaged", 
     axes = FALSE, 
     box = FALSE)
plot(IUCN, 
     col = "white", 
     main = "IUCN Range Map") ### cut out canada

plot(condec_ras,
     col = c("red","blue"),
     add = TRUE) ### Is this figure necessary?
plot(sites_sf,
     pch = 20,
     cex = 1,
     col = "black", 
     add = TRUE)

save.image("BAM.RData")

