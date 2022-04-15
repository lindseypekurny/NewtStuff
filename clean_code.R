rm(list = ls())

library(RPresence)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(readxl)

### Counted each dipnet as it's own survey event and each trap as it's own survey event
### So, each dipnet was own survey, even though it was the same day; same for traps. Attempting to 
#increase survey visits.

setwd("C:/Users/14842/Documents/Projects")

# load workspace called "post_models.RData" Do not need to re-run models
load("C:/Users/14842/Documents/Projects/post_models.RData")

#---------------- SKIP ------------------------

## Read in csv containing data
newtdata5 <- read.csv("data_5.24.21.csv", skip = 1, na.strings = "-")

## Detection history
dethist5 <- newtdata5[,2:15]
print(head(dethist5))

## Read in data columns
Site <- newtdata5[,1]

## Read in survey dates and turn them into Julian date, aka day of year
date <- newtdata5[,21:34] 
date <- unlist(date)
d.date <- as.Date(date, format = "%m/%d/%Y")
doy.date <- strftime(d.date, format = "%j")
date <- as.numeric(doy.date)
date <- scale(date) # standardize

### Read in environmental covariates ----------------------------------
# HII - Human Interference Index
### Create 2 HII, one for Deciduous forest and one for Coniferous
HII <- newtdata5[,16]
hii <- HII 
str(HII)
HII <- scale(HII)
unique(HII)
HII <- as.vector(HII)
hii <- na.omit(hii)

# Forest Type
### 1 is coniferous, 0 is deciduous
dc <- newtdata5[,17]
dc <- as.data.frame(dc)

# Now I want to make a new object that replaces the forest type with deciduous or coniferous. The pattern is if X or Y or Z is true (argument 1), then return this value (argument 2), otherwise this value (argument 3).I have it return the words "Deciduous" or "coniferous".
dc.col <- ifelse(dc == "Coniferous", 2, 1)
dc.col <- as.factor(dc.col)

# Now I'm going to put the two columns side-by-side to make sure it looks right
dc <- cbind(dc,dc.col)
head(dc)

## Detection covariates
# Dipnet also includes VES surveys
Dipnet <- newtdata5[,2:7]
Dipnet[Dipnet == "0" | Dipnet == "1"] <- 1
Dipnet <- as.data.frame(Dipnet)
(dv <- unlist(Dipnet))
(dv <- na.omit(dv))

Traps <- newtdata5[,8:15]
Traps[Traps == "0" | Traps == "1"] <- 2
Traps <- as.data.frame(Traps)
(t <- unlist(Traps))
(t <- na.omit(t))

Method = cbind(Dipnet,
               Traps,
               row.names = NULL)
Method <- as.data.frame(Method)
Method = unlist(Method) ### Method needs to look exactly like detection data before it is unlisted
Method <- as.factor(Method)
SurveyCov = as.data.frame(Method, row.names = NULL)

###-- Presence - Occupancy Models ------------------
newtpao = createPao(data = dethist5,
                    unitcov = data.frame(HII = HII,
                                         dc.col =  dc.col), 
                    survcov = data.frame(Method = Method,
                                         date = date),
                    title = "Covariates")

### -------- Occupancy Models
vddc = occMod(data = newtpao,
              type = "so",
              model = list(psi ~ dc.col,
                           p ~ Method),
              modname = "Method and Forest",
              warn = TRUE)
vdHII = occMod(data = newtpao,
               type = "so",
               model = list(psi ~ HII,
                            p ~ Method),
               modname = "Method and HII",
               warn = TRUE)
vdHIIdc = occMod(data = newtpao,
                 type = "so",
                 model = list(psi ~ HII + dc.col,
                              p ~ Method),
                 modname = "Method, Forest, and HII")
# p ~ 1
condc = occMod(data = newtpao,
               type = "so",
               model = list(psi ~ dc.col,
                            p ~ 1),
               modname = "Forest")
conHII = occMod(data = newtpao,
                type = "so",
                model = list(psi ~ HII,
                             p ~ 1),
                modname = "HII")
conHIIdc = occMod(data = newtpao,
                  type = "so",
                  model = list(psi ~ HII + dc.col,
                               p ~ 1),
                  modname = "HII and Forest")
# p ~ Date
datedc = occMod(data = newtpao,
                type = "so",
                model = list(psi ~ dc.col,
                             p ~ date),
                modname = "Date and Forest")
dateHII = occMod(data = newtpao,
                 type = "so",
                 model = list(psi ~ HII,
                              p ~ date),
                 modname = "Date and HII")
dateHIIdc = occMod(data = newtpao,
                   type = "so",
                   model = list(psi ~ HII + dc.col,
                                p ~ date),
                   modname = "Date, HII, and Forest")
# p ~ Date + Method
trapdatedc = occMod(data = newtpao,
                    type = "so",
                    model = list(psi ~ dc.col,
                                 p ~ date + Method),
                    modname = "Method, Date, and Forest")
trapdateHII = occMod(data = newtpao,
                     type = "so",
                     model = list(psi ~ HII,
                                  p ~ date + Method),
                     modname = "Method, Date, and HII")
total = occMod(data = newtpao,
               type = "so",
               model = list(psi ~ HII + dc.col,
                            p ~ date + Method),
               modname = "Total")
# psi ~ 1
conpsi1 = occMod(data = newtpao,
                 type = "so",
                 model = list(psi~ 1,
                              p ~ date + Method),
                 modname = "Date and Method")
conpsi2 = occMod(data = newtpao,type = "so",
                 model = list(psi ~ 1,
                              p ~ date),
                 modname = "Date")
conpsi3 = occMod(data = newtpao,
                 type = "so",
                 model = list(psi ~ 1,
                              p ~ Method),
                 modname = "Method")

# Constant
constantOCC <- occMod(data = newtpao,
                      type = "so",
                      model = list(psi ~ 1,
                                   p ~ 1),
                      modname = "Constant")


### AICc Table
mods = list(condc = condc, 
            conHII,
            conHIIdc,
            vddc,
            vdHII,
            vdHIIdc,
            datedc,
            dateHII,
            dateHIIdc,
            trapdatedc,
            trapdateHII,
            total,
            constantOCC,
            conpsi1, 
            conpsi2,
            conpsi3)
results <- createAicTable(mods, 
                          use.aicc = TRUE) 
print(results$table)
res <- as.data.frame(results$table)

#### Plotting Models & Model Averaging 
library(raster)
library(viridis)

# rasters to read in:
## HII raster cropped to Northeast - still need to crop it to state level
HIIcrop2 <- raster("C:/Users/14842/Documents/Projects/HIIcrop2.tif")
plot(HIIcrop2)

# Top Model -- HII Only
# generate logit values using the two raster layers and the betas above
logit.raster1 <- trapdateHII$beta$psi["psi.Int","est"] + 
                  trapdateHII$beta$psi["psi.psi.HII","est"] *
  HIIcrop2 

# zip up the logits again, right now are un-crunched aka on logit scale
psi_ras_1 <- exp(logit.raster1) /
  (1 + exp(logit.raster1))
plot(psi_ras_1)

writeRaster(psi_ras_1,
            filename = "C:/Users/14842/Documents/Projects/GIS/psi_ras1.tif",
            options = c('TFW = YES'),
            overwrite = TRUE)
psi_ras_1 <- raster("C:/Users/14842/Documents/Projects/GIS/c_psiras1.tif")
plot(psi_ras_1,
     col = viridis(1e3),
     main = "Predicted Occ Model 1",
     )

### Model 2 --- HII and Date 
logit.raster2 <- dateHII$beta$psi["psi.Int","est"] +  
  dateHII$beta$psi["psi.psi.HII","est"] * 
  HIIcrop2
psi_ras_2 <- exp(logit.raster2) /
  (1 + exp(logit.raster2))
plot(psi_ras_2,
     col = viridis(1e3),
     main = "Model 2")

### Model 3 --- HII & Forest Type & Date & Method
logit.raster3 <- total$beta$psi["psi.Int","est"] +  
  total$beta$psi["psi.psi.HII","est"] * 
  HIIcrop2 +
  total$beta$psi["psi.psi.dc.col2", "est"] *
  Condec_res

psi_ras_3 <- exp(logit.raster3) /
  (1 + exp(logit.raster3))
plot(psi_ras_3,
     col = viridis(1e3),
     main = "Model 3, Total")

### Model 4 --- HII & Forest Type & Date
logit.raster4 <- dateHIIdc$beta$psi["psi.Int","est"] +  
  dateHIIdc$beta$psi["psi.psi.HII","est"] *
  HIIcrop2 +
  dateHIIdc$beta$psi["psi.psi.dc.col2","est"] *
  Condec_res
psi_ras_4 <- exp(logit.raster3) /
  (1 + exp(logit.raster3))
plot(psi_ras_4,
     col = viridis(1e3),
     main = "Model 4")

# Now must convert to the occupancy scale; can use expit function
library(clusterPower) # expit is in library clusterPower
occ.raster <- expit(logit.raster)
plot(occ.raster,
     col = viridis(1e5),
     main = "Predicted Occupancy")

### Model Averaged
### Plot psi_ras_2, psi_ras_3, psi_ras_4 same way as above. Give weights and multiply - make sure everything equals 1

test <- psi_ras_1 * results$table["11","wgt"] + 
  psi_ras_2 * results$table["8","wgt"] +
  psi_ras_3 * results$table["12","wgt"] +
  psi_ras_4 * results$table["9","wgt"]
writeRaster(test, filename = "C:/Users/14842/Documents/Projects/GIS/ModelAvg.tif", options = c('TFW = YES'))
### clipped in QGIS
test <- raster("C:/Users/14842/Documents/Projects/GIS/clip_modavg.tif")
plot(test, col = viridis(1e3), main = "Model Averaged")

