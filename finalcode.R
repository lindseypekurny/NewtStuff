
rm(list = ls())

library(RPresence)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(readxl)

### Counted each dipnet as it's own survey event and each trap as it's own survey event
### So, each dipnet was own survey, even though it was the same day; same for traps. This allowed for multiple surveys during one visit

setwd("C:/Users/14842/Documents/Projects")

# load workspace called "post_models.RData" Do not need to re-run models
load("C:/Users/14842/Documents/Projects/post_models.RData")

#---SKIP - Loaded in post_models.RData ---------------------

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

### Read in environmental covariates
# Human Interference Index (HII)
HII <- newtdata5[,16]
hii <- HII ### Create 2 HII, one for Deciduous and one for Coniferous
str(HII)
HII <- scale(HII)
unique(HII)
HII <- as.vector(HII)

hii <- na.omit(hii)
max(hii)
min(hii)
mean(hii)
median(hii)

# Forest Type
dc <- newtdata5[,17] ### 1 is coniferous, 0 is deciduous
dc <- as.data.frame(dc)
unique(dc)

# Now I want to make a new object that replaces the forest type with deciduous or coniferous
# the pattern is if X or Y or Z is true (argument 1), then return this value (argument 2), otherwise this value (argument 3)
# I have it return the words "Deciduous" or "coniferous", but you could put numbers instead (0 and 1, or 1 and 2)
dc.col <- ifelse(dc == "Coniferous", 2, 1)
dc.col <- as.factor(dc.col)
head(dc.col)
class(dc.col)
str(dc.col)
dim(dc.col) # NULL means no longer row and column

# Now I'm going to put the two columns side-by-side to make sure it looks right
dc <- cbind(dc,dc.col)
head(dc) # looks good!
dc$dc.col
(deciduous_count <- sum(dc$dc.col == 1))
(coniferous_count <- sum(dc$dc.col == 2))

## Detection covariates

Dipnet <- newtdata5[,2:7]
Dipnet[Dipnet == "0" | Dipnet == "1"] <- 1
Dipnet <- as.data.frame(Dipnet)
glimpse(Dipnet)
unique(Dipnet)
(dv <- unlist(Dipnet))
(dv <- na.omit(dv))

Traps <- newtdata5[,8:15]
Traps[Traps == "0" | Traps == "1"] <- 2
Traps <- as.data.frame(Traps)
glimpse(Traps)
(t <- unlist(Traps))
(t <- na.omit(t))

Method = cbind(Dipnet, Traps, row.names = NULL)
Method <- as.data.frame(Method)
Method = unlist(Method) ### Method needs to look exactly like detection data before it is unlisted
Method <- as.factor(Method)
SurveyCov = as.data.frame(Method, row.names = NULL)

###-- Presence - Occupancy Models ------------------
newtpao = createPao(data = dethist5,
                  unitcov = data.frame(HII = HII, dc.col = dc.col),
                  survcov = data.frame(Method = Method,date = date),
                  title = "Covariates")

### ------------------- Occupancy Models
vddc = occMod(data = newtpao,
              type = "so",
              model = list(psi~dc.col,p~Method),
              modname = "Method and Forest",
              warn = TRUE)
vdHII = occMod(data = newtpao,
        type = "so",
        model = list(psi~HII,p~Method),
        modname = "Method and HII",
        warn = TRUE)
vdHIIdc = occMod(data = newtpao,
                 type = "so",
                 model = list(psi~HII+dc.col,p~Method),
                 modname = "Method, Forest, and HII")
# p ~ 1
condc = occMod(data = newtpao,
               type = "so",
               model = list(psi~dc.col,p~1),
               modname = "Forest")
conHII = occMod(data = newtpao,
                type = "so",
                model = list(psi~HII,p~1),
                modname = "HII")
conHIIdc = occMod(data = newtpao,
                type = "so",
                model = list(psi~HII+dc.col,p~1),
                modname = "HII and Forest")
# p ~ Date
datedc = occMod(data = newtpao,
                type = "so",
                model = list(psi~dc.col,p~date),
                modname = "Date and Forest")
dateHII = occMod(data = newtpao,
                 type = "so",
                 model = list(psi~HII,p~date),
                 modname = "Date and HII")
dateHIIdc = occMod(data = newtpao,
                   type = "so",
                   model = list(psi~HII+dc.col,p~date),
                   modname = "Date, HII, and Forest")
# p ~ Date + Method
trapdatedc = occMod(data = newtpao,
                    type = "so",
                    model = list(psi~dc.col,p~date+Method),
                    modname = "Trap, Date, and Forest")
trapdateHII = occMod(data = newtpao,
                     type = "so",
                     model = list(psi~HII,p~date+Method),
                     modname = "Trap, Date, and HII")
total = occMod(data = newtpao,
               type = "so",
               model = list(psi~HII+dc.col,p~date+Method),
               modname = "Total")
# psi ~ 1
conpsi1 = occMod(data = newtpao,
                 type = "so",
                 model = list(psi~ 1,p~date+Method),
                 modname = "Date and Method")
conpsi2 = occMod(data = newtpao,type = "so",
                 model = list(psi~ 1,p~date),
                 modname = "Date")
conpsi3 = occMod(data = newtpao,
                 type = "so",
                 model = list(psi~ 1,p~ Method),
                 modname = "Method")

# Constant
constantOCC <- occMod(data = newtpao,
                      type = "so",
                      model = list(psi~1,p~1),
                      modname = "Constant")


### AICc Table
## Name list so row names become list
mods = list(condc = condc, conHII,conHIIdc,vddc,vdHII,vdHIIdc,datedc,dateHII
            ,dateHIIdc,trapdatedc,trapdateHII,total,constantOCC, conpsi1, 
            conpsi2, conpsi3)  #conpsi1, conpsi2, conpsi3
results <- createAicTable(mods, use.aicc = TRUE) 
print(results$table)

res <- as.data.frame(results$table)
write.csv(res,
          file = "C:/Users/14842/Documents/Projects/Visuals/results.csv")
(total$beta)
(dateHIIdc$beta)
(trapdateHII$beta)

### Saving Workspace
save.image("post_models.RData") # save workspace

#---------- START RUNNING CODE AGAIN
### ---------------------------- Model Averaging

### Need same number of predictions as raster values?
new <- expand.grid(HII = HII, dc.col = dc.col)
pred.avg <- modAvg(results, param = "psi", index = 1:3, conf = 0.95, predict = TRUE, newdata = new)  ### Model avg values for top three models
pred.avg <- cbind(new, pred.avg) # combine expanded HII and dc values with predicted averaged results for each value combination
avg <- as.data.frame(pred.avg)
write.csv(avg, file = "C:/Users/14842/Documents/Projects/Visuals/pred.avg.csv")

#### Put spatial values into modAvg()
#### Double check DC values match
#### DELETE BELOW? MANUALLY MODEL AVG
new2 <- expand.grid(HII = df_HII, dc.col = df_CD) ### creates all possible combos
new2 <- data.frame(HII = df_HII, dc.col = df_CD)
head(new2)
tail(new2)

pred.avg <- modAvg(results, param = "psi", index = 1:3, conf = 0.95, predict = TRUE, newdata = new2)  ### Model avg values for top three models
pred.avg <- cbind(new2, pred.avg)
head(pred.avg)
summary(pred.avg)
max(pred.avg$est, na.rm = TRUE)

### top three models -------------------------------------------------------

model1betapsi <- as.data.frame(trapdateHII$beta$psi)
model1betap <- as.data.frame(trapdateHII$beta$p)
write.csv(model1betapsi, file = "C:/Users/14842/Documents/Projects/Visuals/model1betapsi.csv")
write.csv(model1betap, file = "C:/Users/14842/Documents/Projects/Visuals/model1betap.csv")

model2betapsi <- as.data.frame(dateHII$beta$psi)
model2betap <- as.data.frame(dateHII$beta$p)
write.csv(model2betapsi, file = "C:/Users/14842/Documents/Projects/Visuals/model2betapsi.csv")
write.csv(model2betap, file = "C:/Users/14842/Documents/Projects/Visuals/model2betap.csv")

model3betapsi <- as.data.frame(total$beta$psi)
model3betap <- as.data.frame(total$beta$p)
write.csv(model3betapsi, file = "C:/Users/14842/Documents/Projects/Visuals/model3betapsi.csv")
write.csv(model3betap, file = "C:/Users/14842/Documents/Projects/Visuals/model3betap.csv")

save.image("CompleteModels.RData")
load("CompleteModels.RData")

### Visual

### ------------------- HII Visual ----------------------------

HIIplot2 <- seq(min(HII, na.rm = T),max(HII, na.rm = T),length.out = 100)
#newdata7 = data.frame(HII = HIIplot2, dc=rep(c(0,1), each=100))
newdata7 = data.frame(HII = HIIplot2, dc = rep(c(0,1), each = 100))

HII <- HIIplot2
dc = (c(0,1))

newdata3 <- expand.grid(HII = HII)
#newdata3 <- expand.grid(HII = HII, dc = dc)
#newdata2<-newdata2[-c(601:800),]

HIIdcpred.psi = predict(total, newdata3, param = 'psi')
HIIpred.psi = predict(trapdateHII, newdata3, param = 'psi') ### Top Model
pred.HII.date = predict(dateHII, newdata3, param = 'psi')  ### Second Model

HIIdcpred.psi$forest = rep(c("deciduous","conifer"),
                           each = 100)
HIIpred.psi$HII = newdata3$HII

### Visual for HII and occupancy estimates

a <- ggplot(data = HIIpred.psi,aes(x = HII,y = est,)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = lower_0.95,
                    ymax = upper_0.95, 
                    alpha = 0.5))
b <- a + 
    xlab("Human Interference Index") + 
    ylab("Occupancy Estimate") + 
    theme_bw() + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
b + theme(legend.position = "none") 


#### HII & DC Visual ---------------------------------------

HIIdcpred.psi = predict(total, new, param = 'psi')
pred.avg$forest = rep(c("deciduous","conifer"),each = 100)
pred.avg$HII = new$HII

a <- ggplot(pred.avg,aes(x = HII,
                        y = est,
                        fill = as.factor(dc.col))) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_0.95,
                  ymax = upper_0.95, 
                  fill = as.factor(dc.col)),
                  alpha = .5)

b <- a +
  scale_fill_manual(values = alpha(c("#56B4E9","#CC79A7"),
                                   .3)) +
  labs(title = "Forest and HII Effect on N. viridescens   
       Occupancy") +
  xlab("Human Interference Index") +
  ylab("Occupancy Estimate") + 
  labs(fill = "Forest Type") +
  theme_classic()

(c <- b +
    scale_fill_manual(name = "Forest Type",
                      labels = c("Deciduous", "Coniferous"),
                      values = c("#CC79A7", "#999999"),
                      guide = guide_legend(reverse = TRUE)))


#### Date and Method  Visual---------------------------

dateplot2 <- seq(min(date, na.rm = T),
                 max(date, na.rm = T),
                 length.out = 100)
#newdata7 = data.frame(HII = HIIplot2,
                     # dc = rep(c(0,1),
                      #each = 100))

date <- dateplot2
Method = (c(0,1))

newdata3 <- expand.grid(date = date,
                        Method = Method)
pred.p = predict(trapdateHII,
                 newdata3,
                 param = 'p')

pred.p$Method = rep(c("Dipnet",
                      "Trap"),
                    each = 100)
pred.p$d.date = newdata3$date
summary(pred.p)


a <-
  ggplot(data = pred.p,
         aes(x = d.date, 
          y = est, 
          fill = as.factor(Method))) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_0.95, 
                  ymax = upper_0.95,
                  fill = as.factor(Method)), 
                  alpha = .5)


b <- a + 
  scale_fill_manual(values = alpha(c("black","gray"),
                                   .3)) + 
  labs(title = "Effect of Date and Method Type on Detection") + 
  xlab("Day of Year") + 
  ylab("Detection Estimate") + 
  labs(fill = "Method") + 
  theme_classic()

b + scale_x_continuous(breaks = c(-2,-1,0,1,2,3), 
                       labels = c("March",
                                  "April",
                                  "May",
                                  "June",
                                  "July" ,
                                  "August"))

### Mapping --------------------------------------------------------

rm(list = ls())


library(raster)
library(sp)
library(rgdal)
library(sf)
library(ggplot2)
#library(broom)
#library(mapview)
#library(leaflet)
library(plyr)
library(maptools)
library(dplyr)
library(units)
library(rgeos)

setwd("C:/Users/14842/Documents/Projects")

#HII <- raster("C:/Users/14842/Documents/Projects/GIS/HII/hii-n-america-geo-grid/hii-n-america-geo-grid/hii_n_america_grid/hii_n_amer")
#HII <- projectRaster(HII,crs = crs("+proj=longlat +datum=WGS84 +no_defs"))
#plot(HII)
#cropbox1 <- drawExtent()
#HII_Sites <- crop(HII, cropbox1)
#writeRaster(HII_Sites, filename = "HII_Sites.tif", options=c('TFW=YES'))
HII_Sites <- raster("HII_sites.tif")
plot(HII_Sites)
crs(HII_Sites)

### Standardize to HII sample
HIIcrop2 <- scale(HII_Sites)
plot(HIIcrop2)

### Write a raster layer to pull so you do not need to repeat all of the above steps - above steps are now unecessary, but I left them in case I need to edit

writeRaster(HIIcrop2, filename = "HIIcrop2.tif",
            options = c('TFW=YES'))
HIIcrop2 <- raster("C:/Users/14842/Documents/Projects/GIS/HIIcrop2.tif")
plot(HIIcrop2)

shp.r <- raster("C:/Users/14842/Documents/Projects/shp.r.tif") ### This has values taken from HII - no good
#hasValues(shp.r)
#plot(shp.r, "Deciduous")

#Create HII layer
maskHII <- mask(HIIcrop2,shp.r)
plot(maskHII)

HIIreproject <- projectRaster(HIIcrop2,
                              ras,
                              method = 'ngb')

join <- stack(ras, HIIcrop2)
plot(join)  ### Plots ID, not deciduous vs coniferous

#### Drawing DC and HII out of layers to get point values
##--------------------------- HII 

sites <- read.csv("C:/Users/14842/Documents/Projects/GIS/Coordinates_updated.csv")
sites_sf <- st_as_sf(sites,
                     coords = c("Longitude", "Latitude"),
                     crs = "+proj=longlat +datum=WGS84 +no_defs")
sitename <- sites$ï..Site

### Taking raster HII values
rasValue = extract(HII_Sites, sites_sf)
plot(rasValue)
Site.HII <- cbind(sitename, rasValue) ### Combine raster values with point and save as a CSV file.

combinePointValue = cbind(sites,rasValue)
write.csv(Site.HII,file = "Sites_HII_060121.csv")


crs(HII_Sites)
crs(sites_sf)
res(HII_Sites) == res(sites_sf)
extent(Condec_res) == extent(HII_crop)
origin(Condec_res) == origin(HII_crop)


#------------- Plotting DC & Sites -----------------#

library(RColorBrewer)
library(viridis)
c <- colorRampPalette(c("maroon","black"))


plot(Condec_res, col = pal(40))
plot(sites_sf, pch = 16,
     col = "coral2", add = TRUE)



#### ---------------- Frequency Chart
library(ggplot2)
str(HII)
count <- 1:264
freq <- cbind(count, hii)
freq <- as.data.frame(freq)
a <- ggplot(freq, aes(x = hii))
HIIfreq <- a + 
      geom_density(aes(y = ..count..), fill = "lightgray") + 
      theme_minimal()  +
      ylab("Frequency") + 
      xlab("HII Value")

HIIfreq <- HIIfreq +
  geom_rug(
    mapping = NULL,
    data = NULL,
    stat = "identity",
    position = "identity",
    outside = FALSE,
    sides = "bl",
    length = unit(0.03, "npc"),
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE
  )

HIIfreq <- HIIfreq + theme(
  # Hide panel borders and remove grid lines
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change axis line
  axis.line = element_line(colour = "black")
)



str(dc)
count <- 1:264
dc.col

save.image("finalcode.RData")


### Updating R

#install.packages('installr')
#require(installr)
#updateR()



