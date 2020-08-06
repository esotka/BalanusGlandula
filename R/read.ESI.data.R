### pulls in several ESI geodatabases from 
### https://response.restoration.noaa.gov/esi_download

rm(list=ls())
library(sp)
library(sf)
library(rgdal)

### ESI dataset - central California
fgdb = "data/CentralCal_2006_GDB/CentralCaliforniaESI.gdb/"
fc = readOGR(dsn=fgdb,verbose = T,layer = 'esil_arc')
latlon <- c()
ID <- c()
for (i in 1:length(fc@lines))
  {
  tmp <- data.frame(coordinates(fc@lines[[i]]))
  latlon <- rbind(latlon,tmp)
  ID <- c(ID,rep(fc@lines[[i]]@ID,dim(tmp)[1]))
  }
out <- data.frame(lon=latlon[,1],lat=latlon[,2],ID=as.numeric(ID))
out$Landward_shoretype <- fc@data[match(out$ID,rownames(fc@data)),"Landward_shoretype"]
out$ESI <- fc@data[match(out$ID,rownames(fc@data)),"ESI"]
### remove "LINES" with estuarine codes
#out2 <- out[!(out$Landward_shoretype%in%c(
#  "10A: Salt and Brackish Water Marshes",
#  "10B: Freshwater Marshes",
#  "10D: Scrub and Shrub Wetlands",
#  "8A: Sheltered, Impermeable, Rocky Shores",
#  "8B: Sheltered, Solid Man-Made Structures",
#  "8C: Sheltered Riprap",
#  "9A: Sheltered Tidal Flats",
#  "9B: Vegetated Low Banks")),]
out2 <- out[!(out$ESI%in%c("7","8A","8B","8C","8F","9A","9B")),]
### write
write.csv(out2,file="output/centralCal.csv",quote = F,row.names = F)

### ESI dataset - Northern California
fgdb = "data/NorthernCal_2008_GDB/NorthernCaliforniaESI.gdb/"
fc = readOGR(dsn=fgdb,verbose = T,layer = 'esil')
latlon <- c()
ID <- c()
for (i in 1:length(fc@lines))
{
  tmp <- data.frame(coordinates(fc@lines[[i]]))
  latlon <- rbind(latlon,tmp)
  ID <- c(ID,rep(fc@lines[[i]]@ID,dim(tmp)[1]))
}
out <- data.frame(lon=latlon[,1],lat=latlon[,2],ID=as.numeric(ID))
out$Landward_shoretype <- fc@data[match(out$ID,rownames(fc@data)),"Landward_shoretype"]
out$ESI <- fc@data[match(out$ID,rownames(fc@data)),"ESI"]
### remove "LINES" with estuarine codes
#out2 <- out[!(out$Landward_shoretype%in%c(
#  "10A: Salt and Brackish Water Marshes",
#  "10B: Freshwater Marshes",
#  "10D: Scrub and Shrub Wetlands",
#  "8A: Sheltered, Impermeable, Rocky Shores",
#  "8B: Sheltered, Solid Man-Made Structures",
#  "8C: Sheltered Riprap",
#  "9A: Sheltered Tidal Flats",
#  "9B: Vegetated Low Banks")),]
out2 <- out[!(out$ESI%in%c("7","8A","8B","8C","8F","9A","9B")),]
### write
write.csv(out2,file="output/northernCal.csv",quote = F,row.names = F)

### ESI dataset - Southern California
fgdb = "data/SouthernCal_2010_GDB/SouthernCaliforniaESI.gdb/"
fc = readOGR(dsn=fgdb,verbose = T,layer = 'esi_arc')
latlon <- c()
ID <- c()
for (i in 1:length(fc@lines))
{
  tmp <- data.frame(coordinates(fc@lines[[i]]))
  latlon <- rbind(latlon,tmp)
  ID <- c(ID,rep(fc@lines[[i]]@ID,dim(tmp)[1]))
}
out <- data.frame(lon=latlon[,1],lat=latlon[,2],ID=as.numeric(ID))
out$Landward_shoretype <- fc@data[match(out$ID,rownames(fc@data)),"Landward_shoretype"]
out$ESI <- fc@data[match(out$ID,rownames(fc@data)),"ESI"]
### remove "LINES" with estuarine codes
#out2 <- out[!(out$Landward_shoretype%in%c(
#  "10A: Salt and Brackish Water Marshes",
#  "10B: Freshwater Marshes",
#  "10D: Scrub and Shrub Wetlands",
#  "8A: Sheltered, Impermeable, Rocky Shores",
#  "8B: Sheltered, Solid Man-Made Structures",
#  "8C: Sheltered Riprap",
#  "9A: Sheltered Tidal Flats",
#  "9B: Vegetated Low Banks")),]
out2 <- out[!(out$ESI%in%c("7","8A","8B","8C","8F","9A","9B")),]
### write
write.csv(out2,file="output/southernCal.csv",quote = F,row.names = F)






