# print coastline

rm(list=ls())
library(maps)
### habitat
cCal <- read.csv('output/centralCal.csv')
nCal <- read.csv('output/northernCal.csv')
sCal <- read.csv('output/southernCal.csv')
all <- rbind(nCal,cCal,sCal)
all$rocky <- all$ESI%in%c(
  "1A",
  "1B",
  "2A")
all$lon <- as.numeric(as.character(all$lon))
all <- all[complete.cases(all),]

### estimate cline density and proportion rocky
### slice widths in deg latitude for histogram
slice.width = .2
#how wide is the slice (e.g., ± 0.1 on each side of the midpoint yields 0.2)
# degrees latitude (69 miles is 1 deg lat) 
slice.dist <- 0.1

### clines
clines <- read.csv('output/May19.snps-CoCo.cline-parameters.csv')
lat.slices <- hist(clines$mid,col="grey",breaks = seq(30,60,slice.width)) # latitudinal slices 
out <- data.frame(mids=lat.slices$mids,clines.density=lat.slices$density)

prop.rocky <- c()
for (i in 1:dim(out)[1])
{
  targetlat <- out$mids[i] # e.g., hopkins
  targetSlice <- all[all$lat>(targetlat-slice.dist) & all$lat<(targetlat+slice.dist),]
  ### assume that the Lines are approximately equal in size. what proportion of the Lines are rocky
  targetSlice$rocky <- targetSlice$ESI%in%c(
    "1A",
    "1B",
    "2A")
  tmp <- as.data.frame(table(targetSlice$ID,targetSlice$rocky))
  tmp <- tmp[tmp$Freq>0,]
  prop.rocky <- c(prop.rocky,sum(tmp$Var2==TRUE)/dim(tmp)[1])
}

out$prop.rocky <- prop.rocky
out <- out[!(is.na(out$prop.rocky) | out$clines.density==0),]

### add percent cover from PISCO
### output from PISCOdata.R

xbar <- read.csv("output/PISCOdata-Balanus.csv")


library(maps)
pdf('output/analysis.ESI.dataMAP.pdf',width=5,height=8)
par(mar=c(1,1,1,1))
map("state","California",xlim=c(-125,-114),ylim=c(23,44),fill=T,col="lightgrey")
rockyOnly <- all[all$rocky,]
rockyOnlyNC <- rockyOnly[rockyOnly$lat>34.5,]
rockyOnlyS <- rockyOnly[rockyOnly$lat<=34.5,]
points(rockyOnly$lon,rockyOnly$lat,pch=20,cex=.1,col="red")
segments(x0=-130,y0=seq(32,43,1),x1=-114,y1=seq(32,43,1),col="black",lwd=.05)
mtext(side=2,at=32:43,paste(as.character(32:43),"º",sep=""),las=1,cex=0.75)
segments(x0=seq(-130,-110,1),y0=32,x1=seq(-130,-110,1),y1=31.5,col="black",lwd=0.05)
text(x=seq(-116,-124,-2),y=31,paste(as.character(seq(116,124,2)),"º",sep=""),cex=0.75)
rect(-125,32,-114,43)
points(x=rep(-115,dim(out)[1]),y=out$mids,pch=20,cex=out$prop.rocky*5) 
text(-115,42.5,"rocky",col="black",cex=0.5)
#points(x=rep(-119,3),y=c(42.5,41.5,40.5),cex=c(0.5,1,1.5),pch=20)
#text(x=rep(-118.5,3),y=c(42.5,41.5,40.5),c("16%","33%","50%"),pos=4)
points(x=rep(-116,dim(out)[1]),y=out$mids,pch=20,cex=out$clines.density*7,co="blue") 
text(-116,42.5,"clines",col="blue",cex=0.5)
points(x=rep(-117,dim(xbar)[1]),y=xbar$latitude,pch=20,cex=xbar$MeanPercCover/3,col="red")
text(-117,42.5,"% Cover",col="red",cex=0.5)
dev.off()
