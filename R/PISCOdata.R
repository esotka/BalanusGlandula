# ok working with the PISCO/MARINE photoplot data for Balanus 2020

piscodat<-read.table("data/photo_plots_raw_2020_0604.txt",header = TRUE)

#ok now piscodat is what remains in that .txt extracted from .xlsx
#need to sort by latitude / species / time period to get % cover 
# then to aggregate by latitude / species to get data for boxplot, mean, variance
#likely a lapply/sapply lesson...

#str(piscodat)
# 'data.frame':	44952 obs. of  16 variables:
#subset a dataframe

bal <- piscodat[piscodat$species_code=="BALGLA",]
xbar <- aggregate(bal$percent_cover,by=list(bal$latitude),mean)
names(xbar) <- c("latitude","MeanPercCover")
write.csv(xbar,"output/PISCOdata-Balanus.csv",row.names = F)

