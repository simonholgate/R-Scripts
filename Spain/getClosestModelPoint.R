# Scripts to find the closest model points to Spanish stations
# so we can extract those points from the POLCOMS model data

library(spam)

# Spanish station lats and lons
load("~/diskx/polcoms/spain/spanishLatLons.RData")

lenlatLons <- length(latLons[,1])
stnLats <- latLons[,1]
stnLons <- latLons[,2]

for (i in 1:lenlatLons){
  if (stnLons[i] > 180) {
    stnLons[i] <- stnLons[i] - 360
  }
}

# POLCOMS grid data
load("~/diskx/polcoms/iseajseanpsea.Rdata")

lat<-seq(from=40.11111,by=1/9,length=224)
lon<-seq(from=-19.83333,by=1/6,length=198)

minlat <- min(lat)
minlon <- min(lon)
maxlat <- max(lat)
maxlon <- max(lon)

l<-198
m<-224
lm<-l*m
n<-40

# Now find nearest grid point to each 
x <- cbind(stnLons, stnLats)

y <- expand.grid(lon,lat)

#library(spam)
#source("~/workspace/RScripts/greatcircle.R")
near <- nearest.dist(x,y)
#near <- array(NA, dim=c(length(lon), length(lat), length(latLons[,2])))
#for (i in 1:length(lon)){
#  for (j in 1:length(lat)){
#    junk <- vector(length=length(latLons[,2]), mode="numeric")
#    for (k in 1:length(latLons[,2])){
#      junk[k] <- greatcircle(x[k,2], x[k,1], lon[i], lat[j])
#    }
#  }
#}
#nearest <- which(near==min(near),arr.ind=T)

gc()

nearest <- vector(length=lenlatLons, mode="integer")
for (i in 1:lenlatLons){
  gtzero <- which(near[i,]>0)
# There are some altimetry points outside of the model region so these will
# return all zeros
  if (length(gtzero>0)){
    nearest[i] <- gtzero[which.min(near[i,gtzero])]
  } else {
    nearest[i] <- NA
  }
}
# Now extract these grid points from the model
spanishStns <- y[nearest,]

# Change the lats and lons into x and y
spanishStnsXY <- array(NA, dim=c(lenlatLons,2))
for (i in 1:lenlatLons){
  junk <- which(lon==spanishStns[i,1])
  if(length(junk)>0){
    spanishStnsXY[i,1] <- junk
  }
  
  junk <- which(lat==spanishStns[i,2])
  if(length(junk)>0){
    spanishStnsXY[i,2] <- junk
  }
}

# Stations 11 & 12 (Vigo I & II) are NAs in this method so replace empirically
spanishStnsXY[11:12,1] <- 67

save(file="spanishStnsLatLon.RData", spanishStns, nearest, spanishStnsXY)

