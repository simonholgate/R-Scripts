# Scripts to extract the altimetry tracks 61, 70, 239 and 248 provided by
# Laurent Testut from the POLCOMS model data

# Altimetry track lats and lons
load("~/diskx/altimetry/netcdf_brest_newlin/trackLatLon.RData")
lenlatvec <- length(lat_vec)
lenlatvec61 <- length(lat_vec61)
lenlatvec70 <- length(lat_vec70)
lenlatvec239 <- length(lat_vec239)
lenlatvec248 <- length(lat_vec248)

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
x <- cbind(lon_vec, lat_vec)
x61 <- cbind(lon_vec61, lat_vec61)
x70 <- cbind(lon_vec70, lat_vec70)
x239 <- cbind(lon_vec239, lat_vec239)
x248 <- cbind(lon_vec248, lat_vec248)
y <- expand.grid(lon,lat)

#library(spam)
#source("~/workspace/RScripts/greatcircle.R")
near <- nearest.dist(x,y)
#near <- array(NA, dim=c(length(lon), length(lat), length(lon_vec)))
#for (i in 1:length(lon)){
#  for (j in 1:length(lat)){
#    junk <- vector(length=length(lon_vec), mode="numeric")
#    for (k in 1:length(lon_vec)){
#      junk[k] <- greatcircle(x[k,2], x[k,1], lon[i], lat[j])
#    }
#  }
#}
#nearest <- which(near==min(near),arr.ind=T)


near61 <- nearest.dist(x61,y)
near70 <- nearest.dist(x70,y)
near239 <- nearest.dist(x239,y)
near248 <- nearest.dist(x248,y)
gc()

nearest <- vector(length=lenlatvec, mode="integer")
for (i in 1:lenlatvec){
  gtzero <- which(near[i,]>0)
# There are some altimetry pointsoutside of the model region so these will
# return all zeros
  if (length(gtzero>0)){
    nearest[i] <- gtzero[which.min(near[i,gtzero])]
  } else {
    nearest[i] <- NA
  }
}
# Now extract these grid points from the model
altimTrackPoints <- y[nearest,]

nearest61 <- vector(length=lenlatvec61, mode="integer")
for (i in 1:lenlatvec61){
  gtzero <- which(near61[i,]>0)
# There are some altimetry pointsoutside of the model region so these will
# return all zeros
  if (length(gtzero>0)){
    nearest61[i] <- gtzero[which.min(near61[i,gtzero])]
  } else {
    nearest61[i] <- NA
  }
}
# Now extract these grid points from the model
altimTrackPoints61 <- y[nearest61,]

nearest70 <- vector(length=lenlatvec70, mode="integer")
for (i in 1:lenlatvec70){
  gtzero <- which(near70[i,]>0)
# There are some altimetry pointsoutside of the model region so these will
# return all zeros
  if (length(gtzero>0)){
    nearest70[i] <- gtzero[which.min(near70[i,gtzero])]
  } else {
    nearest70[i] <- NA
  }
}
# Now extract these grid points from the model
altimTrackPoints70 <- y[nearest70,]

nearest239 <- vector(length=lenlatvec239, mode="integer")
for (i in 1:lenlatvec239){
  gtzero <- which(near239[i,]>0)
# There are some altimetry pointsoutside of the model region so these will
# return all zeros
  if (length(gtzero>0)){
    nearest239[i] <- gtzero[which.min(near239[i,gtzero])]
  } else {
    nearest239[i] <- NA
  }
}
# Now extract these grid points from the model
altimTrackPoints239 <- y[nearest239,]

nearest248 <- vector(length=lenlatvec248, mode="integer")
for (i in 1:lenlatvec248){
  gtzero <- which(near248[i,]>0)
# There are some altimetry pointsoutside of the model region so these will
# return all zeros
  if (length(gtzero>0)){
    nearest248[i] <- gtzero[which.min(near248[i,gtzero])]
  } else {
    nearest248[i] <- NA
  }
}
# Now extract these grid points from the model
altimTrackPoints248 <- y[nearest248,]

library(fields)
x11()
par(family='HersheySans')
plot(altimTrackPoints)
world(add=T)

x11()
par(family='HersheySans')
plot(altimTrackPoints61, col='blue')
points(altimTrackPoints70, col='red')
points(altimTrackPoints239, col='cyan')
points(altimTrackPoints248, col='magenta')
world(add=T)

save(file="altimTrackPoints.RData", altimTrackPoints, nearest, 
   altimTrackPoints61, nearest61, altimTrackPoints70, nearest70, 
   altimTrackPoints239, nearest239, altimTrackPoints248, nearest248)

