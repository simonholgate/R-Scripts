##****************************************************************************************************##
## Calculate EOFs over the N Atlantic using ACRE wind stress data to give some feel for the data
## Just focus on N Atlantic 100W-15E and 5S-80N
##****************************************************************************************************##

##****************************************************************************************************##
########################################################################################################
## Functions for use below                                                                            ##
########################################################################################################
##****************************************************************************************************##

## monthly.ws.ncdf <- function(path, nlon, nlat, nyr){
##   ## Read the sea level data into an array
##   nmon <- 12

##   nc <- open.ncdf(path)
##   slp.id <- nc$var[[2]]
##   slpMonthlyArray <- get.var.ncdf(nc, slp.id)
##   close.ncdf(nc)

##   dim(slpMonthlyArray) <- c(nlon,nlat,nmon,nyr)
##   slpMonthlyArray
## }

##*******************************************************************************************************

## natl.ws <- function(slpArray, lon, lat){
##   ## Extract the N Atl region from the SLP array
  
##   nlon <- dim(slpArray)[1]
##   nlat <- dim(slpArray)[2]
##   nyr <- dim(slpArray)[4]
##   nmon <- 12
  

##   nNAtlLon <- length(lon)
##   nNAtlLat <- length(lat)
  
##   dim(slpArray)<-c(nlon,nlat,nmon*nyr)
  
##   slpNAtlantic <- slpArray[lon,lat,]
##   dim(slpNAtlantic) <- c(nNAtlLon,nNAtlLat,nmon,nyr)
  
##   slpNAtlantic
## }
##*****************************************************************************

natl.ws <- function(dataArray, lon, lat){
  ## Extract the N Atl region from the annual mean  wind stress array
  dataNAtlantic <- dataArray[lon,lat,]
  dataNAtlantic
}

##*****************************************************************************

natl.lon.lat.ws <- function(xlon, ylat){

  dlon <- xlon[2]-xlon[1]
  dlat <- ylat[2]-ylat[1]

  # Atlantic covers 100W to 15E and 5S to 80N
  # We need to see whether we are working in 0 to 360 or -180 to 180
  # co-ordinate system
  if(min(xlon)<0){
    nAtlLon <- match(seq(from=-101.250,to=15,by=dlon), xlon)
  } else {
    nAtlLon <- c(match(seq(from=261.250,to=359,by=dlon), xlon),
                 match(seq(from=0,to=15,by=dlon), xlon))
  }
  # windstress lats aren't evenly spaced so just write this in
    nAtlLat <- c(5:50)
  
  list(lon=nAtlLon,lat=nAtlLat)
}

## ##*******************************************************************************************************

## annual.slp <- function(slpArray){
##   ## Convert a 4D array of monthly data into a 3D annual array
  
##   nlon <- dim(slpArray)[1]
##   nlat <- dim(slpArray)[2]
##   nyr <- dim(slpArray)[4]
  
##   slpAnnualArray <- array(NA, dim=c(nlon,nlat,nyr))
##   for(i in 1:nlon){
##     for(j in 1:nlat){
##       slpAnnualArray[i,j,] <- colMeans(slpArray[i,j,,])
##     }
##   }
  
##   slpAnnualArray
## }

## ##*******************************************************************************************************

## natl.lon.lat <- function(xlon, ylat){

##   dlon <- xlon[2]-xlon[1]
##   dlat <- ylat[2]-ylat[1]

##   # Atlantic covers 100W to 15E and 5S to 80N
##   # We need to see whether we are working in 0 to 360 or -180 to 180 co-ordinate system
##   if(min(xlon)<0){
##     nAtlLon <- match(seq(from=-100,to=15,by=dlon), xlon)
##   } else {
##     nAtlLon <- c(match(seq(from=260,to=359,by=dlon), xlon), match(seq(from=0,to=15,by=dlon), xlon))
##   }
##   # Do we go from N to S or S to N?
##   if(dlat>0){ # S to N
##     nAtlLat <- match(seq(from=-5,to=80, by=dlat), ylat)
##   } else { # N to S
##     nAtlLat <- match(seq(from=80,to=-5, by=dlat), ylat)
##   }
  
##   list(lon=nAtlLon,lat=nAtlLat)
## }

##*******************************************************************************************************

#########################
## Non-functional part ##
#########################

library(fields)
library(maps)
library(mapdata)
library(ncdf)
library(robust)
source("/home/simonh/workspace/RScripts/matrixMethods.R")

nlon<-192
nlat<-94
nyr<-138

## lon<-seq(from=0,by=2,length=nlon)
## lon2 <- c(seq(from=0,by=2,to=180),seq(from=-178,by=2,to=-2))
lon<-seq(from=0,by=1.875,length=nlon)
lon2 <- c(seq(from=0,by=1.875,to=180),seq(from=-178.125,by=1.875,to=-1.875))

## lat<-seq(from=90,by=-2,length=nlat)     
wslon <- seq(from=-180, to=178.125, by=1.875)
lat <- c(88.542, 86.6531, 84.7532, 82.8508, 80.9473, 79.0435, 77.1394, 75.2351, 
    73.3307, 71.4262, 69.5217, 67.6171, 65.7125, 63.8079, 61.9033, 59.9986, 
    58.0939, 56.1893, 54.2846, 52.3799, 50.4752, 48.5705, 46.6658, 44.7611, 
    42.8564, 40.9517, 39.047, 37.1422, 35.2375, 33.3328, 31.4281, 29.5234, 
    27.6186, 25.7139, 23.8092, 21.9044, 19.9997, 18.095, 16.1902, 14.2855, 
    12.3808, 10.47604, 8.57131, 6.66657, 4.76184, 2.8571, 0.952368, 
    -0.952368, -2.8571, -4.76184, -6.66657, -8.57131, -10.47604, -12.3808, 
    -14.2855, -16.1902, -18.095, -19.9997, -21.9044, -23.8092, -25.7139, 
    -27.6186, -29.5234, -31.4281, -33.3328, -35.2375, -37.1422, -39.047, 
    -40.9517, -42.8564, -44.7611, -46.6658, -48.5705, -50.4752, -52.3799, 
    -54.2846, -56.1893, -58.0939, -59.9986, -61.9033, -63.8079, -65.7125, 
    -67.6171, -69.5217, -71.4262, -73.3307, -75.2351, -77.1394, -79.0435, 
    -80.9473, -82.8508, -84.7532, -86.6531, -88.542)

ws.yrs <- c(1871:2008)

newlyn.start.year.partial <- 1953
newlyn.end.year.partial <- 2008
newlyn.start.year.pred <- 1916
newlyn.end.year.pred <- 1943
newlyn.start.year.full <- 1916
newlyn.end.year.full <- 2008

##*****************************************************************************
## N Atlantic wind stress

load("~/data/ACRE/wsACRE.RData")

nAtlLonLatIndex <- natl.lon.lat.ws(wslon,flip.matrix(lat))
wsNAtlArrayE <- natl.ws(wsEAnnualMean, nAtlLonLatIndex$lon,
                         nAtlLonLatIndex$lat)
wsNAtlArrayN <- natl.ws(wsNAnnualMean, nAtlLonLatIndex$lon,
                         nAtlLonLatIndex$lat)

## Mean wind stress
lenLon <- length(nAtlLonLatIndex$lon)
lenLat <- length(nAtlLonLatIndex$lat)
mwsNAtlArrayE <- array(NA, dim=c(lenLon, lenLat))
mwsNAtlArrayN <- array(NA, dim=c(lenLon, lenLat))
for(i in 1:lenLon){
  for(j in 1:lenLat){
    mwsNAtlArrayE[i,j] <- mean(wsNAtlArrayE[i,j,], na.rm=T)
    mwsNAtlArrayN[i,j] <- mean(wsNAtlArrayN[i,j,], na.rm=T)
  }
}

## ##*******************************************************************************************************
## ## N Atlantic pressure

## slpMonthlyArray <- monthly.slp.ncdf("~/data/ACRE/prmsl.mon.mean.nc", nlon, nlat, nyr)
## # Convert Pa to Mb
## slpMonthlyArray <- slpMonthlyArray/100

##nAtlLonLatIndex <- natl.lon.lat.ws(lon,lat)
nAtlLon <- wslon[nAtlLonLatIndex$lon]
nAtlLat <- lat[nAtlLonLatIndex$lat]

## lnAtlLon <- length(nAtlLon)
## lnAtlLat <- length(nAtlLat)
nstns <- lenLon*lenLat

## EwsNAtlAnnArray <- natl.ws(wsNatlArray, nAtlLonLatIndex$lon, nAtlLonLatIndex$lat)
##  <- annual.(slpNAtlArray)

## ##*******************************************************************************************************

## Prepare annual data for EOF caculations
dmEwsAnnual <- wsNAtlArrayE
dmNwsAnnual <- wsNAtlArrayN

## Reshape to 2D array
dim(dmEwsAnnual) <- c(nstns,nyr)
dim(dmNwsAnnual) <- c(nstns,nyr)

## Rotate so that we have nstns columns
dmEwsAnnual <- t(dmEwsAnnual)
dmNwsAnnual <- t(dmNwsAnnual)

## Remove column means
cEmeans <- colMeans(dmEwsAnnual, na.rm=T)
cNmeans <- colMeans(dmNwsAnnual, na.rm=T)
for(i in 1:nstns){
  dmEwsAnnual[,i] <- dmEwsAnnual[,i]-cEmeans[i]
  dmNwsAnnual[,i] <- dmNwsAnnual[,i]-cNmeans[i]
}


#######################################################################################################
## SVD method                                                                                        ##
#######################################################################################################

## Spatial loading pattern
svdDmEwsAnnual <- svd(dmEwsAnnual)
svdDmNwsAnnual <- svd(dmNwsAnnual)

## Amount of variance explained is in d^2
svdVarExpE <- svdDmEwsAnnual$d^2/sum(svdDmEwsAnnual$d^2)*100
svdVarExpN <- svdDmNwsAnnual$d^2/sum(svdDmNwsAnnual$d^2)*100

## PCA1
PCAann1E <- svdDmEwsAnnual$v[,1]
dim(PCAann1E) <- c(lenLon, lenLat)
PCAann1N <- svdDmNwsAnnual$v[,1]
dim(PCAann1N) <- c(lenLon, lenLat)

## PCA2
PCAann2E <- svdDmEwsAnnual$v[,2]
dim(PCAann2E) <- c(lenLon, lenLat)
PCAann2N <- svdDmNwsAnnual$v[,2]
dim(PCAann2N) <- c(lenLon, lenLat)

## PCA3
PCAann3E <- svdDmEwsAnnual$v[,3]
dim(PCAann3E) <- c(lenLon, lenLat)
PCAann3N <- svdDmNwsAnnual$v[,3]
dim(PCAann3N) <- c(lenLon, lenLat)

## Find expansion coefficients
ECann1E <- dmEwsAnnual %*% svdDmEwsAnnual$v[,1]
ECann2E <- dmEwsAnnual %*% svdDmEwsAnnual$v[,2]
ECann3E <- dmEwsAnnual %*% svdDmEwsAnnual$v[,3]
ECann1N <- dmNwsAnnual %*% svdDmNwsAnnual$v[,1]
ECann2N <- dmNwsAnnual %*% svdDmNwsAnnual$v[,2]
ECann3N <- dmNwsAnnual %*% svdDmNwsAnnual$v[,3]


#######################################################################################################
## Eigenvector method                                                                                ##
#######################################################################################################

## Form the co-variance matrix
covDmEwsAnnual <- t(dmEwsAnnual) %*% dmEwsAnnual
covDmNwsAnnual <- t(dmNwsAnnual) %*% dmNwsAnnual

## Calculate eigenvalues and eigenvectors
eigCovDmEwsAnnual <- eigen(covDmEwsAnnual)
eigCovDmNwsAnnual <- eigen(covDmNwsAnnual)

## % variance explained
eigVarExpE <- eigCovDmEwsAnnual$values/sum(eigCovDmEwsAnnual$values)*100
eigVarExpN <- eigCovDmNwsAnnual$values/sum(eigCovDmNwsAnnual$values)*100

## EOF1
EOFann1E <- eigCovDmEwsAnnual$vectors[,1]
dim(EOFann1E) <- c(lenLon, lenLat)
EOFann1N <- eigCovDmNwsAnnual$vectors[,1]
dim(EOFann1N) <- c(lenLon, lenLat)

## EOF2
EOFann2E <- eigCovDmEwsAnnual$vectors[,2]
dim(EOFann2E) <- c(lenLon, lenLat)
EOFann2N <- eigCovDmNwsAnnual$vectors[,2]
dim(EOFann2N) <- c(lenLon, lenLat)

## EOF3
EOFann3E <- eigCovDmEwsAnnual$vectors[,3]
dim(EOFann3E) <- c(lenLon, lenLat)
EOFann3N <- eigCovDmNwsAnnual$vectors[,3]
dim(EOFann3N) <- c(lenLon, lenLat)

## Find expansion coefficients
CannE <- diag(x=eigCovDmEwsAnnual$values,nrow=nstns, ncol=nstns)
CannN <- diag(x=eigCovDmNwsAnnual$values,nrow=nstns, ncol=nstns)

PCann1E <- dmEwsAnnual %*% CannE[,1]
PCann2E <- dmEwsAnnual %*% CannE[,2]
PCann3E <- dmEwsAnnual %*% CannE[,3]
PCann1N <- dmNwsAnnual %*% CannN[,1]
PCann2N <- dmNwsAnnual %*% CannN[,2]
PCann3N <- dmNwsAnnual %*% CannN[,3]

## PCann1 looks very odd - see plots below. Percentage variance explained and spatial patterns are
## identical. However, prefer SVD method on basis o time series reconstruction.

##***************************************************************************************************##
##***************************************************************************************************##
## Plots                                                                                             ##
##***************************************************************************************************##
##***************************************************************************************************##

#######################################################################################################
## Basic stress map                                                                                ##
#######################################################################################################
## x11()
## filled.contour(nAtlLon,flip.matrix(nAtlLat),mirror.matrix(wsNAtlArrayE[,,1]),
##                color=tim.colors, nlevels=20, plot.axes = { world(add=T) }, zlim=c(-0.2, 0.2))
x11()
filled.contour(nAtlLon,flip.matrix(nAtlLat),wsNAtlArrayE[,,1], main="Wind Stress E, 1871",
               color=tim.colors, nlevels=20, plot.axes = { world(add=T) }, zlim=c(-0.15, 0.15))
x11()
filled.contour(nAtlLon,flip.matrix(nAtlLat),wsNAtlArrayN[,,1], main="Wind Stress N, 1871",
               color=tim.colors, nlevels=20, plot.axes = { world(add=T) }, zlim=c(-0.15, 0.15))

#######################################################################################################
## SVD method                                                                                        ##
#######################################################################################################

## Map
## x11()
## ## PCA1
## contour(nAtlLon,flip.matrix(nAtlLat), mirror.matrix(PCAann1E), lwd=2, col='red', main="Wind stress PCAs E")
## ## PCA2
## contour(nAtlLon,flip.matrix(nAtlLat), mirror.matrix(PCAann2E), lty='dotdash', add=T, col='blue')
## ## PCA3
## contour(nAtlLon,flip.matrix(nAtlLat), mirror.matrix(PCAann3E), lty='dotted', add=T, col='green')
## world(add=T)
x11()
## PCA1
contour(nAtlLon,flip.matrix(nAtlLat), PCAann1E, lwd=2, col='red', main="Wind stress PCAs E")
## PCA2
contour(nAtlLon,flip.matrix(nAtlLat), PCAann2E, lty='dotdash', add=T, col='blue')
## PCA3
contour(nAtlLon,flip.matrix(nAtlLat), PCAann3E, lty='dotted', add=T, col='green')
world(add=T)

## x11()
## ## PCA1
## contour(nAtlLon,flip.matrix(nAtlLat), mirror.matrix(PCAann1N), lwd=2, col='red', main="Wind stress PCAs N")
## ## PCA2
## contour(nAtlLon,flip.matrix(nAtlLat), mirror.matrix(PCAann2N), lty='dotdash', add=T, col='blue')
## ## PCA3
## contour(nAtlLon,flip.matrix(nAtlLat), mirror.matrix(PCAann3N), lty='dotted', add=T, col='green')
## world(add=T)
x11()
## PCA1
contour(nAtlLon,flip.matrix(nAtlLat), PCAann1N, lwd=2, col='red', main="Wind stress PCAs N")
## PCA2
contour(nAtlLon,flip.matrix(nAtlLat), PCAann2N, lty='dotdash', add=T, col='blue')
## PCA3
contour(nAtlLon,flip.matrix(nAtlLat), PCAann3N, lty='dotted', add=T, col='green')
world(add=T)

## Time series pattern
#quartz()
x11()
plot(ws.yrs,ECann1E, type='l', col='red', main="SVD method time series pattern, E")
lines(ws.yrs,ECann2E, col='blue')
lines(ws.yrs,ECann3E, col='magenta')
x11()
plot(ws.yrs,ECann1N, type='l', col='red', main="SVD method time series pattern, N")
lines(ws.yrs,ECann2N, col='blue')
lines(ws.yrs,ECann3N, col='magenta')


#######################################################################################################
## Eigenvector method                                                                                ##
#######################################################################################################

## Map 
#quartz()
x11()
#PCA1
# png("pca1NAtl.png")
contour(nAtlLon,flip.matrix(nAtlLat), EOFann1E, lwd=2, col='red', labcex=1, main="Wind stress EOFs E")
# world(add=T)
# dev.off()
##PCA2
# png("pca2NAtl.png")
contour(nAtlLon,flip.matrix(nAtlLat), EOFann2E, lty='dotdash', col='blue', lwd=2, labcex=1, add=T)
# world(add=T)
# dev.off()
##PCA3
# png("pca3NAtl.png")
contour(nAtlLon,flip.matrix(nAtlLat), EOFann3E, lty='dotted', col='green', lwd=2, labcex=1, add=T)
world(add=T)
# dev.off()
x11()
#PCA1
# png("pca1NAtl.png")
contour(nAtlLon,flip.matrix(nAtlLat), EOFann1N, lwd=2, col='red', labcex=1, main="Wind stress EOFs N")
# world(add=T)
# dev.off()
##PCA2
# png("pca2NAtl.png")
contour(nAtlLon,flip.matrix(nAtlLat), EOFann2N, lty='dotdash', col='blue', lwd=2, labcex=1, add=T)
# world(add=T)
# dev.off()
##PCA3
# png("pca3NAtl.png")
contour(nAtlLon,flip.matrix(nAtlLat), EOFann3N, lty='dotted', col='green', lwd=2, labcex=1, add=T)
world(add=T)
# dev.off()

## Time series pattern
#quartz()
x11()
plot(ws.yrs,PCann1E, type='l', col='red', main="Eigenvector method time series pattern, E")
lines(ws.yrs,PCann2E, col='blue')
lines(ws.yrs,PCann3E, col='magenta')
x11()
plot(ws.yrs,PCann1N, type='l', col='red', main="Eigenvector method time series pattern, N")
lines(ws.yrs,PCann2N, col='blue')
lines(ws.yrs,PCann3N, col='magenta')
