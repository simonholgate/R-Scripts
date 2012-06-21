##****************************************************************************************************##
## Calculate EOFs over the N Atlantic using ACRE SLP data to give some feel for the data
## Just focus on N Atlantic 100W-15E and 5S-80N
##****************************************************************************************************##

##****************************************************************************************************##
########################################################################################################
## Functions for use below                                                                            ##
########################################################################################################
##****************************************************************************************************##

monthly.slp.ncdf <- function(path, nlon, nlat, nyr){
  ## Read the sea level data into an array
  nmon <- 12

  nc <- open.ncdf(path)
  slp.id <- nc$var[[2]]
  slpMonthlyArray <- get.var.ncdf(nc, slp.id)
  close.ncdf(nc)

  dim(slpMonthlyArray) <- c(nlon,nlat,nmon,nyr)
  slpMonthlyArray
}

##*******************************************************************************************************

natl.slp <- function(slpArray, lon, lat){
  ## Extract the N Atl region from the SLP array
  
  nlon <- dim(slpArray)[1]
  nlat <- dim(slpArray)[2]
  nyr <- dim(slpArray)[4]
  nmon <- 12
  

  nNAtlLon <- length(lon)
  nNAtlLat <- length(lat)
  
  dim(slpArray)<-c(nlon,nlat,nmon*nyr)
  
  slpNAtlantic <- slpArray[lon,lat,]
  dim(slpNAtlantic) <- c(nNAtlLon,nNAtlLat,nmon,nyr)
  
  slpNAtlantic
}

##*******************************************************************************************************

annual.slp <- function(slpArray){
  ## Convert a 4D array of monthly data into a 3D annual array
  
  nlon <- dim(slpArray)[1]
  nlat <- dim(slpArray)[2]
  nyr <- dim(slpArray)[4]
  
  slpAnnualArray <- array(NA, dim=c(nlon,nlat,nyr))
  for(i in 1:nlon){
    for(j in 1:nlat){
      slpAnnualArray[i,j,] <- colMeans(slpArray[i,j,,])
    }
  }
  
  slpAnnualArray
}

##*******************************************************************************************************

natl.lon.lat <- function(xlon, ylat){

  dlon <- xlon[2]-xlon[1]
  dlat <- ylat[2]-ylat[1]

  # Atlantic covers 100W to 15E and 5S to 80N
  # We need to see whether we are working in 0 to 360 or -180 to 180 co-ordinate system
  if(min(xlon)<0){
    nAtlLon <- match(seq(from=-100,to=15,by=dlon), xlon)
  } else {
    nAtlLon <- c(match(seq(from=260,to=359,by=dlon), xlon), match(seq(from=0,to=15,by=dlon), xlon))
  }
  # Do we go from N to S or S to N?
  if(dlat>0){ # S to N
    nAtlLat <- match(seq(from=-5,to=80, by=dlat), ylat)
  } else { # N to S
    nAtlLat <- match(seq(from=80,to=-5, by=dlat), ylat)
  }
  
  list(lon=nAtlLon,lat=nAtlLat)
}

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

nlon<-180
nlat<-91
nyr<-138

lon<-seq(from=0,by=2,length=nlon)
lon2 <- c(seq(from=0,by=2,to=180),seq(from=-178,by=2,to=-2))
lat<-seq(from=90,by=-2,length=nlat)
slp.yrs <- c(1871:2008)

##*******************************************************************************************************
## N Atlantic pressure

slpMonthlyArray <- monthly.slp.ncdf("~/data/ACRE/prmsl.mon.mean.nc", nlon, nlat, nyr)
# Convert Pa to Mb
slpMonthlyArray <- slpMonthlyArray/100

nAtlLonLatIndex <- natl.lon.lat(lon,lat)
nAtlLon <- lon2[nAtlLonLatIndex$lon]
nAtlLat <- lat[nAtlLonLatIndex$lat]

lnAtlLon <- length(nAtlLon)
lnAtlLat <- length(nAtlLat)
nstns <- lnAtlLon*lnAtlLat

slpNAtlArray <- natl.slp(slpMonthlyArray, nAtlLonLatIndex$lon, nAtlLonLatIndex$lat)
slpNAtlAnnArray <- annual.slp(slpNAtlArray)

##*******************************************************************************************************

## Prepare annual data for EOF caculations
dmSlpAnnual <- slpNAtlAnnArray

## Reshape to 2D array
dim(dmSlpAnnual) <- c(nstns,nyr)

## Rotate so that we have nstns columns
dmSlpAnnual <- t(dmSlpAnnual)

## Remove column means
cmeans <- colMeans(dmSlpAnnual, na.rm=T)
for(i in 1:nstns){
  dmSlpAnnual[,i] <- dmSlpAnnual[,i]-cmeans[i]
}


#######################################################################################################
## SVD method                                                                                        ##
#######################################################################################################

## Spatial loading pattern
svdDmSlpAnnual <- svd(dmSlpAnnual)

## Amount of variance explained is in d^2
svdVarExp <- svdDmSlpAnnual$d^2/sum(svdDmSlpAnnual$d^2)*100

## PCA1
PCAann1 <- svdDmSlpAnnual$v[,1]
dim(PCAann1) <- c(lnAtlLon, lnAtlLat)

## PCA2
PCAann2 <- svdDmSlpAnnual$v[,2]
dim(PCAann2) <- c(lnAtlLon, lnAtlLat)

## PCA3
PCAann3 <- svdDmSlpAnnual$v[,3]
dim(PCAann3) <- c(lnAtlLon, lnAtlLat)

## Find expansion coefficients
ECann1 <- dmSlpAnnual %*% svdDmSlpAnnual$v[,1]
ECann2 <- dmSlpAnnual %*% svdDmSlpAnnual$v[,2]
ECann3 <- dmSlpAnnual %*% svdDmSlpAnnual$v[,3]


#######################################################################################################
## Eigenvector method                                                                                ##
#######################################################################################################

## Form the co-variance matrix
covDmSlpAnnual <- t(dmSlpAnnual) %*% dmSlpAnnual

## Calculate eigenvalues and eigenvectors
eigCovDmSlpAnnual <- eigen(covDmSlpAnnual)

## % variance explained
eigVarExp <- eigCovDmSlpAnnual$values/sum(eigCovDmSlpAnnual$values)*100

## EOF1
EOFann1 <- eigCovDmSlpAnnual$vectors[,1]
dim(EOFann1) <- c(lnAtlLon, lnAtlLat)

## EOF2
EOFann2 <- eigCovDmSlpAnnual$vectors[,2]
dim(EOFann2) <- c(lnAtlLon, lnAtlLat)

## EOF3
EOFann3 <- eigCovDmSlpAnnual$vectors[,3]
dim(EOFann3) <- c(lnAtlLon, lnAtlLat)

## Find expansion coefficients
Cann <- diag(x=eigCovDmSlpAnnual$values,nrow=nstns, ncol=nstns)

PCann1 <- dmSlpAnnual %*% Cann[,1]
PCann2 <- dmSlpAnnual %*% Cann[,2]
PCann3 <- dmSlpAnnual %*% Cann[,3]

## PCann1 looks very odd - see plots below. Percentage variance explained and spatial patterns are
## identical. However, prefer SVD method on basis o time series reconstruction.

##***************************************************************************************************##
##***************************************************************************************************##
## Plots                                                                                             ##
##***************************************************************************************************##
##***************************************************************************************************##

#######################################################################################################
## Basic pressure map                                                                                ##
#######################################################################################################
x11()
contour(nAtlLon,flip.matrix(nAtlLat),mirror.matrix(slpNAtlAnnArray[,,1]))
world(add=T)

#######################################################################################################
## SVD method                                                                                        ##
#######################################################################################################

## Map
#quartz()
x11()
#map("worldHires", xlim=c(-100,15), ylim=c(-5,80), interior=F, fill=F, col="grey50", resolution=0)
#map.axes()
## PCA1
contour(nAtlLon,flip.matrix(nAtlLat), mirror.matrix(PCAann1), lwd=2, col='red')
## PCA2
contour(nAtlLon,flip.matrix(nAtlLat), mirror.matrix(PCAann2), lty='dotdash', add=T, col='blue')
## PCA3
contour(nAtlLon,flip.matrix(nAtlLat), mirror.matrix(PCAann3), lty='dotted', add=T, col='green')

world(add=T)

## Time series pattern
#quartz()
x11()
plot(slp.yrs,ECann1, type='l', col='red')
lines(slp.yrs,ECann2, col='blue')
lines(slp.yrs,ECann3, col='magenta')


#######################################################################################################
## Eigenvector method                                                                                ##
#######################################################################################################

## Map 
#quartz()
x11()
#map("worldHires", xlim=c(-100,15), ylim=c(-5,80), interior=F, fill=F, col="grey50", resolution=0)
#map.axes()
#PCA1
# png("pca1NAtl.png")
contour(nAtlLon,flip.matrix(nAtlLat), mirror.matrix(EOFann1), lwd=2, col='red', labcex=1)
world(add=T)
# dev.off()
##PCA2
# png("pca2NAtl.png")
contour(nAtlLon,flip.matrix(nAtlLat), mirror.matrix(EOFann2), lty='dotdash', col='blue', lwd=2, labcex=1)
world(add=T)
# dev.off()
##PCA3
# png("pca3NAtl.png")
contour(nAtlLon,flip.matrix(nAtlLat), mirror.matrix(EOFann3), lty='dotted', col='green', lwd=2, labcex=1)
world(add=T)
# dev.off()

## Time series pattern
#quartz()
x11()
plot(slp.yrs,PCann1, type='l', col='red')
lines(slp.yrs,PCann2, col='blue')
lines(slp.yrs,PCann3, col='magenta')
