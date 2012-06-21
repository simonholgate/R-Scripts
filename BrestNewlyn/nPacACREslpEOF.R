##****************************************************************************************************##
## Calculate EOFs over the N Pacific using ACRE SLP data to give some feel for the data
## Just focus on N Atlantic 90E-250E and 5S-80N
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

region.slp <- function(slpArray, lon, lat){
  ## Extract the N Atl region from the SLP array
  
  nlon <- dim(slpArray)[1]
  nlat <- dim(slpArray)[2]
  nyr <- dim(slpArray)[4]
  nmon <- 12
  

  nRegionLon <- length(lon)
  nRegionLat <- length(lat)
  
  dim(slpArray)<-c(nlon,nlat,nmon*nyr)
  
  slpRegion <- slpArray[lon,lat,]
  dim(slpRegion) <- c(nRegionLon,nRegionLat,nmon,nyr)
  
  slpRegion
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

npac.lon.lat <- function(xlon, ylat){
#  xlon<-seq(from=-180,by=5,length=nlon)
#  ylat<-seq(from=90,by=-5,length=nlat)
  dlon <- xlon[2]-xlon[1]
  dlat <- ylat[2]-ylat[1]

  # N Pacific covers 90E to 260E and 5S to 80N
  # We need to see whether we are working in 0 to 360 or -180 to 180 co-ordinate system
  if(min(xlon)<0){
    nPacLon <- c(match(seq(from=90,to=180,by=dlon), xlon), match(seq(from=-180,to=-100,by=dlon), xlon))
  } else {
    nPacLon <- match(seq(from=90,to=260,by=dlon), xlon)
  }
  # Do we go from N to S or S to N?
  if(dlat>0){ # S to N
    nPacLat <- match(seq(from=-5,to=80, by=dlat), ylat)
  } else { # N to S
    nPacLat <- match(seq(from=80,to=-5, by=dlat), ylat)
  }
  
  list(lon=nPacLon,lat=nPacLat)
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
## N Pacific pressure

slpMonthlyArray <- monthly.slp.ncdf("~/data/ACRE/prmsl.mon.mean.nc", nlon, nlat, nyr)
# Convert Pa to Mb
slpMonthlyArray <- slpMonthlyArray/100

nPacLonLatIndex <- npac.lon.lat(lon,lat)
nPacLon <- lon2[nPacLonLatIndex$lon]
nPacLat <- lat[nPacLonLatIndex$lat]

lnPacLon <- length(nPacLon)
lnPacLat <- length(nPacLat)
nstns <- lnPacLon*lnPacLat

slpNPacArray <- region.slp(slpMonthlyArray, nPacLonLatIndex$lon, nPacLonLatIndex$lat)
slpNPacAnnArray <- annual.slp(slpNPacArray)

##*******************************************************************************************************

## Prepare annual data for EOF caculations
dmSlpAnnual <- slpNPacAnnArray

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
dim(PCAann1) <- c(lnPacLon, lnPacLat)

## PCA2
PCAann2 <- svdDmSlpAnnual$v[,2]
dim(PCAann2) <- c(lnPacLon, lnPacLat)

## PCA3
PCAann3 <- svdDmSlpAnnual$v[,3]
dim(PCAann3) <- c(lnPacLon, lnPacLat)

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
dim(EOFann1) <- c(lnPacLon, lnPacLat)

## EOF2
EOFann2 <- eigCovDmSlpAnnual$vectors[,2]
dim(EOFann2) <- c(lnPacLon, lnPacLat)

## EOF3
EOFann3 <- eigCovDmSlpAnnual$vectors[,3]
dim(EOFann3) <- c(lnPacLon, lnPacLat)

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
lon2Pac <- nPacLon
lon2Pac[which(lon2Pac<0)] <- lon2Pac[which(lon2Pac<0)]+360

#######################################################################################################
## Basic pressure map                                                                                ##
#######################################################################################################
x11()
contour(lon2Pac,flip.matrix(nPacLat),mirror.matrix(slpNPacAnnArray[,,1]))
world(add=T, shift=T)

#######################################################################################################
## SVD method                                                                                        ##
#######################################################################################################

## Map
#quartz()
x11()
#map("worldHires", xlim=c(-100,15), ylim=c(-5,80), interior=F, fill=F, col="grey50", resolution=0)
#map.axes()
## PCA1
# png("pca1NPac.png")
contour(lon2Pac,flip.matrix(nPacLat), mirror.matrix(PCAann1), lwd=2, col='red', labcex=1)
world(add=T, shift=T)
# dev.off()
## PCA2
# png("pca2NPac.png")
contour(lon2Pac,flip.matrix(nPacLat), mirror.matrix(PCAann2), lty='dotdash',
        col='blue', lwd=2, labcex=1)
world(add=T, shift=T)
# dev.off()
## PCA3
# png("pca3NPac.png")
contour(lon2Pac,flip.matrix(nPacLat), mirror.matrix(PCAann3), lty='dotted',
        col='green', lwd=2, labcex=1)
world(add=T, shift=T)
# dev.off()

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
contour(lon2Pac,flip.matrix(nPacLat), mirror.matrix(EOFann1), lwd=2, col='red')
#PCA2
contour(lon2Pac,flip.matrix(nPacLat), mirror.matrix(EOFann2), lty='dotdash', add=T, col='blue')
#PCA3
contour(lon2Pac,flip.matrix(nPacLat), mirror.matrix(EOFann3), lty='dotted', add=T, col='green')

world(add=T, shift=T)

## Time series pattern
#quartz()
x11()
plot(slp.yrs,PCann1, type='l', col='red')
lines(slp.yrs,PCann2, col='blue')
lines(slp.yrs,PCann3, col='magenta')
