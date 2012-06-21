##****************************************************************************************************##
## Calculate EOFs over the N Atlantic using dynamic height calculated from Doug Smith's data.
## Just focus on N Atlantic 100W-15E and 5S-80N
##****************************************************************************************************##

##****************************************************************************************************##
########################################################################################################
## Functions for use below                                                                            ##
########################################################################################################
##****************************************************************************************************##

## monthly.dh.ncdf <- function(path, nlon, nlat, nyr){
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

## natl.dh <- function(slpArray, lon, lat){
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

natl.dh <- function(dataArray, lon, lat){
  ## Extract the N Atl region from the annual mean  wind stress array
  dataNAtlantic <- dataArray[lon,lat,]
  dataNAtlantic
}

##*****************************************************************************

natl.lon.lat.dh <- function(xlon, ylat){

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
##library(maps)
##library(mapdata)
##library(ncdf)
library(robust)
source("/home/simonh/workspace/RScripts/matrixMethods.R")

full.p.yrs <- c(1916:2008)

nlonds<-93
nlatds<-70
nyrds<-60

nlonen<-116
nlaten<-86
nyren<-45

nlonis <- 117
nlatis <- 87
nyris <- 66

dhds.yrs <- c(1950:2009)
dhen.yrs <- c(1966:2010)
dhis.yrs <- c(1945:2010)

## nlon<-192
## nlat<-94
## nyr<-138

## lon<-seq(from=0,by=2,length=nlon)
## lon2 <- c(seq(from=0,by=2,to=180),seq(from=-178,by=2,to=-2))
## lon<-seq(from=0,by=1.875,length=nlon)
## lon2 <- c(seq(from=0,by=1.875,to=180),seq(from=-178.125,by=1.875,to=-1.875))

## lat<-seq(from=90,by=-2,length=nlat)     
## dhlon <- seq(from=-180, to=178.125, by=1.875)
## lat <- c(88.542, 86.6531, 84.7532, 82.8508, 80.9473, 79.0435, 77.1394, 75.2351, 
##     73.3307, 71.4262, 69.5217, 67.6171, 65.7125, 63.8079, 61.9033, 59.9986, 
##     58.0939, 56.1893, 54.2846, 52.3799, 50.4752, 48.5705, 46.6658, 44.7611, 
##     42.8564, 40.9517, 39.047, 37.1422, 35.2375, 33.3328, 31.4281, 29.5234, 
##     27.6186, 25.7139, 23.8092, 21.9044, 19.9997, 18.095, 16.1902, 14.2855, 
##     12.3808, 10.47604, 8.57131, 6.66657, 4.76184, 2.8571, 0.952368, 
##     -0.952368, -2.8571, -4.76184, -6.66657, -8.57131, -10.47604, -12.3808, 
##     -14.2855, -16.1902, -18.095, -19.9997, -21.9044, -23.8092, -25.7139, 
##     -27.6186, -29.5234, -31.4281, -33.3328, -35.2375, -37.1422, -39.047, 
##     -40.9517, -42.8564, -44.7611, -46.6658, -48.5705, -50.4752, -52.3799, 
##     -54.2846, -56.1893, -58.0939, -59.9986, -61.9033, -63.8079, -65.7125, 
##     -67.6171, -69.5217, -71.4262, -73.3307, -75.2351, -77.1394, -79.0435, 
##     -80.9473, -82.8508, -84.7532, -86.6531, -88.542)

## dh.yrs <- c(1871:2008)

## newlyn.start.year.partial <- 1953
## newlyn.end.year.partial <- 2008
## newlyn.start.year.pred <- 1916
## newlyn.end.year.pred <- 1943
## newlyn.start.year.full <- 1916
## newlyn.end.year.full <- 2008

## Doug Smith
newlyn.start.year.partial.ds <- 1965
newlyn.end.year.partial.ds <- 2008
newlyn.start.year.pred.ds <- 1950
newlyn.end.year.pred.ds <- 1964
newlyn.start.year.full.ds <- 1950
newlyn.end.year.full.ds <- 2008
## Consistent dates across all sets - last 20 years
newlyn.start.year.cons.ds <- 1988
newlyn.end.year.cons.ds <- 2008

full.p.start.partial.ds <- which(full.p.yrs==newlyn.start.year.partial.ds)
full.p.end.partial.ds <- which(full.p.yrs==newlyn.end.year.partial.ds)
full.p.start.pred.ds <- which(full.p.yrs==newlyn.start.year.pred.ds)
full.p.end.pred.ds <- which(full.p.yrs==newlyn.end.year.pred.ds)
full.p.start.full.ds <- which(full.p.yrs==newlyn.start.year.full.ds)
full.p.end.full.ds <- which(full.p.yrs==newlyn.end.year.full.ds)
## Consistent dates across all sets - last 20 years
full.p.start.cons.ds <- which(full.p.yrs==newlyn.start.year.cons.ds)
full.p.end.cons.ds <- which(full.p.yrs==newlyn.end.year.cons.ds)

## EN3c
newlyn.start.year.partial.en <- 1980
newlyn.end.year.partial.en <- 2008
newlyn.start.year.pred.en <- 1966
newlyn.end.year.pred.en <- 1979
newlyn.start.year.full.en <- 1966
newlyn.end.year.full.en <- 2008
## Consistent dates across all sets - last 20 years
newlyn.start.year.cons.en <- 1988
newlyn.end.year.cons.en <- 2008

full.p.start.partial.en <- which(full.p.yrs==newlyn.start.year.partial.en)
full.p.end.partial.en <- which(full.p.yrs==newlyn.end.year.partial.en)
full.p.start.pred.en <- which(full.p.yrs==newlyn.start.year.pred.en)
full.p.end.pred.en <- which(full.p.yrs==newlyn.end.year.pred.en)
full.p.start.full.en <- which(full.p.yrs==newlyn.start.year.full.en)
full.p.end.full.en <- which(full.p.yrs==newlyn.end.year.full.en)
## Consistent dates across all sets - last 20 years
full.p.start.cons.en <- which(full.p.yrs==newlyn.start.year.cons.en)
full.p.end.cons.en <- which(full.p.yrs==newlyn.end.year.cons.en)

## Ishii
newlyn.start.year.partial.is <- 1965
newlyn.end.year.partial.is <- 2008
newlyn.start.year.pred.is <- 1945
newlyn.end.year.pred.is <- 1964
newlyn.start.year.full.is <- 1945
newlyn.end.year.full.is <- 2008
## Consistent dates across all sets - last 20 years
newlyn.start.year.cons.is <- 1988
newlyn.end.year.cons.is <- 2008

full.p.start.partial.is <- which(full.p.yrs==newlyn.start.year.partial.is)
full.p.end.partial.is <- which(full.p.yrs==newlyn.end.year.partial.is)
full.p.start.pred.is <- which(full.p.yrs==newlyn.start.year.pred.is)
full.p.end.pred.is <- which(full.p.yrs==newlyn.end.year.pred.is)
full.p.start.full.is <- which(full.p.yrs==newlyn.start.year.full.is)
full.p.end.full.is <- which(full.p.yrs==newlyn.end.year.full.is)
## Consistent dates across all sets - last 20 years
full.p.start.cons.is <- which(full.p.yrs==newlyn.start.year.cons.is)
full.p.end.cons.is <- which(full.p.yrs==newlyn.end.year.cons.is)

##*****************************************************************************
## N Atlantic steric height

## load("~/data/ACRE/dhACRE.RData")
load("~/Dropbox/DynamicHeight/dh_DS_EN3c_IS.RData")

## nAtlLonLatIndex <- natl.lon.lat.dh(dhlon,flip.matrix(lat))
## dhNAtlArrayE <- natl.dh(dhEAnnualMean, nAtlLonLatIndex$lon,
##                          nAtlLonLatIndex$lat)
## dhNAtlArrayN <- natl.dh(dhNAnnualMean, nAtlLonLatIndex$lon,
##                          nAtlLonLatIndex$lat)

## Mean steric height
## lenLon <- length(nAtlLonLatIndex$lon)
## lenLat <- length(nAtlLonLatIndex$lat)
## mdhNAtlArrayE <- array(NA, dim=c(lenLon, lenLat))
## mdhNAtlArrayN <- array(NA, dim=c(lenLon, lenLat))
## for(i in 1:lenLon){
##   for(j in 1:lenLat){
##     mdhNAtlArrayE[i,j] <- mean(dhNAtlArrayE[i,j,], na.rm=T)
##     mdhNAtlArrayN[i,j] <- mean(dhNAtlArrayN[i,j,], na.rm=T)
##   }
## }

## ##*******************************************************************************************************
## ## N Atlantic pressure

## slpMonthlyArray <- monthly.slp.ncdf("~/data/ACRE/prmsl.mon.mean.nc", nlon, nlat, nyr)
## # Convert Pa to Mb
## slpMonthlyArray <- slpMonthlyArray/100

##nAtlLonLatIndex <- natl.lon.lat.dh(lon,lat)
## nAtlLon <- dhlon[nAtlLonLatIndex$lon]
## nAtlLat <- lat[nAtlLonLatIndex$lat]

## lnAtlLon <- length(nAtlLon)
## lnAtlLat <- length(nAtlLat)

nstns.ds <- nlonds*nlatds
nstns.en <- nlonen*nlaten
nstns.is <- nlonis*nlatis

## EdhNAtlAnnArray <- natl.dh(dhNatlArray, nAtlLonLatIndex$lon, nAtlLonLatIndex$lat)
##  <- annual.(slpNAtlArray)

## ##*******************************************************************************************************

## Prepare annual data for EOF caculations
dmdhdsAnnual <- dhdsAnnualMean
dmdhenAnnual <- dhenAnnualMean
dmdhisAnnual <- dhisAnnualMean

## Reshape to 2D array
dim(dmdhdsAnnual) <- c(nstns.ds,nyrds)
dim(dmdhenAnnual) <- c(nstns.en,nyren)
dim(dmdhisAnnual) <- c(nstns.is,nyris)

## Rotate so that we have nstns columns
dmdhdsAnnual <- t(dmdhdsAnnual)
dmdhenAnnual <- t(dmdhenAnnual)
dmdhisAnnual <- t(dmdhisAnnual)

## Remove column means
cdsmeans <- colMeans(dmdhdsAnnual, na.rm=T)
cenmeans <- colMeans(dmdhenAnnual, na.rm=T)
cismeans <- colMeans(dmdhisAnnual, na.rm=T)

for(i in 1:nstns.ds){
  dmdhdsAnnual[,i] <- dmdhdsAnnual[,i]-cdsmeans[i]
}
for(i in 1:nstns.en){
  dmdhenAnnual[,i] <- dmdhenAnnual[,i]-cenmeans[i]
}
for(i in 1:nstns.is){
  dmdhisAnnual[,i] <- dmdhisAnnual[,i]-cismeans[i]
}

## Remove columns that are all NA
fincols.ds <- which(is.finite(dmdhdsAnnual[1,]))
dmdhdsAnnual <- dmdhdsAnnual[,fincols.ds]

fincols.en <- which(is.finite(dmdhenAnnual[1,]))
dmdhenAnnual <- dmdhenAnnual[,fincols.en]

fincols.is <- which(is.finite(dmdhisAnnual[1,]))
dmdhisAnnual <- dmdhisAnnual[,fincols.is]

#######################################################################################################
## SVD method                                                                                        ##
#######################################################################################################

## Spatial loading pattern
svdDmdhdsAnnual <- svd(dmdhdsAnnual)
svdDmdhenAnnual <- svd(dmdhenAnnual)
svdDmdhisAnnual <- svd(dmdhisAnnual)

## Amount of variance explained is in d^2
svdVarExpDS <- svdDmdhdsAnnual$d^2/sum(svdDmdhdsAnnual$d^2)*100
svdVarExpEN <- svdDmdhenAnnual$d^2/sum(svdDmdhenAnnual$d^2)*100
svdVarExpIS <- svdDmdhisAnnual$d^2/sum(svdDmdhisAnnual$d^2)*100

## PCA1
PCAann1DS <- vector(mode="numeric", length=nstns.ds)
PCAann1DS[fincols.ds] <- svdDmdhdsAnnual$v[,1]
dim(PCAann1DS) <- c(nlonds, nlatds)
PCAann1EN <- vector(mode="numeric", length=nstns.en)
PCAann1EN[fincols.en] <- svdDmdhenAnnual$v[,1]
dim(PCAann1EN) <- c(nlonen, nlaten)
PCAann1IS <- vector(mode="numeric", length=nstns.is)
PCAann1IS[fincols.is] <- svdDmdhisAnnual$v[,1]
dim(PCAann1IS) <- c(nlonis, nlatis)

## PCA2
PCAann2DS <- vector(mode="numeric", length=nstns.ds)
PCAann2DS[fincols.ds] <- svdDmdhdsAnnual$v[,2]
dim(PCAann2DS) <- c(nlonds, nlatds)
PCAann2EN <- vector(mode="numeric", length=nstns.en)
PCAann2EN[fincols.en] <- svdDmdhenAnnual$v[,2]
dim(PCAann2EN) <- c(nlonen, nlaten)
PCAann2IS <- vector(mode="numeric", length=nstns.is)
PCAann2IS[fincols.is] <- svdDmdhisAnnual$v[,2]
dim(PCAann2IS) <- c(nlonis, nlatis)

## PCA3
PCAann3DS <- vector(mode="numeric", length=nstns.ds)
PCAann3DS[fincols.ds] <- svdDmdhdsAnnual$v[,3]
dim(PCAann3DS) <- c(nlonds, nlatds)
PCAann3EN <- vector(mode="numeric", length=nstns.en)
PCAann3EN[fincols.en] <- svdDmdhenAnnual$v[,3]
dim(PCAann3EN) <- c(nlonen, nlaten)
PCAann3IS <- vector(mode="numeric", length=nstns.is)
PCAann3IS[fincols.is] <- svdDmdhisAnnual$v[,3]
dim(PCAann3IS) <- c(nlonis, nlatis)

## Find expansion coefficients
ECann1DS <- dmdhdsAnnual %*% svdDmdhdsAnnual$v[,1]
ECann2DS <- dmdhdsAnnual %*% svdDmdhdsAnnual$v[,2]
ECann3DS <- dmdhdsAnnual %*% svdDmdhdsAnnual$v[,3]

ECann1EN <- dmdhenAnnual %*% svdDmdhenAnnual$v[,1]
ECann2EN <- dmdhenAnnual %*% svdDmdhenAnnual$v[,2]
ECann3EN <- dmdhenAnnual %*% svdDmdhenAnnual$v[,3]

ECann1IS <- dmdhisAnnual %*% svdDmdhisAnnual$v[,1]
ECann2IS <- dmdhisAnnual %*% svdDmdhisAnnual$v[,2]
ECann3IS <- dmdhisAnnual %*% svdDmdhisAnnual$v[,3]


## #######################################################################################################
## ## Eigenvector method                                                                                ##
## #######################################################################################################

## ## Form the co-variance matrix
## covDmEdhAnnual <- t(dmEdhAnnual) %*% dmEdhAnnual
## covDmNdhAnnual <- t(dmNdhAnnual) %*% dmNdhAnnual

## ## Calculate eigenvalues and eigenvectors
## eigCovDmEdhAnnual <- eigen(covDmEdhAnnual)
## eigCovDmNdhAnnual <- eigen(covDmNdhAnnual)

## ## % variance explained
## eigVarExpE <- eigCovDmEdhAnnual$values/sum(eigCovDmEdhAnnual$values)*100
## eigVarExpN <- eigCovDmNdhAnnual$values/sum(eigCovDmNdhAnnual$values)*100

## ## EOF1
## EOFann1E <- eigCovDmEdhAnnual$vectors[,1]
## dim(EOFann1E) <- c(lenLon, lenLat)
## EOFann1N <- eigCovDmNdhAnnual$vectors[,1]
## dim(EOFann1N) <- c(lenLon, lenLat)

## ## EOF2
## EOFann2E <- eigCovDmEdhAnnual$vectors[,2]
## dim(EOFann2E) <- c(lenLon, lenLat)
## EOFann2N <- eigCovDmNdhAnnual$vectors[,2]
## dim(EOFann2N) <- c(lenLon, lenLat)

## ## EOF3
## EOFann3E <- eigCovDmEdhAnnual$vectors[,3]
## dim(EOFann3E) <- c(lenLon, lenLat)
## EOFann3N <- eigCovDmNdhAnnual$vectors[,3]
## dim(EOFann3N) <- c(lenLon, lenLat)

## ## Find expansion coefficients
## CannE <- diag(x=eigCovDmEdhAnnual$values,nrow=nstns, ncol=nstns)
## CannN <- diag(x=eigCovDmNdhAnnual$values,nrow=nstns, ncol=nstns)

## PCann1E <- dmEdhAnnual %*% CannE[,1]
## PCann2E <- dmEdhAnnual %*% CannE[,2]
## PCann3E <- dmEdhAnnual %*% CannE[,3]
## PCann1N <- dmNdhAnnual %*% CannN[,1]
## PCann2N <- dmNdhAnnual %*% CannN[,2]
## PCann3N <- dmNdhAnnual %*% CannN[,3]

## PCann1 looks very odd - see plots below. Percentage variance explained and spatial patterns are
## identical. However, prefer SVD method on basis o time series reconstruction.

##***************************************************************************************************##
##***************************************************************************************************##
## Plots                                                                                             ##
##***************************************************************************************************##
##***************************************************************************************************##

#######################################################################################################
## Basic height map                                                                                  ##
#######################################################################################################
x11()
filled.contour(londs,latds,dhdsAnnualMean[,,1], main="Steric Height DS, 1871",
               color=tim.colors, nlevels=20, plot.axes = { world(add=T) }, zlim=c(-16, 29))

x11()
filled.contour(lonen,laten,dhenAnnualMean[,,1], main="Steric Height EN, 1871",
               color=tim.colors, nlevels=20, plot.axes = { world(add=T) }, zlim=c(-16, 29))

x11()
filled.contour(lonis,latis,dhisAnnualMean[,,1], main="Steric Height IS, 1871",
               color=tim.colors, nlevels=20, plot.axes = { world(add=T) }, zlim=c(-16, 29))

#######################################################################################################
## SVD method                                                                                        ##
#######################################################################################################

## Maps
x11()
## PCA1
contour(londs,latds, PCAann1DS, lwd=2, col='red', main="Steric Height PCAs DS")
## PCA2
contour(londs,latds, PCAann2DS, lty='dotdash', add=T, col='blue')
## PCA3
contour(londs,latds, PCAann3DS, lty='dotted', add=T, col='green')
world(add=T)

x11()
## PCA1
contour(lonen,laten, PCAann1EN, lwd=2, col='red', main="Steric Height PCAs EN")
## PCA2
contour(lonen,laten, PCAann2EN, lty='dotdash', add=T, col='blue')
## PCA3
contour(lonen,laten, PCAann3EN, lty='dotted', add=T, col='green')
world(add=T)

x11()
## PCA1
contour(lonis,latis, PCAann1IS, lwd=2, col='red', main="Steric Height PCAs IS")
## PCA2
contour(lonis,latis, PCAann2IS, lty='dotdash', add=T, col='blue')
## PCA3
contour(lonis,latis, PCAann3IS, lty='dotted', add=T, col='green')
world(add=T)


## Time series pattern
x11()
plot(dhds.yrs,ECann1DS, type='l', col='red', main="SVD method time series pattern, DS")
lines(dhds.yrs,ECann2DS, col='blue')
lines(dhds.yrs,ECann3DS, col='magenta')

x11()
plot(dhen.yrs,ECann1EN, type='l', col='red', main="SVD method time series pattern, EN")
lines(dhen.yrs,ECann2EN, col='blue')
lines(dhen.yrs,ECann3EN, col='magenta')

x11()
plot(dhis.yrs,ECann1IS, type='l', col='red', main="SVD method time series pattern, IS")
lines(dhis.yrs,ECann2IS, col='blue')
lines(dhis.yrs,ECann3IS, col='magenta')


## #######################################################################################################
## ## Eigenvector method                                                                                ##
## #######################################################################################################

## ## Map 
## #quartz()
## x11()
## #PCA1
## # png("pca1NAtl.png")
## contour(nAtlLon,flip.matrix(nAtlLat), EOFann1E, lwd=2, col='red', labcex=1, main="Wind stress EOFs E")
## # world(add=T)
## # dev.off()
## ##PCA2
## # png("pca2NAtl.png")
## contour(nAtlLon,flip.matrix(nAtlLat), EOFann2E, lty='dotdash', col='blue', lwd=2, labcex=1, add=T)
## # world(add=T)
## # dev.off()
## ##PCA3
## # png("pca3NAtl.png")
## contour(nAtlLon,flip.matrix(nAtlLat), EOFann3E, lty='dotted', col='green', lwd=2, labcex=1, add=T)
## world(add=T)
## # dev.off()
## x11()
## #PCA1
## # png("pca1NAtl.png")
## contour(nAtlLon,flip.matrix(nAtlLat), EOFann1N, lwd=2, col='red', labcex=1, main="Wind stress EOFs N")
## # world(add=T)
## # dev.off()
## ##PCA2
## # png("pca2NAtl.png")
## contour(nAtlLon,flip.matrix(nAtlLat), EOFann2N, lty='dotdash', col='blue', lwd=2, labcex=1, add=T)
## # world(add=T)
## # dev.off()
## ##PCA3
## # png("pca3NAtl.png")
## contour(nAtlLon,flip.matrix(nAtlLat), EOFann3N, lty='dotted', col='green', lwd=2, labcex=1, add=T)
## world(add=T)
## # dev.off()

## ## Time series pattern
## #quartz()
## x11()
## plot(dh.yrs,PCann1E, type='l', col='red', main="Eigenvector method time series pattern, E")
## lines(dh.yrs,PCann2E, col='blue')
## lines(dh.yrs,PCann3E, col='magenta')
## x11()
## plot(dh.yrs,PCann1N, type='l', col='red', main="Eigenvector method time series pattern, N")
## lines(dh.yrs,PCann2N, col='blue')
## lines(dh.yrs,PCann3N, col='magenta')
