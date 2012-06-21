##****************************************************************************************************##
## Calculate EOFs over the Brest-Newlyn area using ACRE SLP data
## Just focus on 20W-10E and 40N-60N
##****************************************************************************************************##

##****************************************************************************************************##
########################################################################################################
## Functions for use below                                                                            ##
########################################################################################################
##****************************************************************************************************##

annual.2d.slp <- function(slpArray){
  ## Convert a 2D array of monthly data into a 2D annual array.
  ## There should be nstns columns and nmon rows
  
  nstns <- dim(slpArray)[2]
  nmon <- dim(slpArray)[1]
  nyr <- nmon/12


  slpAnnualArray <- array(NA, dim=c(nyr,nstns))

  for(i in 1:nstns){
      ## Make a temporray vector of each station
      temp <- slpArray[,i]
      ## Reshape to 3d array with nstns*nyr*12
      dim(temp) <- c(12, nyr)
      ## Place the nyr column means from the temporary array into the column for stn i
      slpAnnualArray[,i] <- colMeans(temp)
  }
  
  slpAnnualArray
}

##*******************************************************************************************************

#########################
## Non-functional part ##
#########################

library(fields)
#library(robust)
#source("/home/simonh/workspace/RScripts/matrixMethods.R")

nyr<-138

#lat <- c(48.38, 50.1,
#                  40, 40, 40, 40, 40, 40, 40,
#                  45, 45, 45, 45, 45, 45, 45,
#                  50, 50, 50, 50, 50, 50, 50, 
#                  55, 55, 55, 55, 55, 55, 55, 
#                  60, 60, 60, 60, 60, 60, 60)
#lon <- c(-4.5, -5.55,
#                  -20, -15, -10, -5, 0, 5, 10,
#                  -20, -15, -10, -5, 0, 5, 10,
#                  -20, -15, -10, -5, 0, 5, 10,
#                  -20, -15, -10, -5, 0, 5, 10,
#                  -20, -15, -10, -5, 0, 5, 10)

lat <- seq(from=40, to=60, by = 5)
lon <- seq(from=-20, to=10, by=5)

slp.yrs <- c(1871:2008)

##*******************************************************************************************************
## Brest-Newlyn area pressure

## Monthly data 37 stations. First two are Brest and Newlyn. Variable name is slpNewlynStns.
load("/home/simonh/Dropbox/ACRE/brestNewlynACREslp.RData")

# Convert Pa to Mb
slpNewlynStns <- slpNewlynStns/100

lLon <- length(lon)
lLat <- length(lat)
nstns <- 37

#slpNAtlArray <- natl.slp(slpMonthlyArray, nAtlLonLatIndex$lon, nAtlLonLatIndex$lat)
slpAnnArray <- annual.2d.slp(slpNewlynStns)


##*******************************************************************************************************

## Prepare annual data for EOF calculations
## Drop Brest and Newlyn from the array
dmSlpAnnual <- slpAnnArray[,3:nstns]

## Remove column means
cmeans <- colMeans(dmSlpAnnual, na.rm=T)
for(i in 1:(nstns-2)){
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
dim(PCAann1) <- c(lLon, lLat)

## PCA2
PCAann2 <- svdDmSlpAnnual$v[,2]
dim(PCAann2) <- c(lLon, lLat)

## PCA3
PCAann3 <- svdDmSlpAnnual$v[,3]
dim(PCAann3) <- c(lLon, lLat)

## Find expansion coefficients
ECann1 <- dmSlpAnnual %*% svdDmSlpAnnual$v[,1]
ECann2 <- dmSlpAnnual %*% svdDmSlpAnnual$v[,2]
ECann3 <- dmSlpAnnual %*% svdDmSlpAnnual$v[,3]

## Save the data
save(ECann1, ECann2, ECann3, PCAann1, PCAann2, PCAann3, dmSlpAnnual, svdDmSlpAnnual, file="brestNewlynACREslpEOF.RData")
save(slpAnnArray, file="brestNewlynACREannualSLP.RData")
     
##***************************************************************************************************##
##***************************************************************************************************##
## Plots                                                                                             ##
##***************************************************************************************************##
##***************************************************************************************************##

#######################################################################################################
## Basic pressure map                                                                                ##
#######################################################################################################
x11()
slpAnnArray1 <- slpAnnArray[1,3:nstns]
dim(slpAnnArray1) <- c(lLon,lLat)
contour(lon,lat,slpAnnArray1)
world(add=T)

#######################################################################################################
## SVD method                                                                                        ##
#######################################################################################################

## Map
#quartz()
x11()
## PCA1
contour(lon,lat, PCAann1, lwd=2, col='red')
## PCA2
contour(lon,lat, PCAann2, lty='dotdash', add=T, col='blue')
## PCA3
contour(lon,lat, PCAann3, lty='dotted', add=T, col='green')

world(add=T)

## Time series pattern
#quartz()
x11()
plot(slp.yrs,ECann1, type='l', col='red')
lines(slp.yrs,ECann2, col='blue')
lines(slp.yrs,ECann3, col='magenta')
