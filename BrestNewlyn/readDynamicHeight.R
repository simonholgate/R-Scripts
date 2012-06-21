###############################################################################
## readDynamicHeight.R
## 
## Reads the fields of annual mean dynamic height (calculated from the monthly 
## temperature and salinity fields) for the North Atlantic. 
## Dynamic heights are calculated from three datasets: Ishii, EN3c and Doug Smith's
## The data here is stored as float64 (real*8) written out from Matlab.
## Note: The values here are in m2/s2.
##
## Author: simonh
###############################################################################

# Newlyn 50 06 N  05 33 W => 50.1 -5.55
# Brest 48 23 N  04 30 W

library(fields)
source("~/Dropbox/BrestNewlyn/matrixMethods.R")

## Doug Smith data. Reference depth is level 15 at 2116.15m. Runs 1950-2009 (60 years)
lds <- 93
mds <- 70
lmds <- lds*mds
numYearsDS <- 60

latds <- seq(from=-5.625, to=80.625, by=1.25)
londs <- seq(from=-100, to=15, by=1.25)

dhds <- file("dhDS_NAtlantic.bin", "rb")
dhdsAnnualMean <- array(NA, dim=c(lds, mds, numYearsDS))

## EN3c data. Reference depth is level L31 at 2133.5m. Runs 1966-2010 (45 years)
len <- 116
men <- 86
lmen <- len*men
numYearsEN <- 45

laten <- seq(from=-5, length=men, by=1)
lonen <- seq(from=-100, length=len, by=1)

dhen <- file("dhEN3c_NAtlantic.bin", "rb")
dhenAnnualMean <- array(NA, dim=c(len, men, numYearsEN))

## Ishii data. Reference depth is level L24 at 1500m. Runs 1945-2010 (66 years)
lis <- 117
mis <- 87
lmis <- lis*mis
numYearsIS <- 66

latis <- seq(from=-5.5, length=mis, by=1)
lonis <- seq(from=-100.5, length=lis, by=1)

dhis <- file("dhIshii_NAtlantic.bin", "rb")
dhisAnnualMean <- array(NA, dim=c(lis, mis, numYearsIS))

## Read DS data
for (i in 1:numYearsDS){
  junk <- readBin(dhds, "numeric", size=8, n = lmds)
  dim(junk) <- c(lds,mds)
  dhdsAnnualMean[,,i] <- junk
}
close(dhds)

## Plot DS data
x11()
filled.contour(londs,latds,dhdsAnnualMean[,,1], color=tim.colors, nlevels=20,
 plot.axes = { world(shift=F, add=T) }, zlim=c(-17, 25))


## Read EN3c data
for (i in 1:numYearsEN){
  junk <- readBin(dhen, "numeric", size=8, n = lmen)
  dim(junk) <- c(len,men)
  dhenAnnualMean[,,i] <- junk
}
close(dhen)

## Plot EN3c data
x11()
filled.contour(lonen,laten,dhenAnnualMean[,,1], color=tim.colors, nlevels=20,
 plot.axes = { world(shift=F, add=T) }, zlim=c(-17, 25))

## Read Ishii data
for (i in 1:numYearsIS){
  junk <- readBin(dhis, "numeric", size=8, n = lmis)
  dim(junk) <- c(lis,mis)
  dhisAnnualMean[,,i] <- junk
}
close(dhis)

## Replace zeros with NAs
nulls <- which(dhisAnnualMean==0)
dhisAnnualMean[nulls] <- NA

## Plot Ishii data
x11()
filled.contour(lonis,latis,dhisAnnualMean[,,1], color=tim.colors, nlevels=20,
 plot.axes = { world(shift=F, add=T) }, zlim=c(-17, 25))

# Save the data
save(file="~/Dropbox/DynamicHeight/dh_DS_EN3c_IS.RData", lonen,laten,dhenAnnualMean,
     londs,latds,dhdsAnnualMean, lonis,latis,dhisAnnualMean)
