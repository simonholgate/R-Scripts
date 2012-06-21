## Function to correlate sea level corrected for pressure at Key West (total pressure) with sea level
## pressure everywhere

source("~/Dropbox/BrestNewlyn/correlateFunctionsACRE.R")

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
#lat2<-seq(from=-90,by=2,length=nlat)
slp.yrs <- c(1871:2009)

keywest.start.year.partial <- 1953
keywest.end.year.partial <- 2008
keywest.start.year.pred <- 1916
keywest.end.year.pred <- 1952
keywest.start.year.full <- 1916
keywest.end.year.full <- 2008
keywest.lon <- -81.8
keywest.lat <- 24.55

##*******************************************************************************************************
## N Atlantic pressure

slpMonthlyArray <- monthly.slp.ncdf("~/data/ACRE/prmsl.mon.mean.nc", nlon, nlat, nyr)
# Convert Pa to Mb
slpMonthlyArray <- slpMonthlyArray/100

nAtlLonLatIndex <- natl.lon.lat(lon,lat)
slpNAtlArray <- region.slp(slpMonthlyArray, nAtlLonLatIndex$lon, nAtlLonLatIndex$lat)
slpNAtlAnnArray <- annual.slp(slpNAtlArray)

##*******************************************************************************************************
## N Pacific pressure

nPacLonLatIndex <- npac.lon.lat(lon,lat)
slpNPacArray <- region.slp(slpMonthlyArray, nPacLonLatIndex$lon, nPacLonLatIndex$lat)
slpNPacAnnArray <- annual.slp(slpNPacArray)

##*******************************************************************************************************
## Key West

## Partial
keywest.sl.partial <- load.tg.data(keywest.start.year.partial,keywest.end.year.partial,"keywest")
keywest.pa.partial <-
  interp.local.pressure(keywest.lon,keywest.lat,slpNAtlAnnArray,
                        range(lon2[nAtlLonLatIndex$lon]),range(lat[nAtlLonLatIndex$lat]),
                         keywest.start.year.partial, keywest.end.year.partial, slp.yrs)

keywest.tot.p.partial <- total.pressure(keywest.sl.partial$Height,keywest.pa.partial)

keywest.corr.partial <- corr.pressure(keywest.tot.p.partial, slpNAtlAnnArray,
                             keywest.start.year.partial, keywest.end.year.partial, slp.yrs)

keywest.corr.dt.partial <- corr.pressure.dt(keywest.tot.p.partial, slpNAtlAnnArray,
                             keywest.start.year.partial, keywest.end.year.partial, slp.yrs)

## Prediction
keywest.sl.pred <- load.tg.data(keywest.start.year.pred,keywest.end.year.pred,"keywest")
keywest.pa.pred <-
  interp.local.pressure(keywest.lon,keywest.lat,slpNAtlAnnArray,range(lon2[nAtlLonLatIndex$lon]),
                        range(lat[nAtlLonLatIndex$lat]), keywest.start.year.pred,
                        keywest.end.year.pred, slp.yrs)

keywest.tot.p.pred <- total.pressure(keywest.sl.pred$Height,keywest.pa.pred)

keywest.corr.pred <- corr.pressure(keywest.tot.p.pred, slpNAtlAnnArray,
                             keywest.start.year.pred, keywest.end.year.pred, slp.yrs)

## Full
keywest.sl.full <- load.tg.data(keywest.start.year.full,keywest.end.year.full,"keywest")
keywest.pa.full <-
  interp.local.pressure(keywest.lon,keywest.lat,slpNAtlAnnArray,
                        range(lon2[nAtlLonLatIndex$lon]),range(lat[nAtlLonLatIndex$lat]),
                         keywest.start.year.full, keywest.end.year.full, slp.yrs)

keywest.tot.p.full <- total.pressure(keywest.sl.full$Height,keywest.pa.full)

keywest.corr.full <- corr.pressure(keywest.tot.p.full, slpNAtlAnnArray,
                             keywest.start.year.full, keywest.end.year.full, slp.yrs)

keywest.corr.dt.full <- corr.pressure.dt(keywest.tot.p.full, slpNAtlAnnArray,
                             keywest.start.year.full, keywest.end.year.full, slp.yrs)

##*******************************************************************************************************
## Sort correlated data

lmRobControl <- lmRob.control(mxr=100,mxf=100,trace=F)

## Partial
keywest.sorted.pa.partial <- sort.corr.pa(keywest.corr.partial[,,1], slpNAtlAnnArray,
                                         keywest.start.year.partial, keywest.end.year.partial, slp.yrs)
data.keywest.partial <- data.frame(msl=keywest.tot.p.partial,
                                  t=seq(from=keywest.start.year.partial,to=keywest.end.year.partial),
                                  keywest.sorted.pa.partial)
tg.lmRob.keywest.partial <- lmRob(msl ~ ., x=T, y=T, data=data.keywest.partial, control=lmRobControl,
                                 na.action=na.exclude)
## Variance reduction
tg.lmRob.keywest.partial$r.sq
#[1] 0.7279121

## Prediction
keywest.sorted.pa.pred <- sort.corr.pa(keywest.corr.pred[,,1],
                                      slpNAtlAnnArray,
                                      keywest.start.year.pred,
                                      keywest.end.year.pred, slp.yrs)
## Singluar matrix encountered with more than 13 pressures here
data.keywest.pred <- data.frame(msl=keywest.tot.p.pred,
                                  t=seq(from=keywest.start.year.pred,to=keywest.end.year.pred),
                               keywest.sorted.pa.pred[,1:13])
tg.lmRob.keywest.pred <- lmRob(msl ~ ., x=T, y=T, data=data.keywest.pred, control=lmRobControl,
                              na.action=na.exclude)
## Variance reduction
tg.lmRob.keywest.pred$r.sq
#[1] 0.6592839

## Full
keywest.sorted.pa.full <- sort.corr.pa(keywest.corr.full[,,1],
                                      slpNAtlAnnArray,
                                      keywest.start.year.full,
                                      keywest.end.year.full, slp.yrs)
## Singluar matrix encountered with more than 13 pressures here
data.keywest.full <- data.frame(msl=keywest.tot.p.full,
                                  t=seq(from=keywest.start.year.full,to=keywest.end.year.full),
                               keywest.sorted.pa.full[,1:13])
tg.lmRob.keywest.full <- lmRob(msl ~ ., x=T, y=T, data=data.keywest.full, control=lmRobControl,
                              na.action=na.exclude)
## Variance reduction
tg.lmRob.keywest.full$r.sq
#[1] 0.7420201

##
## Most correlated
##
keywest.corr.pa.partial <- slpNAtlAnnArray[43,18,83:138]
data.keywest.corr.partial <- data.frame(msl=keywest.tot.p.partial,
                                  t=seq(from=keywest.start.year.partial,to=keywest.end.year.partial),
                                  keywest.corr.pa.partial)
tg.lmRob.keywest.corr.partial <- lmRob(msl ~ ., x=T, y=T, data=data.keywest.corr.partial,
                                 control=lmRobControl, na.action=na.exclude)
## Variance reduction
tg.lmRob.keywest.corr.partial$r.sq
#[1] 0.6592032

## Prediction
keywest.corr.pa.pred <- slpNAtlAnnArray[43,18,46:82]
data.keywest.corr.pred <- data.frame(msl=keywest.tot.p.pred,
                                  t=seq(from=keywest.start.year.pred,to=keywest.end.year.pred),
                               keywest.corr.pa.pred)
tg.lmRob.keywest.corr.pred <- lmRob(msl ~ ., x=T, y=T, data=data.keywest.corr.pred, control=lmRobControl,
                              na.action=na.exclude)
## Variance reduction
tg.lmRob.keywest.corr.pred$r.sq
#[1] 0.7157227

##*******************************************************************************************************
##*******************************************************************************************************
## Sanfran
sanfran.start.year.partial <- 1953
sanfran.end.year.partial <- 2008
sanfran.start.year.pred <- 1916
sanfran.end.year.pred <- 1952
sanfran.start.year.full <- 1916
sanfran.end.year.full <- 2008
sanfran.start.year.early <- 1854
sanfran.end.year.early <- 2008
sanfran.lon <- -122.467
sanfran.lat <- 37.8

sanfran.sl.partial <- load.tg.data(sanfran.start.year.partial,sanfran.end.year.partial,"sanfran")
sanfran.pa.partial <-
  interp.local.pressure(sanfran.lon,sanfran.lat,slpNPacAnnArray,range(lon2[nPacLonLatIndex$lon]),
                        range(lat[nPacLonLatIndex$lat]),
                         sanfran.start.year.partial, sanfran.end.year.partial, slp.yrs)

sanfran.tot.p.partial <- total.pressure(sanfran.sl.partial$Height,sanfran.pa.partial)

sanfran.corr.partial <- corr.pressure(sanfran.tot.p.partial, slpNPacAnnArray,
                             sanfran.start.year.partial, sanfran.end.year.partial, slp.yrs)

sanfran.corr.dt.partial <- corr.pressure.dt(sanfran.tot.p.partial, slpNPacAnnArray,
                             sanfran.start.year.partial, sanfran.end.year.partial, slp.yrs)

## Prediction
sanfran.sl.pred <- load.tg.data(sanfran.start.year.pred,sanfran.end.year.pred,"sanfran")
sanfran.pa.pred <-
  interp.local.pressure(sanfran.lon,sanfran.lat,slpNPacAnnArray,range(lon2[nPacLonLatIndex$lon]),
                        range(lat[nPacLonLatIndex$lat]), sanfran.start.year.pred,
                        sanfran.end.year.pred, slp.yrs)

sanfran.tot.p.pred <- total.pressure(sanfran.sl.pred$Height,sanfran.pa.pred)

sanfran.corr.pred <- corr.pressure(sanfran.tot.p.pred, slpNPacAnnArray,
                             sanfran.start.year.pred, sanfran.end.year.pred, slp.yrs)

## Full
sanfran.sl.full <- load.tg.data(sanfran.start.year.full,sanfran.end.year.full,"sanfran")
sanfran.pa.full <-
  interp.local.pressure(sanfran.lon,sanfran.lat,slpNPacAnnArray,
                        range(lon2[nPacLonLatIndex$lon]),range(lat[nPacLonLatIndex$lat]),
                         sanfran.start.year.full, sanfran.end.year.full, slp.yrs)

sanfran.tot.p.full <- total.pressure(sanfran.sl.full$Height,sanfran.pa.full)

sanfran.corr.full <- corr.pressure(sanfran.tot.p.full, slpNPacAnnArray,
                             sanfran.start.year.full, sanfran.end.year.full, slp.yrs)

sanfran.corr.dt.full <- corr.pressure.dt(sanfran.tot.p.full, slpNPacAnnArray,
                             sanfran.start.year.full, sanfran.end.year.full, slp.yrs)

## Start of SLP records: 1871
sanfran.sl.1871 <- load.tg.data(1871,sanfran.end.year.pred,"sanfran")
sanfran.pa.1871 <-
  interp.local.pressure(sanfran.lon,sanfran.lat,slpNPacAnnArray,
                        range(lon2[nPacLonLatIndex$lon]),range(lat[nPacLonLatIndex$lat]),
                         1871, sanfran.end.year.pred, slp.yrs)

sanfran.sl.1871.full <- load.tg.data(1871,sanfran.end.year.full,"sanfran")
sanfran.pa.1871.full <-
  interp.local.pressure(sanfran.lon,sanfran.lat,slpNPacAnnArray,
                        range(lon2[nPacLonLatIndex$lon]),range(lat[nPacLonLatIndex$lat]),
                         1871, sanfran.end.year.full, slp.yrs)

sanfran.tot.p.1871 <- total.pressure(sanfran.sl.1871$Height,sanfran.pa.1871)
sanfran.tot.p.1871.full <- total.pressure(sanfran.sl.1871.full$Height,sanfran.pa.1871.full)

sanfran.corr.1871 <- corr.pressure(sanfran.tot.p.1871.full, slpNPacAnnArray,
                             1871, sanfran.end.year.full, slp.yrs)

sanfran.corr.dt.1871 <- corr.pressure.dt(sanfran.tot.p.1871, slpNPacAnnArray,
                             1871, sanfran.end.year.pred, slp.yrs)

##*******************************************************************************************************
## Sort correlated data

lmRobControl <- lmRob.control(mxr=100,mxf=100,trace=F)

## Partial
sanfran.sorted.pa.partial <- sort.corr.pa(sanfran.corr.partial[,,1], slpNPacAnnArray,
                                         sanfran.start.year.partial, sanfran.end.year.partial, slp.yrs)
data.sanfran.partial <- data.frame(msl=sanfran.tot.p.partial,
                                  t=seq(from=sanfran.start.year.partial,to=sanfran.end.year.partial),
                                  sanfran.sorted.pa.partial)
tg.lmRob.sanfran.partial <- lmRob(msl ~ ., x=T, y=T, data=data.sanfran.partial, control=lmRobControl,
                                 na.action=na.exclude)
## Variance reduction
tg.lmRob.sanfran.partial$r.sq
#[1] 0.5960196

## Prediction

sanfran.sorted.pa.pred <- sort.corr.pa(sanfran.corr.pred[,,1],
                                      slpNPacAnnArray,
                                      sanfran.start.year.pred,
                                      sanfran.end.year.pred, slp.yrs)
## Singluar matrix encountered with more than 13 pressures here
data.sanfran.pred <- data.frame(msl=sanfran.tot.p.pred,
                                  t=seq(from=sanfran.start.year.pred,to=sanfran.end.year.pred),
                               sanfran.sorted.pa.pred[,1:13])
tg.lmRob.sanfran.pred <- lmRob(msl ~ ., x=T, y=T, data=data.sanfran.pred, control=lmRobControl,
                              na.action=na.exclude)
## Variance reduction
tg.lmRob.sanfran.pred$r.sq
#[1] 0.677276

## Full
sanfran.sorted.pa.full <- sort.corr.pa(sanfran.corr.full[,,1],
                                      slpNPacAnnArray,
                                      sanfran.start.year.full,
                                      sanfran.end.year.full, slp.yrs)
## Singluar matrix encountered with more than 13 pressures here
data.sanfran.full <- data.frame(msl=sanfran.tot.p.full,
                                  t=seq(from=sanfran.start.year.full,to=sanfran.end.year.full),
                               sanfran.sorted.pa.full[,1:13])
tg.lmRob.sanfran.full <- lmRob(msl ~ ., x=T, y=T, data=data.sanfran.full, control=lmRobControl,
                              na.action=na.exclude)
## Variance reduction
tg.lmRob.sanfran.full$r.sq
#[1] 0.720166

## 1871 full
sanfran.sorted.pa.1871.full <- sort.corr.pa(sanfran.corr.1871[,,1],
                                      slpNPacAnnArray,
                                      1871,
                                      sanfran.end.year.full, slp.yrs)
## Singluar matrix encountered with more than 13 pressures here
data.sanfran.1871.full <- data.frame(msl=sanfran.tot.p.1871.full,
                                  t=seq(from=1871,to=sanfran.end.year.full),
                               sanfran.sorted.pa.1871.full[,1:13])
tg.lmRob.sanfran.1871.full <- lmRob(msl ~ ., x=T, y=T, data=data.sanfran.1871.full, control=lmRobControl,
                              na.action=na.exclude)
## Variance reduction
tg.lmRob.sanfran.1871.full$r.sq
#[1] 0.6216469

##
## Most correlated
##
sanfran.corr.pa.partial <- slpNPacAnnArray[33,40,83:138]
data.sanfran.corr.partial <- data.frame(msl=sanfran.tot.p.partial,
                                  t=seq(from=sanfran.start.year.partial,to=sanfran.end.year.partial),
                                  sanfran.corr.pa.partial)
tg.lmRob.sanfran.corr.partial <- lmRob(msl ~ ., x=T, y=T, data=data.sanfran.corr.partial,
                                 control=lmRobControl, na.action=na.exclude)
## Variance reduction
tg.lmRob.sanfran.corr.partial$r.sq
#[1] 0.5334922

## Prediction
sanfran.corr.pa.pred <- slpNPacAnnArray[33,40,46:82]
data.sanfran.corr.pred <- data.frame(msl=sanfran.tot.p.pred,
                                  t=seq(from=sanfran.start.year.pred,to=sanfran.end.year.pred),
                               sanfran.corr.pa.pred)
tg.lmRob.sanfran.corr.pred <- lmRob(msl ~ ., x=T, y=T, data=data.sanfran.corr.pred, control=lmRobControl,
                              na.action=na.exclude)
## Variance reduction
tg.lmRob.sanfran.corr.pred$r.sq
#[1] 0.4632224

## Prediction to 1871
sanfran.corr.pa.1871 <- slpNPacAnnArray[71,22,1:82]
data.sanfran.corr.1871 <- data.frame(msl=sanfran.tot.p.1871,
                                  t=seq(from=1871,to=sanfran.end.year.pred),
                               sanfran.corr.pa.1871)
tg.lmRob.sanfran.corr.1871 <- lmRob(msl ~ ., x=T, y=T, data=data.sanfran.corr.1871, control=lmRobControl,
                              na.action=na.exclude)
## Variance reduction
tg.lmRob.sanfran.corr.1871$r.sq
#[1] 0.2041718

##*******************************************************************************************************
##*******************************************************************************************************
## Plots
##
#x11()
#make.image5(lon2[nAtlLonLatIndex$lon],flip.matrix(lat[nAtlLonLatIndex$lat]),keywest.corr.partial[,,1])

#x11()
lon2Pac <- lon2[nPacLonLatIndex$lon]
lon2Pac[which(lon2Pac<0)] <- lon2Pac[which(lon2Pac<0)]+360
#make.image5(lon2Pac,flip.matrix(lat[nAtlLonLatIndex$lat]),sanfran.corr.partial[,,1])

#x11()
#plot(keywest.start.year.partial:keywest.end.year.partial,keywest.tot.p.partial,col='blue',type='l',
#     ylim=c(0,4000), xlim=c(keywest.start.year.pred, keywest.end.year.partial), ann=F)

#lines(tg.lmRob.keywest.partial$x[,2], tg.lmRob.keywest.partial$fitted,
#      col='cyan')
#lines(keywest.start.year.pred:keywest.end.year.pred,keywest.tot.p.pred,col='blue')
#lines(tg.lmRob.keywest.pred$x[,2], tg.lmRob.keywest.pred$fitted, col='orange')

#x11()
#make.image5(lon2[nAtlLonLatIndex$lon],flip.matrix(lat[nAtlLonLatIndex$lat]),keywest.corr.full[,,1])

#x11()
#plot(sanfran.start.year.partial:sanfran.end.year.partial,sanfran.tot.p.partial,col='blue',type='l',
#     ylim=c(0,4000), xlim=c(sanfran.start.year.pred, sanfran.end.year.partial), ann=F)
#
#lines(tg.lmRob.sanfran.partial$x[,2], tg.lmRob.sanfran.partial$fitted,
#      col='cyan')
#lines(sanfran.start.year.pred:sanfran.end.year.pred,sanfran.tot.p.pred,col='blue')
#lines(tg.lmRob.sanfran.pred$x[,2], tg.lmRob.sanfran.pred$fitted, col='orange')
#lines(tg.lmRob.sanfran.full$x[,2],tg.lmRob.sanfran.full$fitted,col='red')

x11()
make.image5(lon2[nAtlLonLatIndex$lon],flip.matrix(lat[nAtlLonLatIndex$lat]),keywest.corr.dt.partial[,,1])
#x11()
#make.image5(lon2[nAtlLonLatIndex$lon],flip.matrix(lat[nAtlLonLatIndex$lat]),keywest.corr.dt.full[,,1])

x11()
plot(keywest.start.year.partial:keywest.end.year.partial,keywest.tot.p.partial,col='blue',type='l',
     ylim=c(0,4000), xlim=c(keywest.start.year.pred, keywest.end.year.partial), ann=F)

lines(tg.lmRob.keywest.corr.partial$x[,2], tg.lmRob.keywest.corr.partial$fitted,
      col='cyan')
lines(keywest.start.year.pred:keywest.end.year.pred,keywest.tot.p.pred,col='blue')
lines(tg.lmRob.keywest.corr.pred$x[,2], tg.lmRob.keywest.corr.pred$fitted, col='orange')

x11()
make.image4(lon2Pac,flip.matrix(lat[nPacLonLatIndex$lat]),sanfran.corr.dt.partial[,,1])

#x11()
#make.image4(lon2Pac,flip.matrix(lat[nPacLonLatIndex$lat]),sanfran.corr.dt.full[,,1])

#x11()
#make.image4(lon2Pac,flip.matrix(lat[nPacLonLatIndex$lat]),sanfran.corr.dt.1871[,,1])

x11()
plot(sanfran.start.year.partial:sanfran.end.year.partial,sanfran.tot.p.partial,col='blue',type='l',
     ylim=c(0,4000), xlim=c(sanfran.start.year.pred, sanfran.end.year.partial), ann=F)

lines(tg.lmRob.sanfran.corr.partial$x[,2], tg.lmRob.sanfran.corr.partial$fitted,
      col='cyan')
lines(sanfran.start.year.pred:sanfran.end.year.pred,sanfran.tot.p.pred,col='blue')
lines(tg.lmRob.sanfran.corr.pred$x[,2], tg.lmRob.sanfran.corr.pred$fitted, col='orange')

#x11()
#plot(sanfran.start.year.partial:sanfran.end.year.partial,sanfran.tot.p.partial,col='blue',type='l',
#     ylim=c(0,4000), xlim=c(sanfran.start.year.pred, sanfran.end.year.partial), ann=F)
#
#lines(tg.lmRob.sanfran.corr.partial$x[,2], tg.lmRob.sanfran.corr.partial.neg$fitted,
#      col='cyan')
#lines(sanfran.start.year.pred:sanfran.end.year.pred,sanfran.tot.p.pred,col='blue')
#lines(tg.lmRob.sanfran.corr.pred$x[,2], tg.lmRob.sanfran.corr.pred.neg$fitted, col='orange')

#x11()
#plot(sanfran.start.year.partial:sanfran.end.year.partial,sanfran.tot.p.partial,col='blue',type='l',
#     ylim=c(0,4000), xlim=c(1871, sanfran.end.year.partial), ann=F)
#
#lines(tg.lmRob.sanfran.corr.partial$x[,2], tg.lmRob.sanfran.corr.partial$fitted,
#      col='cyan')
#lines(1871:sanfran.end.year.pred,sanfran.tot.p.1871,col='blue')
#lines(tg.lmRob.sanfran.corr.1871$x[,2], tg.lmRob.sanfran.corr.1871$fitted, col='orange')

x11()
plot(sanfran.start.year.partial:sanfran.end.year.partial,sanfran.tot.p.partial,col='blue',type='l',
     ylim=c(0,4000), xlim=c(1871, sanfran.end.year.partial), ann=F)

lines(tg.lmRob.sanfran.corr.partial$x[,2], tg.lmRob.sanfran.corr.partial$fitted,
      col='cyan')
lines(1871:sanfran.end.year.pred,sanfran.tot.p.1871,col='blue')
lines(tg.lmRob.sanfran.corr.1871$x[,2], tg.lmRob.sanfran.corr.1871$fitted, col='orange')
