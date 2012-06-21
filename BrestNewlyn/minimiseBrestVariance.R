## Function to minimise the variance between the reconstructed sea level,
## corrected for pressure, at Brest (total pressure) when an offset is added
## to account for the jump across the wars


#############################
## Functions for use below ##
#############################

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

##*****************************************************************************

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

##*****************************************************************************

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

##*****************************************************************************

natl.lon.lat <- function(xlon, ylat){
#  xlon<-seq(from=-180,by=5,length=nlon)
#  ylat<-seq(from=90,by=-5,length=nlat)
  dlon <- xlon[2]-xlon[1]
  dlat <- ylat[2]-ylat[1]

  # Atlantic covers 100W to 15E and 5S to 80N
  # We need to see whether we are working in 0 to 360 or -180 to 180
  # co-ordinate system
  if(min(xlon)<0){
    nAtlLon <- match(seq(from=-100,to=15,by=dlon), xlon)
  } else {
    nAtlLon <- c(match(seq(from=260,to=359,by=dlon), xlon),
                 match(seq(from=0,to=15,by=dlon), xlon))
  }
  # Do we go from N to S or S to N?
  if(dlat>0){ # S to N
    nAtlLat <- match(seq(from=-5,to=80, by=dlat), ylat)
  } else { # N to S
    nAtlLat <- match(seq(from=80,to=-5, by=dlat), ylat)
  }
  
  list(lon=nAtlLon,lat=nAtlLat)
}

##*****************************************************************************

load.tgs <- function(){
  tg <- new.env()
  load("../tg/tg.RData", envir=tg)
  tg
}

##*****************************************************************************

total.pressure <- function(sea.level, local.pressure){
  total.p <- 1025*9.8*sea.level/1000 + local.pressure

  total.p
}

##*****************************************************************************

offset.brest.tg.data <- function(tg.start.year, tg.end.year, offset){
##  tg <- new.env()
##  load("../tg/tg.RData", envir=tg)

     brestx <- tg$brest
      tgs <- which(tg$brest$Year==tg.start.year) # index of start of
                                                 # Brest TG data
      tge <- which(tg$brest$Year==tg.end.year) # index of end of Brest TG data

      ## 'Correct' pre-1943 data by 'offset' mm
      if(tg.start.year<=1943){
        tgl <- which(brestx$Year==1943)
        brestx$Height[1:tgl] <- brestx$Height[1:tgl] - offset
      }
      
      sea.level <- brestx[tgs:tge,]

  sea.level
}


##*****************************************************************************

extract.local.pressure <- function(lon,lat, slpArray, lon.range,
                                   lat.range, start.year, end.year, years){
  lon.values <- seq(from=lon.range[1], to=lon.range[2], by=5)
  lat.values <- seq(from=lat.range[1], to=lat.range[2], by=5)
  
  start <- which(years==start.year)
  end <- which(years==end.year)
  
  x <- which(lon.values == lon)
  y <- which(lat.values == lat)
  local.pressure <- slpArray[x,y,start:end]
  
  local.pressure

}

##*****************************************************************************

interp.local.pressure <- function(lon,lat, slpArray, lon.range,
                                   lat.range, start.year, end.year, years){
  lon.values <- seq(from=lon.range[1], to=lon.range[2], by=2)
  lat.values <- seq(from=lat.range[1], to=lat.range[2], by=2)
  
  start <- which(years==start.year)
  end <- which(years==end.year)
  nyrs <- end - start

  # Check whether lon runs from 0 to 360 or -180 to 180. We want it to be -180
  # to 180.

  if ((lon >= 180) && (lon.range[1]<0)){
    lon <- lon - 360
  } 

  # For each of the nyrs, we need to interpolate the grid to the locations
  # we've chosen
  local.pressure <- vector(length=nyrs, mode="numeric")
  for (i in start:end) {
    obj <- list(x=lon.values, y=lat.values, z=slpArray[,,i])
    local.pressure[i-start+1] <- interp.surface(obj,cbind(lon,lat))
  }

  local.pressure

}

##*****************************************************************************
## Robustly linearly detrend a vector. Missing values must be included as NA,
## but are excluded in the resgression.
detrend <- function(data_vector){
  df.data_vector <- data.frame(p=data_vector,
                               t=seq(from=1,to=length(data_vector)))
  data.lmRob <- lmRob(p ~ t, x=T, y=T, data=df.data_vector,
                      control=lmRobControl, na.action=na.exclude)
## We want the returned vector to be the same length as the original vector
## so use the $x values
  resid.data_vector <- vector(mode="numeric", length=length(data_vector))
  resid.data_vector[data.lmRob$x[,2]] <- data.lmRob$residuals 

  resid.data_vector
}

##*****************************************************************************
## Correlate detrended pressure
corr.pressure.dt <- function(tot.p, slpArray, start.year, end.year, years){
  
  nlon <- dim(slpArray)[1]
  nlat <- dim(slpArray)[2]
  corr.array <- array(NA, dim=c(nlon,nlat,4))
  
  start <- which(years==start.year)
  end <- which(years==end.year)
  
  for(i in 1:nlon){
    for(j in 1:nlat){
       corr <- cor.test(tot.p, detrend(slpArray[i,j,start:end]),
                        alternative="two.sided", method="pearson",
                        na.action=na.exclude)
       if(corr$p.value<=0.05){ # We're only interested in >= 95% confidence
         corr.array[i,j,] <-
          c(corr$estimate,corr$conf.int[1],corr$conf.int[2],corr$p.value)
       }
    }
  }

  corr.array
}

##*****************************************************************************
sort.corr.pa <- function(corr.array, slp.array, start.year, end.year, years){

  ps <- which(years==start.year)
  pe <- which(years==end.year)
  
  sorted.corr <- sort(abs(corr.array),decreasing=T, index=T)
  
  tmp <- slp.array
  i <- dim(tmp)[1]
  j <- dim(tmp)[2]
  k <- dim(tmp)[3]
  
  dim(tmp) <- c((i*j),k)
  
  sorted.pa <- tmp[sorted.corr$ix[1:20],ps:pe] # Only return top 20 correlations

  t(sorted.pa)
}

##*****************************************************************************
most.corr.pa <- function(corr.array, slp.array, start.year, end.year, years){
  
  ps <- which(years==start.year)
  pe <- which(years==end.year)
  
  most.correlated.point <- which(abs(corr.array)==max(abs(corr.array), na.rm=T))
  
  tmp <- slp.array
  i <- dim(tmp)[1]
  j <- dim(tmp)[2]
  k <- dim(tmp)[3]
  
  dim(tmp) <- c((i*j),k)
  
  most.correlated <- tmp[most.correlated.point,ps:pe] # Only return top 20 correlations

  most.correlated

}

##*****************************************************************************

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
slp.yrs <- c(1871:2009)

##*****************************************************************************
## N Atlantic pressure
##*****************************************************************************

## Read in monthly sea level pressures
slpMonthlyArray <- monthly.slp.ncdf("~/data/ACRE/prmsl.mon.mean.nc",
                                    nlon, nlat, nyr)
## Leave pressure as Pa

## Extract N Atlantic region from the SLP
nAtlLonLatIndex <- natl.lon.lat(lon,lat)
slpNAtlArray <- natl.slp(slpMonthlyArray, nAtlLonLatIndex$lon,
                         nAtlLonLatIndex$lat)

## Convert monthly SLP to annual SLP
slpNAtlAnnArray <- annual.slp(slpNAtlArray)

## Set the robust parameters
lmRobControl <- lmRob.control(mxr=100,mxf=100,trace=F)

##*****************************************************************************
## Brestx
##*****************************************************************************
brestx.start.year <- 1916
brestx.end.year <- 2008
offset.low <- 35
offset.high <- 50
num.offsets <- length(offset.low:offset.high)

var.red <- vector(mode="numeric", length=num.offsets)

tg <- load.tgs()

## Get the slp to go with the sl data
brestx.pa <-
  interp.local.pressure(-4.5,48.38,slpNAtlAnnArray,
                        range(lon2[nAtlLonLatIndex$lon]),
                        range(lat[nAtlLonLatIndex$lat]),
                        brestx.start.year,
                        brestx.end.year, slp.yrs)

##*****************************************************************************
## We're going to look at offsets from 18mm to 25mm and see which reduces the
## variance most. Douglas took 22mm and we want to see if that fits us best.
##*****************************************************************************

count <- 1

for (offset in offset.low:offset.high) {
  ## Progress indicator
  message(offset)
  ## Get the sl data over the entire period 1971-2008
  brestx.sl <- offset.brest.tg.data(brestx.start.year, brestx.end.year, offset)
  ## Convert sl into total pressure
  brestx.tot.p <- total.pressure(brestx.sl$Height, brestx.pa)
  
  ## Correlate the detrended total pressure with detrended
  ## slp everywhere.
  ## Returns a 3D array of dimensions (lon,lat,4) with correlation estimate,
  ## conf int [1], conf int [2] and p value. Any estimates with a p value
  ## giving less than 95% confidence is returned as NA
  brestx.corr.dt <- corr.pressure.dt(brestx.tot.p,
                                            slpNAtlAnnArray,
                                            brestx.start.year,
                                            brestx.end.year, slp.yrs)
  ## Sort the correlated data
  brestx.most.corr.pa <- most.corr.pa(brestx.corr.dt[,,1],
                                          slpNAtlAnnArray,
                                          brestx.start.year,
                                          brestx.end.year, slp.yrs)
  ## Produce a data frame of the most correlated points
  data.brestx <- data.frame(msl=brestx.tot.p,
                                  t=seq(from=brestx.start.year,
                                    to=brestx.end.year),
                                  brestx.most.corr.pa)
  ## Build the model
  tg.lmRob.brestx <- lmRob(msl ~ ., x=T, y=T,
                                  data=data.brestx,
                                  control=lmRobControl,
                                  na.action=na.exclude)
  ## Variance reduction
  var.red[count] <- (var(brestx.tot.p, na.rm=T)-var(tg.lmRob.brestx$resid, na.rm=T)) /
    var(brestx.tot.p, na.rm=T)*100
  count <- count + 1
 
} ## Loop over offsets

##*********##
##* Plots *##
##*********##
plot(1916:2008, brestx.tot.p, type='l', col='blue')
lines(tg.lmRob.brestx$x[,2],tg.lmRob.brestx$fitted , col='red')
