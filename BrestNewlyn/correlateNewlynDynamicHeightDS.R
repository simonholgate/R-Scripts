## Function to correlate sea level corrected for pressure at Newlyn (total
## pressure) with dynamic heights (based on Doug Smith - DS - and EN3c - en - data sets) everywhere


#############################
## Functions for use below ##
#############################


load.tg.data <- function(tg.start.year, tg.end.year, station){
  tg <- new.env()
  load("../tg/tg.RData", envir=tg)

  if(station=='brest'){
    tgs <- which(tg$brest$Year==tg.start.year) # index of start of Brest TG data
    tge <- which(tg$brest$Year==tg.end.year) # index of end of Brest TG data
    sea.level <- tg$brest[tgs:tge,]
  } else {
    if(station=='brest22'){
      brest22 <- tg$brest
      tgs <- which(tg$brest$Year==tg.start.year) # index of start of
                                                 # Brest TG data
      tge <- which(tg$brest$Year==tg.end.year) # index of end of Brest TG data

      ## 'Correct' pre-1943 data by 22mm
      if(tg.start.year<=1943){
        tgl <- which(brest22$Year==1943)
        brest22$Height[1:tgl] <- brest22$Height[1:tgl]-22
      }
      
      sea.level <- brest22[tgs:tge,]
      
    } else {
      if(station=='newlyn'){
        tgs <- which(tg$newlyn$Year==tg.start.year) # index of start of
                                                    # Newlyn TG data
        tge <- which(tg$newlyn$Year==tg.end.year) # index of end of
                                                  # Newlyn TG data
        sea.level <- tg$newlyn[tgs:tge,]
      } else {
        error(paste("station not known:", station))
      }
    }
  }

  sea.level
}

##*****************************************************************************

total.pressure <- function(sea.level, local.pressure){
  total.p <- 1025*9.8*sea.level/1000 + local.pressure

  total.p
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

corr.dh <- function(tot.p, dhArray, start.year, end.year, years){
  
  nlon <- dim(dhArray)[1]
  nlat <- dim(dhArray)[2]
  corr.array <- array(NA, dim=c(nlon,nlat,4))
  
  start <- which(years==start.year)
  end <- which(years==end.year)
  
  for(i in 1:nlon){
    for(j in 1:nlat){
       corr <- cor.test(tot.p, dhArray[i,j,start:end], alternative="two.sided",
        method="pearson", na.action=na.exclude)
       if(corr$p.value<=0.05){ # We're only interested in >= 95% confidence
         corr.array[i,j,] <-
          c(corr$estimate,corr$conf.int[1],corr$conf.int[2],corr$p.value)
       }
    }
  }

  corr.array
}

##*****************************************************************************
## Robustly linearly detrend a vector. Missing values must be included as NA,
## but are excluded in the resgression.
detrend <- function(data_vector){
  df.data_vector <- data.frame(p=data_vector,t=seq(from=1,
                                               to=length(data_vector)))
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
corr.dh.dt <- function(tot.p, dhArray, start.year, end.year, years){
  
  nlon <- dim(dhArray)[1]
  nlat <- dim(dhArray)[2]
  corr.array <- array(NA, dim=c(nlon,nlat,4))
  
  start <- which(years==start.year)
  end <- which(years==end.year)
  
  for(i in 1:nlon){
    for(j in 1:nlat){
       corr <- cor.test(tot.p, detrend(dhArray[i,j,start:end]), alternative="two.sided",
        method="pearson", na.action=na.exclude)
       if(corr$p.value<=0.05){ # We're only interested in >= 95% confidence
         corr.array[i,j,] <-
          c(corr$estimate,corr$conf.int[1],corr$conf.int[2],corr$p.value)
       }
    }
  }

  corr.array
}

##*****************************************************************************
corr.dh.dt2 <- function(tot.p, p.yr.range, dhArray, start.year, end.year, years){
  
  nlon <- dim(dhArray)[1]
  nlat <- dim(dhArray)[2]
  corr.array <- array(NA, dim=c(nlon,nlat,4))

  p.yrs <- c(p.yr.range[1]:p.yr.range[2])
  
  pstart <- which(p.yrs==start.year)
  pend <- which(p.yrs==end.year)
  
  start <- which(years==start.year)
  end <- which(years==end.year)
  
  for(i in 1:nlon){
    for(j in 1:nlat){
      ## Land is all NaNs/NAs
      if(!all(is.nan(dhArray[i,j,start:end])) && !all(is.na(dhArray[i,j,start:end]))){ 
        corr <- cor.test(tot.p[pstart:pend], detrend(dhArray[i,j,start:end]), alternative="two.sided",
                         method="pearson", na.action=na.exclude)
       
        if(corr$p.value<=0.05){ # We're only interested in >= 95% confidence
          corr.array[i,j,] <-
            c(corr$estimate,corr$conf.int[1],corr$conf.int[2],corr$p.value)
        }
      }
    }
  }

  corr.array
}

##*****************************************************************************

make.image <- function(lon,lat,data.array){
  image.plot(lon,lat,data.array, zlim=c(-0.6,0.6))
  map("worldHires", xlim=c(-100,15), ylim=c(-5,80), interior=F, fill=F,
      col="grey50", resolution=0, add=T)
  map.axes()
  points(expand.grid(c(-20, -15,-10, -5, 0 , 5, 10), c(40,45,50,55,60)),
         pch=19, col="grey20")
  points( y=50.1, x=-5.55, pch='*', col="grey20", cex=3)
  points( y=48.38, x=-4.5, pch='*', col="grey20", cex=3)
  points(y=c(45, 45, 45, 45, 40, 40, 40, 50, 50), x=c(-10, -15, -20,
  -5, -15, -10, -20, -20, -15), pch=19, col="yellow")
}

##*****************************************************************************

make.image2 <- function(lon,lat,data.array){
  image(lon,lat,data.array, zlim=c(-0.6,0.6), col=tim.colors())
  map("worldHires", xlim=c(-100,15), ylim=c(-5,80), interior=F, fill=F,
      col="grey50", resolution=0, add=T)
  map.axes()
  points(expand.grid(c(-20, -15,-10, -5, 0 , 5, 10), c(40,45,50,55,60)),
         pch=19, col="grey20")
  points( y=50.1, x=-5.55, pch='*', col="grey20", cex=3)
  points( y=48.38, x=-4.5, pch='*', col="grey20", cex=3)
  points(y=c(45, 45, 45, 45, 40, 40, 40, 50, 50), x=c(-10, -15, -20,
  -5, -15, -10, -20, -20, -15), pch=19, col="yellow")
}



##*****************************************************************************

make.image3 <- function(lon,lat,data.array){
  image.plot(lon,lat,data.array, zlim=c(-0.6,0.6))
  world(add=T)
##  points(expand.grid(c(-20, -15,-10, -5, 0 , 5, 10), c(40,45,50,55,60)),
##         pch=19, col="grey20")
  points( y=50.1, x=-5.55, pch='*', col="grey20", cex=3)
  points( y=48.38, x=-4.5, pch='*', col="grey20", cex=3)
##  points(-25.7,35.75, pch=19, col="magenta", cex=1)
##  points(-7.5, 40.9517, pch=19, col="grey20", cex=1)
}

##*****************************************************************************
# As make.image3 but without a fixed zlim
make.image4 <- function(lon,lat,data.array){
  image.plot(lon,lat,data.array)
  world(add=T, shift=F)
##  points(expand.grid(c(-20, -15,-10, -5, 0 , 5, 10), c(40,45,50,55,60)),
##         pch=19, col="grey20")
  points( y=50.1, x=-5.55, pch='*', col="grey20", cex=3)
  points( y=48.38, x=-4.5, pch='*', col="grey20", cex=3)
##  points( y=37.8, x=-122.467, pch='*', col="grey20", cex=3)
##  points(y=c(45, 45, 45, 45, 40, 40, 40, 50, 50), x=c(-10, -15, -20,
##  -5, -15, -10, -20, -20, -15), pch=19, col="yellow")
}

##*****************************************************************************
# Using filled contours
make.image5 <- function(lon,lat,data.array, arrow.array.x, arrow.array.y){
  lenLon <- length(lon)
  lenLat <- length(lat)
  lonLat <- expand.grid(lon[seq(from=2, to=lenLon, by=4)],
                      lat[seq(from=2, to=lenLat, by=4)])
  filled.contour(lon, lat, data.array, color=tim.colors, nlevels=20,
               plot.axes = { world(add=T);
                             ##points(-7.5, 40.9517, pch=19, col="grey60", cex=1);
                             points(-25.7,35.75, pch=19, col="magenta", cex=1);
                             grid(col="black");
                             arrow.plot(lonLat[,1], lonLat[,2],
                                        u=as.vector(arrow.array.x[seq(from=2, to=lenLon, by=4),
                                          seq(from=2, to=lenLat, by=4)]),
                                        v=as.vector(arrow.array.y[seq(from=2, to=lenLon, by=4),
                                          seq(from=2, to=lenLat, by=4)]),
                                        arrow.ex=.2, col='gray50', length=.05, lwd=2,
                                        true.angle=T);
                             axis(1, at=NULL);
                             axis(2, at=NULL) },
               zlim=c(-0.5,0.5))
}

##*****************************************************************************
# Using filled contours but no wind arrows and no zlim
make.image6 <- function(lon,lat,data.array){
  lenLon <- length(lon)
  lenLat <- length(lat)
  lonLat <- expand.grid(lon[seq(from=2, to=lenLon, by=4)],
                      lat[seq(from=2, to=lenLat, by=4)])
  filled.contour(lon, lat, data.array, color=tim.colors, nlevels=20,
               plot.axes = { world(add=T);
                             ##points(-7.5, 40.9517, pch=19, col="grey60", cex=1);
                             points(-25.7,35.75, pch=19, col="magenta", cex=1);
                             grid(col="black");
                             axis(1, at=NULL);
                             axis(2, at=NULL) },
               zlim=c(-0.9,0.9))
}

##*****************************************************************************
sort.corr.dh <- function(corr.array, dh.array, start.year, end.year, years){

  ps <- which(years==start.year)
  pe <- which(years==end.year)
  
  sorted.corr <- sort(abs(corr.array),decreasing=T, index=T)
  
  tmp <- dh.array
  i <- dim(tmp)[1]
  j <- dim(tmp)[2]
  k <- dim(tmp)[3]
  
  dim(tmp) <- c((i*j),k)
  
  sorted.pa <- tmp[sorted.corr$ix[1:20],ps:pe] # Only return top 20 correlations

  t(sorted.pa)
}

##*****************************************************************************
most.corr.dh <- function(corr.array, dh.array, start.year, end.year, years){
  
  ps <- which(years==start.year)
  pe <- which(years==end.year)
  
  most.correlated.point <- which(abs(corr.array)==max(abs(corr.array),
                                      na.rm=T), arr.ind=T)
  
  most.correlated <- dh.array[most.correlated.point[1],
                               most.correlated.point[2], ps:pe] 

  most.correlated
  
}

##*****************************************************************************

#########################
## Non-functional part ##
#########################

library(fields)
library(robust)
source("~/Dropbox/BrestNewlyn/matrixMethods.R")
load("~/Dropbox/brestNewlynData/analysis/paper/correlationACRE/brestNewlyn.tot.ps.RData")

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
## N Atlantic dynamic height

load("~/Dropbox/DynamicHeight/dh_DS_EN3c_IS.RData")


##*****************************************************************************
## Newlyn

lmRobControl <- lmRob.control(mxr=100,mxf=100,trace=F)

## Partial
newlyn.corr.dt.partial.ds <- corr.dh.dt2(newlyn.tot.p.full, c(full.p.start, full.p.end), dhdsAnnualMean,
                             newlyn.start.year.partial.ds,
                                     newlyn.end.year.partial.ds, dhds.yrs)

newlyn.corr.dt.partial.en <- corr.dh.dt2(newlyn.tot.p.full, c(full.p.start, full.p.end), dhenAnnualMean,
                             newlyn.start.year.partial.en,
                                     newlyn.end.year.partial.en, dhen.yrs)

newlyn.corr.dt.partial.is <- corr.dh.dt2(newlyn.tot.p.full, c(full.p.start, full.p.end), dhisAnnualMean,
                             newlyn.start.year.partial.is,
                                     newlyn.end.year.partial.is, dhis.yrs)


## Prediction
newlyn.corr.dt.pred.ds <- corr.dh.dt2(newlyn.tot.p.full, c(full.p.start, full.p.end), dhdsAnnualMean,
                             newlyn.start.year.pred.ds, newlyn.end.year.pred.ds,
                                  dhds.yrs)
newlyn.corr.dt.pred.en <- corr.dh.dt2(newlyn.tot.p.full, c(full.p.start, full.p.end), dhenAnnualMean,
                             newlyn.start.year.pred.en, newlyn.end.year.pred.en,
                                  dhen.yrs)
newlyn.corr.dt.pred.is <- corr.dh.dt2(newlyn.tot.p.full, c(full.p.start, full.p.end), dhisAnnualMean,
                             newlyn.start.year.pred.is, newlyn.end.year.pred.is,
                                  dhis.yrs)

## Full
newlyn.corr.dt.full.ds <- corr.dh.dt2(newlyn.tot.p.full, c(full.p.start, full.p.end), dhdsAnnualMean,
                             newlyn.start.year.full.ds, newlyn.end.year.full.ds,
                                        dhds.yrs)
newlyn.corr.dt.full.en <- corr.dh.dt2(newlyn.tot.p.full, c(full.p.start, full.p.end), dhenAnnualMean,
                             newlyn.start.year.full.en, newlyn.end.year.full.en,
                                        dhen.yrs)
newlyn.corr.dt.full.is <- corr.dh.dt2(newlyn.tot.p.full, c(full.p.start, full.p.end), dhisAnnualMean,
                             newlyn.start.year.full.is, newlyn.end.year.full.is,
                                        dhis.yrs)

## Consistent period - last 20 years
newlyn.corr.dt.cons.ds <- corr.dh.dt2(newlyn.tot.p.full, c(full.p.start, full.p.end), dhdsAnnualMean,
                             newlyn.start.year.cons.ds, newlyn.end.year.cons.ds,
                                        dhds.yrs)
newlyn.corr.dt.cons.en <- corr.dh.dt2(newlyn.tot.p.full, c(full.p.start, full.p.end), dhenAnnualMean,
                             newlyn.start.year.cons.en, newlyn.end.year.cons.en,
                                        dhen.yrs)
newlyn.corr.dt.cons.is <- corr.dh.dt2(newlyn.tot.p.full, c(full.p.start, full.p.end), dhisAnnualMean,
                             newlyn.start.year.cons.is, newlyn.end.year.cons.is,
                                        dhis.yrs)

##
## Plots
##
x11()
make.image6(londs,latds,newlyn.corr.dt.full.ds[,,1])
x11()
make.image6(londs,latds,newlyn.corr.dt.partial.ds[,,1])
x11()
make.image6(londs,latds,newlyn.corr.dt.pred.ds[,,1])
x11()
make.image6(londs,latds,newlyn.corr.dt.cons.ds[,,1])
x11()
make.image6(lonen,laten,newlyn.corr.dt.full.en[,,1])
x11()
make.image6(lonen,laten,newlyn.corr.dt.partial.en[,,1])
x11()
make.image6(lonen,laten,newlyn.corr.dt.pred.en[,,1])
x11()
make.image6(lonen,laten,newlyn.corr.dt.cons.en[,,1])
x11()
make.image6(lonis,latis,newlyn.corr.dt.full.is[,,1])
x11()
make.image6(lonis,latis,newlyn.corr.dt.partial.is[,,1])
x11()
make.image6(lonis,latis,newlyn.corr.dt.pred.is[,,1])
x11()
make.image6(lonis,latis,newlyn.corr.dt.cons.is[,,1])

##
## Most correlated
##

## Partial
newlyn.partial.most.corr.ds <- most.corr.dh(newlyn.corr.dt.full.ds[,,1],
                                          dhdsAnnualMean,
                                          newlyn.start.year.partial.ds,
                                          newlyn.end.year.partial.ds, dhds.yrs)

newlyn.partial.most.corr.en <- most.corr.dh(newlyn.corr.dt.full.en[,,1],
                                          dhenAnnualMean,
                                          newlyn.start.year.partial.en,
                                          newlyn.end.year.partial.en, dhen.yrs)

newlyn.partial.most.corr.is <- most.corr.dh(newlyn.corr.dt.full.is[,,1],
                                          dhisAnnualMean,
                                          newlyn.start.year.partial.is,
                                          newlyn.end.year.partial.is, dhis.yrs)


data.newlyn.corr.partial.ds <- data.frame(msl=newlyn.tot.p.full[full.p.start.partial.ds:full.p.end.partial.ds],
                                  t=seq(from=newlyn.start.year.partial.ds,
                                    to=newlyn.end.year.partial.ds),
                                  newlyn.partial.most.corr.ds)
tg.lmRob.newlyn.corr.partial.ds <- lmRob(msl ~ ., x=T, y=T,
                                      data=data.newlyn.corr.partial.ds,
                                 control=lmRobControl, na.action=na.exclude)
## Variance reduction
tg.lmRob.newlyn.corr.partial.ds$r.sq
#[1] 0.4636512

data.newlyn.corr.partial.en <- data.frame(msl=newlyn.tot.p.full[full.p.start.partial.en:full.p.end.partial.en],
                                  t=seq(from=newlyn.start.year.partial.en,
                                    to=newlyn.end.year.partial.en),
                                  newlyn.partial.most.corr.en)
tg.lmRob.newlyn.corr.partial.en <- lmRob(msl ~ ., x=T, y=T,
                                      data=data.newlyn.corr.partial.en,
                                 control=lmRobControl, na.action=na.exclude)
## Variance reduction
tg.lmRob.newlyn.corr.partial.en$r.sq
#[1] 0.2702787

data.newlyn.corr.partial.is <- data.frame(msl=newlyn.tot.p.full[full.p.start.partial.is:full.p.end.partial.is],
                                  t=seq(from=newlyn.start.year.partial.is,
                                    to=newlyn.end.year.partial.is),
                                  newlyn.partial.most.corr.is)
tg.lmRob.newlyn.corr.partial.is <- lmRob(msl ~ ., x=T, y=T,
                                      data=data.newlyn.corr.partial.is,
                                 control=lmRobControl, na.action=na.exclude)
## Variance reduction
tg.lmRob.newlyn.corr.partial.is$r.sq
#[1] 0.4636512


## Full

newlyn.full.most.corr.ds <- most.corr.dh(newlyn.corr.dt.full.ds[,,1],
                                          dhdsAnnualMean,
                                          newlyn.start.year.full.ds,
                                          newlyn.end.year.full.ds, dhds.yrs)

newlyn.full.most.corr.en <- most.corr.dh(newlyn.corr.dt.full.en[,,1],
                                          dhenAnnualMean,
                                          newlyn.start.year.full.en,
                                          newlyn.end.year.full.en, dhen.yrs)

newlyn.full.most.corr.is <- most.corr.dh(newlyn.corr.dt.full.is[,,1],
                                          dhisAnnualMean,
                                          newlyn.start.year.full.is,
                                          newlyn.end.year.full.is, dhis.yrs)


data.newlyn.corr.full.ds <- data.frame(msl=newlyn.tot.p.full[full.p.start.full.ds:full.p.end.full.ds],
                                  t=seq(from=newlyn.start.year.full.ds,
                                    to=newlyn.end.year.full.ds),
                                  newlyn.full.most.corr.ds)
tg.lmRob.newlyn.corr.full.ds <- lmRob(msl ~ ., x=T, y=T,
                                      data=data.newlyn.corr.full.ds,
                                 control=lmRobControl, na.action=na.exclude)
## Variance reduction
tg.lmRob.newlyn.corr.full.ds$r.sq
#[1] 0.46601

data.newlyn.corr.full.en <- data.frame(msl=newlyn.tot.p.full[full.p.start.full.en:full.p.end.full.en],
                                  t=seq(from=newlyn.start.year.full.en,
                                    to=newlyn.end.year.full.en),
                                  newlyn.full.most.corr.en)
tg.lmRob.newlyn.corr.full.en <- lmRob(msl ~ ., x=T, y=T,
                                      data=data.newlyn.corr.full.en,
                                 control=lmRobControl, na.action=na.exclude)
## Variance reduction
tg.lmRob.newlyn.corr.full.en$r.sq
#[1] 0.4924827

data.newlyn.corr.full.is <- data.frame(msl=newlyn.tot.p.full[full.p.start.full.is:full.p.end.full.is],
                                  t=seq(from=newlyn.start.year.full.is,
                                    to=newlyn.end.year.full.is),
                                  newlyn.full.most.corr.is)
tg.lmRob.newlyn.corr.full.is <- lmRob(msl ~ ., x=T, y=T,
                                      data=data.newlyn.corr.full.is,
                                 control=lmRobControl, na.action=na.exclude)
## Variance reduction
tg.lmRob.newlyn.corr.full.is$r.sq
#[1] 0.4910811

## Prediction
newlyn.pred.most.corr.ds <- most.corr.dh(newlyn.corr.dt.pred.ds[,,1],
                                          dhdsAnnualMean,
                                          newlyn.start.year.pred.ds,
                                          newlyn.end.year.pred.ds, dhds.yrs)

newlyn.pred.most.corr.en <- most.corr.dh(newlyn.corr.dt.pred.en[,,1],
                                          dhenAnnualMean,
                                          newlyn.start.year.pred.en,
                                          newlyn.end.year.pred.en, dhen.yrs)

newlyn.pred.most.corr.is <- most.corr.dh(newlyn.corr.dt.pred.is[,,1],
                                          dhisAnnualMean,
                                          newlyn.start.year.pred.is,
                                          newlyn.end.year.pred.is, dhis.yrs)


data.newlyn.corr.pred.ds <- data.frame(msl=newlyn.tot.p.full[full.p.start.pred.ds:full.p.end.pred.ds],
                                  t=seq(from=newlyn.start.year.pred.ds,
                                    to=newlyn.end.year.pred.ds),
                                  newlyn.pred.most.corr.ds)
tg.lmRob.newlyn.corr.pred.ds <- lmRob(msl ~ ., x=T, y=T,
                                      data=data.newlyn.corr.pred.ds,
                                 control=lmRobControl, na.action=na.exclude)
## Variance reduction
tg.lmRob.newlyn.corr.pred.ds$r.sq
#[1] 0.5941229

data.newlyn.corr.pred.en <- data.frame(msl=newlyn.tot.p.full[full.p.start.pred.en:full.p.end.pred.en],
                                  t=seq(from=newlyn.start.year.pred.en,
                                    to=newlyn.end.year.pred.en),
                                  newlyn.pred.most.corr.en)
tg.lmRob.newlyn.corr.pred.en <- lmRob(msl ~ ., x=T, y=T,
                                      data=data.newlyn.corr.pred.en,
                                 control=lmRobControl, na.action=na.exclude)
## Variance reduction
tg.lmRob.newlyn.corr.pred.en$r.sq
#[1] 0.6844057

data.newlyn.corr.pred.is <- data.frame(msl=newlyn.tot.p.full[full.p.start.pred.is:full.p.end.pred.is],
                                  t=seq(from=newlyn.start.year.pred.is,
                                    to=newlyn.end.year.pred.is),
                                  newlyn.pred.most.corr.is)
tg.lmRob.newlyn.corr.pred.is <- lmRob(msl ~ ., x=T, y=T,
                                      data=data.newlyn.corr.pred.is,
                                 control=lmRobControl, na.action=na.exclude)
## Variance reduction
tg.lmRob.newlyn.corr.pred.is$r.sq
#[1] 0.3304293
