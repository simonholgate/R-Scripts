## Function to correlate sea level corrected for pressure at Newlyn (total pressure) with sea level
## pressure everywhere


#############################
## Functions for use below ##
#############################

monthly.slp <- function(path, nlon, nlat, nyr){
  ## Read the sea level data into an array
  nmon <- 12
  slpMonthlyArray<-array(NA,dim=c(nlon,nlat,nmon,nyr))

  con<-file(path,"rb")

  for (i in 1:nlat){
    for (j in 1:nlon){
      for (k in 1:nyr){
        slp<-readBin(con,what="numeric", size=4, n=nmon, endian='little')
        slpMonthlyArray[j,i,,k]<-slp
      }
    }
  }

  close(con)
  
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

natl.lon.lat <- function(){
  xlon<-seq(from=-180,by=5,length=nlon)
  ylat<-seq(from=90,by=-5,length=nlat)

  nAtlLon <- match(seq(from=-100,to=15,by=5), xlon)
  nAtlLat <- match(seq(from=-5,to=80, by=5), ylat)
  
  list(lon=nAtlLon,lat=nAtlLat)
}

##*******************************************************************************************************

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
      tgs <- which(tg$brest$Year==tg.start.year) # index of start of Brest TG data
      tge <- which(tg$brest$Year==tg.end.year) # index of end of Brest TG data

      ## 'Correct' pre-1943 data by 22mm
      if(tg.start.year<=1943){
        tgl <- which(brest22$Year==1943)
        brest22$Height[1:tgl] <- brest22$Height[1:tgl]-22
      }
      
      sea.level <- brest22[tgs:tge,]
      
    } else {
      if(station=='newlyn'){
        tgs <- which(tg$newlyn$Year==tg.start.year) # index of start of Newlyn TG data
        tge <- which(tg$newlyn$Year==tg.end.year) # index of end of Newlyn TG data
        sea.level <- tg$newlyn[tgs:tge,]
      } else {
        error(paste("station not known:", station))
      }
    }
  }

  sea.level
}

##*******************************************************************************************************

total.pressure <- function(sea.level, local.pressure){
  total.p <- 1025*9.8*sea.level/1000 + local.pressure

  total.p
}

##*******************************************************************************************************

extract.local.pressure <- function(lon,lat, slpArray, lon.range,
                                   lat.range, start.year, end.year){
  lon.values <- seq(from=lon.range[1], to=lon.range[2], by=5)
  lat.values <- seq(from=lat.range[1], to=lat.range[2], by=5)
  
  years <- c(1850:2008)
  start <- which(years==start.year)
  end <- which(years==end.year)
  
  x <- which(lon.values == lon)
  y <- which(lat.values == lat)
  local.pressure <- slpArray[x,y,start:end]
  
  local.pressure

}

##*******************************************************************************************************

corr.pressure <- function(tot.p, slpArray, start.year, end.year){
  
  nlon <- dim(slpArray)[1]
  nlat <- dim(slpArray)[2]
  corr.array <- array(NA, dim=c(nlon,nlat,4))
  
  years <- c(1850:2008)
  start <- which(years==start.year)
  end <- which(years==end.year)
  
  for(i in 1:nlon){
    for(j in 1:nlat){
       corr <- cor.test(tot.p, slpArray[i,j,start:end], alternative="two.sided",
        method="pearson", na.action=na.exclude)
       if(corr$p.value<=0.05){ # We're only interested in >= 95% confidence
         corr.array[i,j,] <-
          c(corr$estimate,corr$conf.int[1],corr$conf.int[2],corr$p.value)
       }
    }
  }

  corr.array
}

##*******************************************************************************************************

make.image <- function(lon,lat,data.array){
  image.plot(lon,lat,data.array, zlim=c(-0.5,0.5))
  map("worldHires", xlim=c(-100,15), ylim=c(-5,80), interior=F, fill=F, col="grey50", resolution=0, add=T)
  map.axes()
  points(expand.grid(c(-20, -15,-10, -5, 0 , 5, 10), c(40,45,50,55,60)), pch=19, col="grey20")
  points( y=50.1, x=-5.55, pch='*', col="grey20", cex=3)
  points( y=48.38, x=-4.5, pch='*', col="grey20", cex=3)
  points(y=c(45, 45, 45, 45, 40, 40, 40, 50, 50), x=c(-10, -15, -20,
  -5, -15, -10, -20, -20, -15), pch=19, col="yellow")
}

##*******************************************************************************************************

make.image2 <- function(lon,lat,data.array){
  image(lon,lat,data.array, zlim=c(-0.5,0.5), col=tim.colors())
  map("worldHires", xlim=c(-100,15), ylim=c(-5,80), interior=F, fill=F, col="grey50", resolution=0, add=T)
  map.axes()
  points(expand.grid(c(-20, -15,-10, -5, 0 , 5, 10), c(40,45,50,55,60)), pch=19, col="grey20")
  points( y=50.1, x=-5.55, pch='*', col="grey20", cex=3)
  points( y=48.38, x=-4.5, pch='*', col="grey20", cex=3)
  points(y=c(45, 45, 45, 45, 40, 40, 40, 50, 50), x=c(-10, -15, -20,
  -5, -15, -10, -20, -20, -15), pch=19, col="yellow")
}

##*******************************************************************************************************

multi.plot.4 <- function(lon,lat, array1, array2, array3, array4, title.str){
##
## Multi-plot of lagged correlations
##

set.panel()
     
## Here is quick but quirky way to add a common legend to several plots. 
## The idea is leave some room in the margin and then over plot in this margin
     
par(oma=c( 0,0,0,4)) # margin of 4 spaces width at right hand side
set.panel( 2,2) # 2X2 matrix of plots
     
## now draw all your plots using usual image command
make.image2(lon,lat,array1)
title(main=title.str[1])
make.image2(lon,lat,array2)
title(main=title.str[2])
make.image2(lon,lat,array3)
title(main=title.str[3])
make.image2(lon,lat,array4)
title(main=title.str[4])

par(oma=c( 0,0,0,1))# reset margin to be much smaller.
image.plot( legend.only=TRUE, zlim=c(-0.5,0.5)) 
     
## image.plot tricked into  plotting in margin of old setting 
     
set.panel() # reset plotting device

}

##*******************************************************************************************************
sort.corr.pa <- function(corr.array, slp.array, start.year, end.year){

  year <- c(1850:2008)

  ps <- which(year==start.year)
  pe <- which(year==end.year)
  
  sorted.corr <- sort(abs(corr.array),decreasing=T, index=T)
  
  tmp <- slp.array
  i <- dim(tmp)[1]
  j <- dim(tmp)[2]
  k <- dim(tmp)[3]
  
  dim(tmp) <- c((i*j),k)
  
  sorted.pa <- tmp[sorted.corr$ix[1:20],ps:pe] # Only return top 20 correlations

  t(sorted.pa)
}

##*******************************************************************************************************

#########################
## Non-functional part ##
#########################

library(fields)
library(maps)
library(mapdata)

nlon<-72
nlat<-37
nyr<-160
lon <- seq(from=-100,to=15,by=5)
lat <- seq(from=-5,to=80, by=5)
newlyn.start.year.partial <- 1953
newlyn.end.year.partial <- 2008
newlyn.start.year.pred <- 1916
newlyn.end.year.pred <- 1943

##*******************************************************************************************************
## N Atlantic pressure

slpMonthlyArray <- monthly.slp("~/diskx/HadSLP2r/hadSLP2_kij_1850-2009.bin", nlon, nlat, nyr)
nAtlLonLatIndex <- natl.lon.lat()
slpNAtlArray <- natl.slp(slpMonthlyArray, nAtlLonLatIndex$lon, nAtlLonLatIndex$lat)
slpNAtlAnnArray <- annual.slp(slpNAtlArray)

##*******************************************************************************************************
## Newlyn

## Partial
newlyn.sl.partial <- load.tg.data(newlyn.start.year.partial,newlyn.end.year.partial,"newlyn")
newlyn.pa.partial <-
  extract.local.pressure(-5,50,slpNAtlAnnArray,range(lon),range(lat),
                         newlyn.start.year.partial, newlyn.end.year.partial)

newlyn.tot.p.partial <- total.pressure(newlyn.sl.partial$Height,newlyn.pa.partial)

newlyn.corr.partial <- corr.pressure(newlyn.tot.p.partial, slpNAtlAnnArray,
                             newlyn.start.year.partial, newlyn.end.year.partial)

newlyn.lag1.corr.partial <- corr.pressure(newlyn.tot.p.partial, slpNAtlAnnArray,
                             (newlyn.start.year.partial-1), (newlyn.end.year.partial-1))

newlyn.lag2.corr.partial <- corr.pressure(newlyn.tot.p.partial, slpNAtlAnnArray,
                             (newlyn.start.year.partial-2), (newlyn.end.year.partial-2))

newlyn.lag10.corr.partial <- corr.pressure(newlyn.tot.p.partial, slpNAtlAnnArray,
                             (newlyn.start.year.partial-10),
                             (newlyn.end.year.partial-10))

## Prediction
newlyn.sl.pred <- load.tg.data(newlyn.start.year.pred,newlyn.end.year.pred,"newlyn")
newlyn.pa.pred <-
  extract.local.pressure(-5,50,slpNAtlAnnArray,range(lon),range(lat),
                         newlyn.start.year.pred, newlyn.end.year.pred)

newlyn.tot.p.pred <- total.pressure(newlyn.sl.pred$Height,newlyn.pa.pred)

newlyn.corr.pred <- corr.pressure(newlyn.tot.p.pred, slpNAtlAnnArray,
                             newlyn.start.year.pred, newlyn.end.year.pred)

##*******************************************************************************************************
## Sort correlated data

lmRobControl <- lmRob.control(mxr=100,mxf=100,trace=F)

## Partial
newlyn.sorted.pa.partial <- sort.corr.pa(newlyn.corr.partial[,,1], slpNAtlAnnArray, newlyn.start.year.partial, newlyn.end.year.partial)
data.newlyn.partial <- data.frame(msl=newlyn.tot.p.partial,
                                  t=seq(from=newlyn.start.year.partial,to=newlyn.end.year.partial), newlyn.sorted.pa.partial)
tg.lmRob.newlyn.partial <- lmRob(msl ~ ., x=T, y=T, data=data.newlyn.partial, control=lmRobControl, na.action=na.exclude)
## Variance reduction
tg.lmRob.newlyn.partial$r.sq

## Prediction

newlyn.sorted.pa.pred <- sort.corr.pa(newlyn.corr.pred[,,1],
                                      slpNAtlAnnArray,
                                      newlyn.start.year.pred,
                                      newlyn.end.year.pred)
## Singluar matrix encountered with more than 13 pressures here
data.newlyn.pred <- data.frame(msl=newlyn.tot.p.pred,
                                  t=seq(from=newlyn.start.year.pred,to=newlyn.end.year.pred), newlyn.sorted.pa.pred[,1:13])
tg.lmRob.newlyn.pred <- lmRob(msl ~ ., x=T, y=T, data=data.newlyn.pred, control=lmRobControl, na.action=na.exclude)
## Variance reduction
tg.lmRob.newlyn.pred$r.sq

##*******************************************************************************************************
## Brest22
brest22.start.year.partial <- 1850
brest22.end.year.partial <- 2008

brest22.sl <- load.tg.data(brest22.start.year.partial,brest22.end.year.partial,"brest22")
brest22.pa <-
  extract.local.pressure(-5,50,slpNAtlAnnArray,range(lon),range(lat),
                         brest22.start.year.partial, brest22.end.year.partial)

brest22.tot.p <- total.pressure(brest22.sl$Height,brest22.pa)

brest22.corr <- corr.pressure(brest22.tot.p, slpNAtlAnnArray,
                             brest22.start.year.partial, brest22.end.year.partial)


##*******************************************************************************************************
## Plots
##
x11()
make.image(lon,lat,newlyn.corr.partial[,,1])

x11()
multi.plot.4(lon,lat,newlyn.corr.partial[,,1],newlyn.lag1.corr.partial[,,1],newlyn.lag2.corr.partial[,,1],newlyn.lag10.corr.partial[,,1],
             c("Zero lag", "Lag 1", "Lag 2", "Lag 10"))

x11()
make.image(lon,lat,brest22.corr[,,1])

x11()
plot(newlyn.start.year.partial:newlyn.end.year.partial,newlyn.tot.p.partial,col='blue',type='l',
     ylim=c(-1000,3000), xlim=c(newlyn.start.year.pred, newlyn.end.year.partial), ann=F)

lines(tg.lmRob.newlyn.partial$x[,2], tg.lmRob.newlyn.partial$fitted,
      col='cyan')
lines(newlyn.start.year.pred:newlyn.end.year.pred,newlyn.tot.p.pred,col='blue')
lines(tg.lmRob.newlyn.pred$x[,2], tg.lmRob.newlyn.pred$fitted, col='orange')
