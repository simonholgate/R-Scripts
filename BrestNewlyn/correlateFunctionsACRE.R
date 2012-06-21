#############################
## Functions for use below ##
#############################

monthly.slp.ncdf <- function(path, nlon, nlat, nyr){
  ## Read the sea level data into an array
  nmon <- 12

#  slpMonthlyArray<-array(NA,dim=c(nlon,nlat,nmon,nyr))

#  con<-file(path,"rb")

#  for (i in 1:nlat){
#    for (j in 1:nlon){
#      for (k in 1:nyr){
#        slp<-readBin(con,what="numeric", size=4, n=nmon, endian='little')
#        slpMonthlyArray[j,i,,k]<-slp
#      }
#    }
#  }
#  close(con)
  
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
#  xlon<-seq(from=-180,by=5,length=nlon)
#  ylat<-seq(from=90,by=-5,length=nlat)
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

  # N Pacific covers 90E to 250E and 5S to 80N
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
        if(station=='keywest'){
          tgs <- which(tg$keywest$Year==tg.start.year) # index of start of Key West TG data
          tge <- which(tg$keywest$Year==tg.end.year) # index of end of Key West TG data
          sea.level <- tg$keywest[tgs:tge,]
        } else {
          if(station=='sanfran'){
            tgs <- which(tg$sanfran$Year==tg.start.year) # index of
                                        # start of San Francsico TG data
            tge <- which(tg$sanfran$Year==tg.end.year) # index of end
                                        # of San Francisco TG data
            sea.level <- tg$sanfran[tgs:tge,]
          } else {
            error(paste("station not known:", station))
          }
        }
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

##*******************************************************************************************************

interp.local.pressure <- function(lon,lat, slpArray, lon.range,
                                   lat.range, start.year, end.year, years){
  lon.values <- seq(from=lon.range[1], to=lon.range[2], by=2)
  lat.values <- seq(from=lat.range[1], to=lat.range[2], by=2)
  
  start <- which(years==start.year)
  end <- which(years==end.year)
  nyrs <- end - start

  # Check whether lon runs from 0 to 360 or -180 to 180. We want it to be -180 to 180.

  if ((lon >= 180) && (lon.range[1]<0)){
    lon <- lon - 360
  } 

  # For each of the nyrs, we need to interpolate the grid to the locations we've chosen
  local.pressure <- vector(length=nyrs, mode="numeric")
  for (i in start:end) {
    obj <- list(x=lon.values, y=lat.values, z=slpArray[,,i])
    local.pressure[i-start+1] <- interp.surface(obj,cbind(lon,lat))
  }

  local.pressure

}

##*******************************************************************************************************

corr.pressure <- function(tot.p, slpArray, start.year, end.year, years){
  
  nlon <- dim(slpArray)[1]
  nlat <- dim(slpArray)[2]
  corr.array <- array(NA, dim=c(nlon,nlat,4))
  
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
## Robustly linearly detrend a vector. Missing values must be included as NA, but are excluded
## in the regression.
detrend <- function(data_vector){
  df.data_vector <- data.frame(p=data_vector,t=seq(from=1,to=length(data_vector)))
  data.lmRob <- lmRob(p ~ t, x=T, y=T, data=df.data_vector, control=lmRobControl,
                              na.action=na.exclude)
## We want the returned vector to be the same length as the original vector so use the $x values
  resid.data_vector <- vector(mode="numeric", length=length(data_vector))
  resid.data_vector[data.lmRob$x[,2]] <- data.lmRob$residuals 

  resid.data_vector
}

##*******************************************************************************************************
## Correlate detrended pressure
corr.pressure.dt <- function(tot.p, slpArray, start.year, end.year, years){
  
  nlon <- dim(slpArray)[1]
  nlat <- dim(slpArray)[2]
  corr.array <- array(NA, dim=c(nlon,nlat,4))
  
  start <- which(years==start.year)
  end <- which(years==end.year)
  
  for(i in 1:nlon){
    for(j in 1:nlat){
       corr <- cor.test(tot.p, detrend(slpArray[i,j,start:end]), alternative="two.sided",
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
  image.plot(lon,lat,data.array, zlim=c(-0.6,0.6))
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
  image(lon,lat,data.array, zlim=c(-0.6,0.6), col=tim.colors())
  map("worldHires", xlim=c(-100,15), ylim=c(-5,80), interior=F, fill=F, col="grey50", resolution=0, add=T)
  map.axes()
  points(expand.grid(c(-20, -15,-10, -5, 0 , 5, 10), c(40,45,50,55,60)), pch=19, col="grey20")
  points( y=50.1, x=-5.55, pch='*', col="grey20", cex=3)
  points( y=48.38, x=-4.5, pch='*', col="grey20", cex=3)
  points(y=c(45, 45, 45, 45, 40, 40, 40, 50, 50), x=c(-10, -15, -20,
  -5, -15, -10, -20, -20, -15), pch=19, col="yellow")
}

##*******************************************************************************************************

make.image3 <- function(lon,lat,data.array){
  image.plot(lon,lat,data.array, zlim=c(-0.6,0.6))
  world(add=T)
  points(c(-20, -15,-10, -5, 0 , 5, 10), c(40,45,50,55,60), pch=19, col="grey20")
  points( y=50.1, x=-5.55, pch='*', col="grey20", cex=3)
  points( y=48.38, x=-4.5, pch='*', col="grey20", cex=3)
  points(y=c(45, 45, 45, 45, 40, 40, 40, 50, 50), x=c(-10, -15, -20,
  -5, -15, -10, -20, -20, -15), pch=19, col="yellow")
}

##*******************************************************************************************************
# As make.image3 but without a fixed zlim
make.image4 <- function(lon,lat,data.array){
  image.plot(lon,lat,data.array)
  world(add=T, shift=T)
  points(expand.grid(c(-20, -15,-10, -5, 0 , 5, 10), c(40,45,50,55,60)), pch=19, col="grey20")
  points( y=50.1, x=-5.55, pch='*', col="grey20", cex=3)
  points( y=48.38, x=-4.5, pch='*', col="grey20", cex=3)
  points( y=37.8, x=360-122.467, pch='*', col="magenta", cex=3)
  points(y=c(45, 45, 45, 45, 40, 40, 40, 50, 50), x=c(-10, -15, -20,
  -5, -15, -10, -20, -20, -15), pch=19, col="yellow")
}

##*******************************************************************************************************

make.image5 <- function(lon,lat,data.array){
  image.plot(lon,lat,data.array, zlim=c(-0.6,0.6))
  world(add=T)
#  points(c(-20, -15,-10, -5, 0 , 5, 10), c(40,45,50,55,60)), pch=19, col="grey20")
  points( y=50.1, x=-5.55, pch='*', col="grey20", cex=3)
  points( y=48.38, x=-4.5, pch='*', col="grey20", cex=3)
  points( y=37.8, x=-122.467, pch='*', col="grey20", cex=3)
  points(-25.7,35.75, pch=19, col="magenta", cex=2)
#  points(y=c(45, 45, 45, 45, 40, 40, 40, 50, 50), x=c(-10, -15, -20,
#  -5, -15, -10, -20, -20, -15), pch=19, col="yellow")
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
image.plot( legend.only=TRUE, zlim=c(-0.6,0.6)) 
     
## image.plot tricked into  plotting in margin of old setting 
     
set.panel() # reset plotting device

}

##*******************************************************************************************************
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

##*******************************************************************************************************
deriv <- function(in.vector, in.dx){
## Function to take the 1st derivative by central difference method

## Input:
## in.vector - a vector of at least length 2 of real valued numbers
## in.dx - the distance between values, either a single value or a vector of length 1 less 
## than the input vector

## Output:
## out.vector - a real valued vector or single value which is the gradient of the 
## input vector

## Example:
## deriv(c(1,2),2) -> 0.5
## deriv(c(1,2,3), 2) -> c(0.5, 0.5)

  len <- length(in.vector)
  lenivm1 <- len - 1

  lendx <- length(in.dx)
  if ((lendx != 1) && (lendx != lenivm1)) {
    stop(paste("length of dx not equal to 1 or len(iv)-1:", lenivm1))
  }

  if(length(in.dx)==1){
    dx <- array(in.dx, dim=c(lenivm1,1))
  } else {
    dx <- in.dx
  }

  dx <- as.vector(dx)

  out.vector <- vector(mode="numeric",length=lenivm1)

  for (i in 1:lenivm1){
    out.vector[i] <- (in.vector[i+1]-in.vector[i])/dx[i]
  }

  out.vector
}
##*******************************************************************************************************

grad <- function(in.array, in.dx){
## Function to take the 1st derivative by central difference method

## Input:
## in.array - a 2D array of at least dimension c(2,2) of real valued numbers
## in.dx - the distance between values, either a single value
## or a vector of length 2, c(dx,dy) 
## or a 2D array with dimensions that are 1 less than the input array in each direction

## Output:
## out.array - a complex valued 2D array or single complex value which is the 2D gradient of the 
## input array

## Example:
## grad(array(c(1,2,1,2), dim=c(2,2)),2) -> 0+0.5i
## grad(array(c(1,1,2,2), dim=c(3,3)),2) -> 0.5+0i
## grad(array(c(1,1,1,2,2,2,3,3,3), 2) -> array(c(0.5+0i,0.5+0i,0.5+0i,0.5+0i), c(2,2))
## grad(array(c(1,1,1,2,2,2,3,3,3), array(c(2,2,2,2), dim=c(2,2))) -> array(c(0.5+0i,0.5+0i,0.5+0i,0.5+0i), c(2,2))

  dim.in.array <- dim(in.array)

  if((is.null(dim.in.array)) || (length(dim.in.array) > 2)){
    stop("2D array required for grad operator")
  } else {
    len.dx <- dim.in.array[1]
    len.dy <- dim.in.array[2]
  } 

  len.in.dx <- length(in.dx)
  if((len.in.dx != 1) && (len.in.dx != 2) && (len.in.dx != ((len.dx -1)*(len.dy-1)))){
    stop("Dimensions of dx not appropriate")
  }

  if(len.in.dx == 1){
    out.dx <- array(data=in.dx, dim=c((len.dx-1),(len.dy-1))
  } else {
    if(len.in.dx == 2){
      dx <- array(data=in.dx[1], dim=c((len.dx-1),1)
      dy <- array(data=in.dx[2], dim=c((len.dy-1),1)
      out.dx <- expand.grid(x=as.vector(dx), y=as.vector(dy))
    } else {
      out.dx <- expand.grid(x=, y=)
    }
  } 
  
  for(i in 1:len.dx){
    x.deriv <- deriv
  }
  for(j in 1:len.dy){
    y.deriv <- deriv
  }



  out.array

}
##*******************************************************************************************************

