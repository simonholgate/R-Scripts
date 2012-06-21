###############################################################################
# Variation of readNewlynBrestModelTimeSeries.R to read the 2D daily mean sea 
# level calculated from POLCOMMS model output along the TOPEX/Jason altimetry
# tracks provided by Laurent Testut (see output from readModelTracks.R). 
# 
# The arrays were written from R with the makeZetDailyMeans.r script and are
# little endian real*4 numbers. Each month contains (day in month - 2) daily
# means.
#
# Simon Holgate, August 2005
#
###############################################################################

#######################
# _*S2 dimensions:    #
#                     #
# *_Xmin:   -19.83333 #
# Xres:   1/6         #
# Xn:     198         #
#                     #
# Ymin:   40.11111    #
# Yres:   1/9         #
# Yn:     224         #
#######################

# Newlyn 50 06 N  05 33 W => 50.1 -5.55
# Newlyn is at x=86, y=90
#> (-5.55-xmin)/xres
#[1] 85.69998
#> (50.1-ymin)/yres
#[1] 89.90001
# Brest 48 23 N  04 30 W
# Brest is at x=92, y=74
#> (-4.5-min)/xres
#[1] 91.99998
#> (48.38-ymin)/yres
#[1] 74.42001
load("~/diskx/polcoms/iseajseanpsea.Rdata")
load("~/diskx/polcoms/brestNewlyn/tracks/trackLatLon.RData")

library(fields)
# Use jet colors for images
source('~/bin/RScripts/jet.colors.R')
lat<-seq(from=40.11111,by=1/9,length=224)
lon<-seq(from=-19.83333,by=1/6,length=198)

l<-198
m<-224
lm<-l*m
n<-40
#
yearsArray<-array(data=
  c(1960,1961,1962,1963,1964,1965,1966,1967,1968,1969,
    1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,
    1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,
    1990,1991,1992,1993,1994,1995,1996,1997,1998,1999),
  dim=c(n,1))

trackMonthlyMean<-array(NA,dim=c(12*n,length(nearest)))
track61MonthlyMean<-array(NA,dim=c(12*n,length(nearest61)))
track70MonthlyMean<-array(NA,dim=c(12*n,length(nearest70)))
track239MonthlyMean<-array(NA,dim=c(12*n,length(nearest239)))
track248MonthlyMean<-array(NA,dim=c(12*n,length(nearest248)))


zettMonthlyMean<-array(NA,dim=c(l,m))
monthCount<-1

for (i in 1:n){
  inMonthlyMean <- file(
    paste("~/diskx/polcoms/S12run405ZetFiles/zett.mm.",yearsArray[i],".dat",sep=""), "rb")

# Read monthly means in
  for (j in 1:12){
    zett<-readBin(inMonthlyMean,n=npsea, what='numeric', size=4)
    for (ip in 1:npsea) {
      zettMonthlyMean[isea[ip],jsea[ip]]<-zett[ip]
    }
    
    trackMonthlyMean[monthCount,]<-zettMonthlyMean[nearest]
    track61MonthlyMean[monthCount,]<-zettMonthlyMean[nearest61]
    track70MonthlyMean[monthCount,]<-zettMonthlyMean[nearest70]  
    track239MonthlyMean[monthCount,]<-zettMonthlyMean[nearest239]
    track248MonthlyMean[monthCount,]<-zettMonthlyMean[nearest248]
      
    monthCount<-monthCount+1
  }
    
  close(inMonthlyMean)
}
#
monthsArray <- seq.Date(from=as.Date("1960/1/15"), to=as.Date("1999/12/15"), by="1 month")
#
zettMonthlyMean[which(zettMonthlyMean==0)]<-NA
x11()
par(family="HersheySans")
image.plot(lon,lat,zettMonthlyMean[c(1:198),c(1:224)],zlim=c(-1,1),
  col=jet.colors(100))
world(add=T)

#plot(monthsArray, brestMonthlyMean, type='l', col='blue', ylim=c(-0.5,0.0))
x11()
split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
par(family="HersheySans")
screen(1)
plot(lon_vec,lat_vec, type='n', ann=F) 
title(xlab='Lon', ylab='Lat', main='Tracks 70, 61, 239 & 248')
ribbon.plot(altimTrackPoints[,1], altimTrackPoints[,2], trackMonthlyMean[1,], 
  ylim=c(-0.5,0.0), lwd=10)
world(add=T)

screen(2)
image.plot(zlim=c(-0.5,0),legend.only=TRUE, smallplot=c(.1,.2, .3,.7),
  col=tim.colors(), legend.lab='Height [mm]')

close.screen( all=TRUE)

x11()
split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
par(family="HersheySans")
screen(1)
plot(lon_vec,lat_vec, type='n', ann=F) 
title(xlab='Lon', ylab='Lat', main='Tracks 70, 61, 239 & 248')
ribbon.plot(altimTrackPoints61[,1], altimTrackPoints61[,2], track61MonthlyMean[1,], 
  ylim=c(-0.5,0.0), lwd=10)
ribbon.plot(altimTrackPoints70[,1], altimTrackPoints70[,2], track70MonthlyMean[1,], 
  ylim=c(-0.5,0.0), lwd=10)
ribbon.plot(altimTrackPoints239[,1], altimTrackPoints239[,2], track239MonthlyMean[1,], 
  ylim=c(-0.5,0.0), lwd=10)
ribbon.plot(altimTrackPoints248[,1], altimTrackPoints248[,2], track248MonthlyMean[1,], 
  ylim=c(-0.5,0.0), lwd=10)
world(add=T)

screen(2)
image.plot(zlim=c(-0.5,0),legend.only=TRUE, smallplot=c(.1,.2, .3,.7),
  col=tim.colors(), legend.lab='Height [mm]')

close.screen( all=TRUE)

save(file="modelAltimTracks.RData", altimTrackPoints, trackMonthlyMean, 
  altimTrackPoints61, track61MonthlyMean, altimTrackPoints70, track70MonthlyMean,
  altimTrackPoints239, track239MonthlyMean, altimTrackPoints248, track248MonthlyMean)