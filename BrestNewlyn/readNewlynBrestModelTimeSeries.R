###############################################################################
# Variation of readzet.r to read the 2D daily mean sea level calculated from
# POLCOMMS model output. 
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

# Also choose 6 boundary points to reflect the deep ocean boundary condition
# Choose points at:
# x=1, y=50, x=1, y=100, x=1, y=150
# x=50, y=50, x=50, y=100, x=50, y=150,

load("~/diskx/polcoms/iseajseanpsea.Rdata")

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

newlynMonthlyMean<-array(NA,dim=c(12*n,1))
brestMonthlyMean<-array(NA,dim=c(12*n,1))
deepSeaMonthlyMean<-array(NA,dim=c(12*n,6))

deepSeaPoints <- c(1,1,1,50,50,50,50,100,150,50,100,150)
dim(deepSeaPoints) <-c(6,2) 

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
    
    newlynMonthlyMean[monthCount]<-zettMonthlyMean[86,90]    
    brestMonthlyMean[monthCount]<-zettMonthlyMean[92,74]
    for (k in 1:6){
      deepSeaMonthlyMean[monthCount,k] <-
        zettMonthlyMean[deepSeaPoints[k,1], deepSeaPoints[k,2]]
    }
    
    monthCount<-monthCount+1
  }
    
  close(inMonthlyMean)
}
#
monthsArray <- seq.Date(from=as.Date("1960/1/15"), to=as.Date("1999/12/15"), by="1 month")
#
zettMonthlyMean[which(zettMonthlyMean==0)]<-NA
image.plot(lon,lat,zettMonthlyMean[c(1:198),c(1:224)],zlim=c(-1,1),
  col=jet.colors(100))
#
x11()
par(family="HersheySans")
plot(monthsArray, brestMonthlyMean, type='l', col='blue', ylim=c(-0.5,0.0))
lines(monthsArray, newlynMonthlyMean, col='red')

#
dmBrestMonthlyMean <- brestMonthlyMean-mean(brestMonthlyMean, na.rm=T)
dmNewlynMonthlyMean <- newlynMonthlyMean-mean(newlynMonthlyMean, na.rm=T)
dmDeepSeaMonthlyMean <- deepSeaMonthlyMean
for (k in 1:6){
 dmDeepSeaMonthlyMean[,k] <- 
   deepSeaMonthlyMean[,k]-mean(deepSeaMonthlyMean[,k])
}

x11()
par(family="HersheySans")
plot(monthsArray, (dmBrestMonthlyMean-dmNewlynMonthlyMean), type='l', 
  col='blue')
filtModelDiff<-filter((dmBrestMonthlyMean-dmNewlynMonthlyMean), 
  trian13ptFilter,method="c",sides=2)
lines(monthsArray, filtModelDiff, col='red')