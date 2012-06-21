###############################################################################
# Variation of readzet.r to read the 2D daily mean sea level calculated from
# POLCOMS model output. 
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

# Delfzijl 53.333 6.933
# Delfzijl is at x=161, y=119
#> (6.933-xmin)/xres
#[1] 160.5980
#> (53.333-ymin)/yres
#[1] 118.997

load("~/diskx/polcoms/iseajseanpsea.Rdata")

library(fields)
# Use jet colors for images
source('~/bin/RScripts/jet.colors.R')
lat<-seq(from=40.11111,by=1/9,length=224)
lon<-seq(from=-19.83333,by=1/6,length=198)

l<-198
m<-224
lm<-l*m
n<-45
#
yearsArray<-array(data=
  c(1960,1961,1962,1963,1964,1965,1966,1967,1968,1969,
    1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,
    1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,
    1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,
    2000,2001,2002,2003,2004),
  dim=c(n,1))

delfzijlMonthlyMean<-array(NA,dim=c(12*n,1))

zettMonthlyMean<-array(NA,dim=c(l,m))
monthCount<-1

for (i in 1:n){
  inMonthlyMean <- file(
    paste("~/diskx/polcoms/S12run408/S12run408ZetFiles/zett.mm.",yearsArray[i],".dat",sep=""), "rb")

# Read monthly means in
  for (j in 1:12){
    zett<-readBin(inMonthlyMean,n=npsea, what='numeric', size=4)
    for (ip in 1:npsea) {
      zettMonthlyMean[isea[ip],jsea[ip]]<-zett[ip]
    }
    delfzijlMonthlyMean[monthCount]<-zettMonthlyMean[86,90]    
    monthCount<-monthCount+1
  }
    
  close(inMonthlyMean)
}
#
monthsArray <- seq.Date(from=as.Date("1960/1/15"), to=as.Date("2004/12/15"), by="1 month")
#
zettMonthlyMean[which(zettMonthlyMean==0)]<-NA
image.plot(lon,lat,zettMonthlyMean[c(1:198),c(1:224)],zlim=c(-1,1),
  col=jet.colors(100))
#
x11()
par(family="HersheySans")
plot(monthsArray, delfzijlMonthlyMean, type='l', col='blue', ylim=c(-0.5,0.0))

#
dmDelfzijlMonthlyMean <- delfzijlMonthlyMean-mean(delfzijlMonthlyMean, na.rm=T)
 
