#######################
# ENACT dimensions:   #
#                     #
# Xmin:   1/3         #
# Xres:   1/3         #
# Xn:    1080         #
#                     #
# Ymin:  -82.0        #
# Yres:  varies w lat #
# Yn:     915         #
#                     #
# Timeslices: 639     #
#                     #
#######################

# Newlyn 50 06 N  05 33 W => 50.1 -5.55
# Newlyn is at x=86, y=90
#> (-5.55-xmin)/xres
#[1] 85.69998
#> (50.1-ymin)/yres
#[1] 89.90001

library(fields)
# Use jet colors for images
#source('~/bin/RScripts/jet.colors.R')
lon<-seq(from=1/3, length=1080, by=1/3)
junk<-read.table('~/diskx/altimetry/lats.lis')
lat<-junk$V1
#
lat180<-c(-89.5:89.5)

# jd (Julian Day) in this context appears to be days since 1950 
daytable<-read.table('~/diskx/altimetry/daytable.lis', 
  col.names=c("seq","jd","year","month","day","yearday"))

#l<-1080
l<-360
#m<-915
m<-180
lm<-l*m
#n<-639
n<-513

trendArray<-array(NA,dim=c(l,m))

# Read the filtered time series that CWH has produced and see which grid
# points have only null values (=9999.0). These can then be flagged as
# land/ice points and removed from the time series.
# The format of hstack_resid (which has annual and semi annual signals removed
# - see the email Altimetry.txt) has time as the first axis, longitude as the
# second and latitude as the third.
inConn <- file("/diskx/users/cwh/msla2/hstack_sm5_1deg_resid.dat", "rb")
for (i in 1:lm) {
#  data<-readBin(inConn, what="numeric", n = n, size = 4)
  data<-readBin(inConn, what="numeric", n = n, size = 4, endian="big")
  if (all(data==9999)) {
  } else {
# Convert cm to mm when calculating the trends
    junk <- lm(data*10 ~ daytable$jd[1:n])
# Convert from mm/day to mm/yr
    trendArray[i] <- junk$coef[2]*365
    if ((i%%100)==0){
      message(paste('Done:', as.character(i)))
    }
  }
}
close(inConn)

# Calculate zonal means
zonalMean <- vector(mode="numeric", length=m)
for ( i in 1:m ){
# Only use values with a reasonable rate (+/-50mm/yr)
  j<-intersect(which(trendArray[,i]>=-50),which(trendArray[,i]<=50))
  zonalMean[i] <- mean(trendArray[j,i],na.rm=TRUE)
}
# Average over latitude bands
#zonalLat <- seq(from=-60, to=60, by=10)
zonalLat <- seq(from=-60, to=60, by=20)
lenZon <- length(zonalLat)
zonalMeanLat <- vector(mode="numeric", length=lenZon)
for (i in 1:lenZon){
#  junk <- intersect(which(lat>=zonalLat[i]-5), which(lat<zonalLat[i]+5))
  junk <- intersect(which(lat180>=zonalLat[i]-10), which(lat180<zonalLat[i]+10))
  zonalMeanLat[i] <- mean(zonalMean[junk], na.rm=TRUE)
}

plot(zonalLat,zonalMeanLat,lwd=2, ylim=c(-2,5))
x11()
image.plot(lon,lat,trendArray,zlim=c(-20,30))

