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

# Use jet colors for images
#source('~/bin/RScripts/jet.colors.R')
lon<-seq(from=1/3, length=1080, by=1/3)
junk<-read.table('~/diskx/altimetry/lats.lis')
lat<-junk$V1
# jd (Julian Day) in this context appears to be days since 1950 
daytable<-read.table('~/diskx/altimetry/daytable.lis', 
  col.names=c("seq","jd","year","month","day","yearday"))

l<-1080
#l<-360
m<-915
#m<-180
lm<-l*m
n<-639
#n<-513

#lonLatArray<-array(NA,dim=c(lm,3))
#lonLatArray[,1]<-rep(lon,m)
#lonLatArray[,2]<-as.vector(t(array(rep(lat,l),dim=c(m,l))))
trendArray<-array(NA,dim=c(l,m))

# Read the filtered time series that CWH has produced and see which grid
# points have only null values (=9999.0). These can then be flagged as
# land/ice points and removed from the time series.
# The format of hstack_resid (which has annual and semi annual signals removed
# - see the email Altimetry.txt) has time as the first axis, longitude as the
# second and latitude as the third.
#inConn <- file("/diskx/users/simonh/altimetry/hstack_sm5_resid.dat", "rb")
inConn <- file("/diskx/users/cwh/msla2/hstack_resid.dat", "rb")
#inConn <- file("/diskx/users/cwh/msla2/hstack_sm5_1deg_resid.dat", "rb")
for (i in 1:lm) {
  data<-readBin(inConn, what="numeric", n = n, size = 4)
#  data<-readBin(inConn, what="numeric", n = n, size = 4, endian="big")
  if (all(data==9999)) {
#    lonLatArray[i,3]<-1
  } else {
#    message(as.character(i))
#    lonLatArray[i,3]<-0
# Convert cm to mm when calculating the trends
    junk <- lm(data*10 ~ daytable$jd[1:n])
# Convert from mm/day to mm/yr
    trendArray[i] <- junk$coef[2]*365
#    message(paste('Done:', as.character(i),'; Trend:',
#      as.character(trendArray[i])))
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
zonalLat <- seq(from=-50, to=50, by=20)
lenZon <- length(zonalLat)
zonalMeanLat <- vector(mode="numeric", length=lenZon)
for (i in 1:lenZon){
#  junk <- intersect(which(lat>=zonalLat[i]-5), which(lat<zonalLat[i]+5))
  junk <- intersect(which(lat>=zonalLat[i]-10), which(lat<zonalLat[i]+10))
  zonalMeanLat[i] <- mean(zonalMean[junk], na.rm=TRUE)
}

load('../wod04/zonalMeansTS9303.RData')
zonalMeanLatMass <- zonalMeanLat-zonalMeanLatTS9303

# Compare zonal means from Pacific and Atlantic
# Pacific is roughly 115 to 260 and Atlantic is 285 to 20
pacLon <- intersect(which(lon>=115),which(lon<=260))
pacific <- trendArray[pacLon,]

atlLon <- union(which(lon>=285),which(lon<=20))
atlantic <- trendArray[atlLon,]

# Calculate zonal means
zonalAtlanticMean <- vector(mode="numeric", length=m)
zonalPacificMean <- vector(mode="numeric", length=m)
for ( i in 1:m ){
  junk <- trendArray[pacLon,i]
# Only use values with a reasonable rate (+/-50mm/yr)
  j<-intersect(which(junk>=-50),which(junk<=50))
  zonalPacificMean[i] <- mean(junk[j],na.rm=TRUE)
  junk <- trendArray[atlLon,i]
  j<-intersect(which(junk>=-50),which(junk<=50))
  zonalAtlanticMean[i] <- mean(junk[j],na.rm=TRUE)
}

zonalAtlMeanLat <- vector(mode="numeric", length=lenZon)
zonalPacMeanLat <- vector(mode="numeric", length=lenZon)
# Average over latitude bands
for (i in 1:lenZon){
  junk <- intersect(which(lat>=zonalLat[i]-10), which(lat<zonalLat[i]+10))
  zonalAtlMeanLat[i] <- mean(zonalAtlanticMean[junk], na.rm=TRUE)
  zonalPacMeanLat[i] <- mean(zonalPacificMean[junk], na.rm=TRUE)
}

zonalMeanLatMass9303 <- zonalMeanLatMass
zonalLat9303 <- zonalLat
save(zonalMeanLatMass9303, zonalLat9303, file='zonalLatMass9303.RData')

