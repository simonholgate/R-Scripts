###############################################################################
# spanishWindStressInterp.R
# 
# Reads the fields of monthly mean winds (calculated from the 6 hourly winds 
# used to force POLCOMS) and extract the time series at Spanish stations of
# the N and E components. 
#
# This version additionally interpolates the 1 degree grid onto the POLCOMS 
# grid
#
# Along with extracting the pressures at Atlantic Spanish stations, following 
# Thompson (1980), I extract 9 pressure stations
#
# These winds are then converted to wind stresses
#
# Author: simonh
###############################################################################


#######################
# Met dimensions:     #
#                     #
# Xmin:   -25E        #
# Xres:   1           #
# Xn:     41          #
#                     #
# Ymin:   40N         #
# Yres:   1           #
# Yn:     26          #
#######################


#######################
# _*S2 dimensions:    #
#                     #
# Xmin:   -19.83333   #
# Xres:   1/6         #
# Xn:     198         #
#                     #
# Ymin:   40.11111    #
# Yres:   1/9         #
# Yn:     224         #
#######################
xmin <- -19.83333
ymin <- 40.11111
xres <- 1/6
yres <- 1/9

# Get station points
load("~/diskx/polcoms/spain/spanishStnsLatLon.RData")

modelLat<-seq(from=40.11111,by=1/9,length=224)
modelLon<-seq(from=-19.83333,by=1/6,length=198)
ll <- 224
mm <- 198

library(fields)
lat<-seq(from=40,by=1,length=26)
lon<-seq(from=-25,by=1,length=41)

l<-26
m<-41
lm<-l*m
#n<-66

leapYrDaysArray<-array(data=
  c(31,29,31,30,31,30,31,31,30,31,30,31),
  dim=c(12,1))
daysArray<-array(data=
  c(31,28,31,30,31,30,31,31,30,31,30,31),
  dim=c(12,1))
yearsArray<-array(data=
  c(1960,1961,1962,1963,1964,1965,1966,1967,1968,1969,
    1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,
    1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,
    1990,1991,1992,1993,1994,1995,1996,1997,1998,1999), dim=c(40,1))
#   c(1960), dim=c(1,1))

monthArray<-array(data=
  c("01","02","03","04","05","06","07","08","09","10","11","12"),
  dim=c(12,1))

numYears <- length(yearsArray)

# Just worry about the 12 Atlantic stations
spanishMetMonthlyMeanInterp <- array(NA,dim=c(12,5,numYears,12))
# And the extra stations in the model domain...
extraMetMonthlyMeanInterp <- array(NA,dim=c(12,5,numYears,17))

nmons<-12
make.surface.grid( list( modelLon,modelLat)) -> loc

# Add extra stations from around the model domain here similar to Thompson (1980)
# St Mawgan - Latitude: 50.454243N  Longitude: 4.99915W
#> (-4.99915-xmin)/xres -> 89
#> (50.454243-ymin)/yres -> 93
# Thorney Island - latitude 50.8166667 longitude -0.9166667
#> (-0.9166667-xmin)/xres -> 114
#> (50.8166667-ymin)/yres -> 96
# Shoeburyness - Latitude: 51.55; Longitude: 0.833
#> (0.833-xmin)/xres -> 124
#> (51.55-ymin)/yres -> 103
# Gorleston - Latitude, 52.5833, Longitude, 1.7167
#> (1.7167-xmin)/xres -> 129
#> (52.5833-ymin)/yres -> 112
# Kilnsea - Latitude, 53.6167, Longitude, 0.1333
#> (0.1333-xmin)/xres -> 120
#> (53.6167-ymin)/yres -> 122
# Eskdalemuir - Latitude: 55.317; Longitude: -3.2
#> (-3.2-xmin)/xres -> 100
#> (55.317-ymin)/yres -> 137
# Kinloss - Latitude: 57.65; Longitude: -3.567
#> (-3.567-xmin)/xres -> 98
#> (57.65-ymin)/yres -> 158
# Ronaldsway - Latitude: 54.07620. Longitude: -4.62333
#> (-4.62333-xmin)/xres -> 91
#> (54.07620-ymin)/yres -> 126
# Mumbles - Latitude: 51.567; Longitude: -3.983
#> (-3.983-xmin)/xres -> 95
#> (51.567-ymin)/yres -> 103
# The next eight stations are slightly aribtrary as I don't have Thompson's
# actually lats and lons for these. I'm also adding a site and making them
# more regular as I can with this dataset.
# Sites 10-18 are then:
# 45N 20W (45-ymin)/yres (-20-xmin)/xres -> 44 1,
# 45N 10W (45-ymin)/yres (-10-xmin)/xres -> 44 59,
# 45N 0W (45-ymin)/yres (0-xmin)/xres -> 44 119,
# 55N 10E (55-ymin)/yres (10-xmin)/xres -> 134 179,
# 65N 10E (65-ymin)/yres (10-xmin)/xres -> 224 179,
# 65N 10W (65-ymin)/yres (-10-xmin)/xres -> 224 59,
# 65N 20W (65-ymin)/yres (-20-xmin)/xres -> 224 1,
# 55N 20W (55-ymin)/yres (-20-xmin)/xres -> 134 1
extraMetStns <- c(89,93,114,96,124,103,129,112,120,122,100,137,98,158,
  91,126,95,103,1,44,59,44,119,44,179,134,179,224,59,224,1,224,1,134)
dim(extraMetStns) <- c(2,17)

for (jj in 1:numYears) {

  inSLPMonthlyMean <- file(paste("~/diskx/polcoms/met_data/monthlyMeans/slp.mm.",
    yearsArray[jj],".dat",sep=""), "rb")

  inWindEMonthlyMean <- file(paste("~/diskx/polcoms/met_data/monthlyMeans/wind.E.mm.",
    yearsArray[jj],".dat",sep=""), "rb")

  inWindNMonthlyMean <- file(paste("~/diskx/polcoms/met_data/monthlyMeans/wind.N.mm.",
    yearsArray[jj],".dat",sep=""), "rb")

  SLPMonthlyMean <- readBin(inSLPMonthlyMean, "numeric", n = nmons*lm)
  windEMonthlyMean <- readBin(inWindEMonthlyMean, "numeric", n = nmons*lm)
  windNMonthlyMean <- readBin(inWindNMonthlyMean, "numeric", n = nmons*lm)
  
  close(inSLPMonthlyMean)
  close(inWindEMonthlyMean)
  close(inWindNMonthlyMean)
  
  dim(SLPMonthlyMean) <- c(m,l,nmons)
  dim(windEMonthlyMean) <- c(m,l,nmons)
  dim(windNMonthlyMean) <- c(m,l,nmons)
  
  SLPMonthlyMeanInterp <- array(NA, dim=c(mm,ll,nmons))
  windEMonthlyMeanInterp <- array(NA, dim=c(mm,ll,nmons))
  windNMonthlyMeanInterp <- array(NA, dim=c(mm,ll,nmons))
  
#  Interpolation here
  
  for (kk in 1:nmons) {
    obj<- list( x= lon, y=lat, z= windEMonthlyMean[,,kk])
    interp.surface( obj, loc) -> windEMonthlyMeanInterp[,,kk]
    
    obj<- list( x= lon, y=lat, z= windNMonthlyMean[,,kk])
    interp.surface( obj, loc) -> windNMonthlyMeanInterp[,,kk]
    
    obj<- list( x= lon, y=lat, z= SLPMonthlyMean[,,kk])
    interp.surface( obj, loc) -> SLPMonthlyMeanInterp[,,kk]
  }
  
  wsp <- sqrt (windEMonthlyMeanInterp^2 + windNMonthlyMeanInterp^2)
  cdw <- as.double(0.63) + as.double(0.066)*wsp
  row <- 1027
  roa <- as.double(1.25)
  fsmet <- roa * cdw * 1.0e-3 * windEMonthlyMeanInterp * wsp / row
  gsmet <- roa * cdw * 1.0e-3 * windNMonthlyMeanInterp * wsp / row
  
#  newlynMetMonthlyMeanInterp[,,jj] <- 
#    cbind(SLPMonthlyMeanInterp[86,90,], windEMonthlyMeanInterp[86,90,], windNMonthlyMeanInterp[86,90,], 
#      fsmet[86,90,], gsmet[86,90,])


  for (k in 1:12) {
    spanishMetMonthlyMeanInterp[,,jj,k] <- 
      cbind(SLPMonthlyMeanInterp[spanishStnsXY[k,1],spanishStnsXY[k,2],], 
        windEMonthlyMeanInterp[spanishStnsXY[k,1],spanishStnsXY[k,2],], 
        windNMonthlyMeanInterp[spanishStnsXY[k,1],spanishStnsXY[k,2],], 
        fsmet[spanishStnsXY[k,1],spanishStnsXY[k,2],], 
        gsmet[spanishStnsXY[k,1],spanishStnsXY[k,2],])
  }

  for (k in 1:17) {
    extraMetMonthlyMeanInterp[,,jj,k] <- 
      cbind(SLPMonthlyMeanInterp[extraMetStns[1,k],extraMetStns[2,k],], 
        windEMonthlyMeanInterp[extraMetStns[1,k],extraMetStns[2,k],], 
        windNMonthlyMeanInterp[extraMetStns[1,k],extraMetStns[2,k],], 
        fsmet[extraMetStns[1,k],extraMetStns[2,k],], 
        gsmet[extraMetStns[1,k],extraMetStns[2,k],])
  }
 }
 save(file="spanishMetMonthlyMeanInterp.RData", spanishMetMonthlyMeanInterp, extraMetMonthlyMeanInterp)

