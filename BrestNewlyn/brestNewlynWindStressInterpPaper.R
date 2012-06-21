###############################################################################
# brestNewlynWindStressInterp.R
# 
# Reads the fields of monthly mean winds (calculated from the 6 hourly winds 
# used to force POLCOMS) and extract the time series at Brest and Newlyn of
# the N and E components. 
#
# This version additionally interpolates the 1 degree grid onto the POLCOMS 
# grid
#
# Along with extracting the pressures at Brest and Newlyn
# a regular grid of 16 pressure stations are extracted
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
# Newlyn 50 06 N  05 33 W => 50.1 -5.55
# Newlyn is at x=21, y=11 
# Brest 48 23 N  04 30 W
# Brest is at x=21, y=9

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

newlynMetMonthlyMeanInterp <- array(NA,dim=c(12,5,numYears))
brestMetMonthlyMeanInterp <- array(NA,dim=c(12,5,numYears))
extraMetMonthlyMeanInterp <- array(NA,dim=c(12,5,numYears,16))

nmons<-12
make.surface.grid( list( modelLon,modelLat)) -> loc

# Add extra stations from around the model domain here similar to Thompson (1980)
# The 16 stations are on a regular 5 deg. grid covering 40N-55N, 10W-5E
# as I can do that with this dataset.
# Sites 1-16 are then:
# 40N 10W (40-ymin)/yres (-10-xmin)/xres -> 1 59
# 40N 5W (40-ymin)/yres (-5-xmin)/xres -> 1 89
# 40N 0W (40-ymin)/yres (0-xmin)/xres -> 1 119
# 40N 5E (40-ymin)/yres (5-xmin)/xres -> 1 149
# 45N 10W (45-ymin)/yres (-10-xmin)/xres -> 44 59
# 45N 5W (45-ymin)/yres (-5-xmin)/xres -> 44 89
# 45N 0W (45-ymin)/yres (0-xmin)/xres -> 44 119
# 45N 5E (45-ymin)/yres (5-xmin)/xres -> 44 149
# 50N 10W (50-ymin)/yres (-10-xmin)/xres -> 89 59
# 50N 5W (50-ymin)/yres (-5-xmin)/xres -> 89 89
# 50N 0W (50-ymin)/yres (0-xmin)/xres -> 89 119
# 50N 5E (50-ymin)/yres (5-xmin)/xres -> 89 149
# 55N 10W (55-ymin)/yres (-10-xmin)/xres -> 134 59
# 55N 5W (55-ymin)/yres (-5-xmin)/xres -> 134 89
# 55N 0W (55-ymin)/yres (0-xmin)/xres -> 134 119
# 55N 5E (55-ymin)/yres (5-xmin)/xres -> 134 149
extraMetStns <- c(1,59,1,89,1,119,1,149,44,59,44,89,44,119,44,149,
                  89,59,89,89,89,119,89,149,134,59,134,89,134,119,134,149)
dim(extraMetStns) <- c(2,16)

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
  
  newlynMetMonthlyMeanInterp[,,jj] <- 
    cbind(SLPMonthlyMeanInterp[86,90,], windEMonthlyMeanInterp[86,90,], windNMonthlyMeanInterp[86,90,], 
      fsmet[86,90,], gsmet[86,90,])
  brestMetMonthlyMeanInterp[,,jj] <- 
    cbind(SLPMonthlyMeanInterp[92,74,], windEMonthlyMeanInterp[92,74,], windNMonthlyMeanInterp[92,74,],
      fsmet[92,74,], gsmet[92,74,])


  for (k in 1:16) {
    extraMetMonthlyMeanInterp[,,jj,k] <- 
      cbind(SLPMonthlyMeanInterp[extraMetStns[1,k],extraMetStns[2,k],], 
        windEMonthlyMeanInterp[extraMetStns[1,k],extraMetStns[2,k],], 
        windNMonthlyMeanInterp[extraMetStns[1,k],extraMetStns[2,k],], 
        fsmet[extraMetStns[1,k],extraMetStns[2,k],], 
        gsmet[extraMetStns[1,k],extraMetStns[2,k],])
  }
 }
 save(file="brestNewlynMetMonthlyMeanInterpPaper.RData", brestMetMonthlyMeanInterp, 
   newlynMetMonthlyMeanInterp, extraMetMonthlyMeanInterp)

