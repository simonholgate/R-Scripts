# brestNewlynWindStress.R
# 
# Reads the fields of monthly mean winds (calculated from the 6 hourly winds 
# used to force POLCOMS) and extract the time series at Brest and Newlyn of
# the N and E components. 
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

newlynMetMonthlyMean <- array(NA,dim=c(12,5,numYears))
brestMetMonthlyMean <- array(NA,dim=c(12,5,numYears))

nmons<-12

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
  
  dim(SLPMonthlyMean) <- c(m,l,nmons)
  dim(windEMonthlyMean) <- c(m,l,nmons)
  dim(windNMonthlyMean) <- c(m,l,nmons)
  
  wsp <- sqrt (windEMonthlyMean^2 + windNMonthlyMean^2)
  cdw <- as.double(0.63) + as.double(0.066)*wsp
  row <- 1027
  roa <- as.double(1.25)
  fsmet <- roa * cdw * 1.0e-3 * windEMonthlyMean * wsp / row
  gsmet <- roa * cdw * 1.0e-3 * windNMonthlyMean * wsp / row
  
  newlynMetMonthlyMean[,,jj] <- 
    cbind(SLPMonthlyMean[21,11,], windEMonthlyMean[21,11,], windNMonthlyMean[21,11,], 
      fsmet[21,11,], gsmet[21,11,])
  brestMetMonthlyMean[,,jj] <- 
    cbind(SLPMonthlyMean[21,9,], windEMonthlyMean[21,9,], windNMonthlyMean[21,9,],
      fsmet[21,9,], gsmet[21,9,])
  
 }
 save(file="brestNewlynMetMonthlyMean.RData", brestMetMonthlyMean, newlynMetMonthlyMean)