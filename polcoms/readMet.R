###############################################################################
# Variation of readzet.r which calculates the 2D daily mean winds & SLP from the
# meterological data. 
# 
# The met arrays are 6 hourly (4 times a day) and are wind east, wind north,
# pressure, air temp and relative humdity.
#
# Simon Holgate, March 2007
#
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
# Newlyn is at x=31, y=11 

library(fields)
# Use jet colors for images
source('~/bin/RScripts/jet.colors.R')
lat<-seq(from=40,by=1,length=26)
lon<-seq(from=-25,by=1,length=41)

l<-26
m<-41
lm<-l*m
n<-66

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

globalSLPMonthlyMean <- vector(mode="numeric", length=(numYears*12))
globalWindEMonthlyMean <- vector(mode="numeric", length=(numYears*12))
globalWindNMonthlyMean <- vector(mode="numeric", length=(numYears*12))

newlynMonthlyMean <- array(NA,dim=c(12,numYears))

for (jj in 1:numYears) {
# Test for leap year
  if ((yearsArray[jj] %% 4)==0) {
    if ((yearsArray[jj] %% 100)==0) {
      if ((yearsArray[jj] %% 400)==0) {
        leapFlag<-1
      } else {
        leapFlag<-0
      }
    } else {
      leapFlag<-1
    }
    leapFlag<-1
  } else {
    leapFlag<-0
  }

  nmons<-12

  SLPMonthlyMean<-array(NA,dim=c(lm,12))
  SLPDailyMean<-array(NA,dim=c(lm,31))

  windEMonthlyMean<-array(NA,dim=c(lm,12))
  windEDailyMean<-array(NA,dim=c(lm,31))

  windNMonthlyMean<-array(NA,dim=c(lm,12))
  windNDailyMean<-array(NA,dim=c(lm,31))

  outSLPMonthlyMean <- file(paste("monthlyMeans/slp.mm.",
    yearsArray[jj],".dat",sep=""), "wb")

  outWindEMonthlyMean <- file(paste("monthlyMeans/wind.E.mm.",
    yearsArray[jj],".dat",sep=""), "wb")

  outWindNMonthlyMean <- file(paste("monthlyMeans/wind.N.mm.",
    yearsArray[jj],".dat",sep=""), "wb")

  for (ii in 1:nmons) {

    fname <- paste("met", as.character(yearsArray[jj]), monthArray[ii], 
      "_WPTRh.txt", sep="")
    inFile <- file(fname, "r")

    outSLPDailyMean <- file(paste("dailyMeans/slp.dm.",monthArray[ii],".",   
      yearsArray[jj],".dat",sep=""), "wb")

    outWindEDailyMean <- file(paste("dailyMeans/wind.E.dm.", 
      monthArray[ii],".",   yearsArray[jj],".dat",sep=""), "wb")

    outWindNDailyMean <- file(paste("dailyMeans/wind.N.dm.", 
      monthArray[ii],".",   yearsArray[jj],".dat",sep=""), "wb")

    if ( leapFlag == 1) {
  # Leap year
      days<-leapYrDaysArray[ii]
    } else {
  # Regular year
      days<-daysArray[ii]
    }

# Read data into a table
    firstRead <- TRUE
    for (ll in 1:days) {
      inDataSLP <- array(NA, dim=c(lm,4))
      inDataWindE <- array(NA, dim=c(lm,4))
      inDataWindN <- array(NA, dim=c(lm,4))
# With 6 hourly data we need to do 4 reads to get 1 daily mean
# We also need to skip over the air temp and humidity data, once we have read
# the 2 wind components and slp for the first time
# The number of lines to skip is 107 per data variable as there are 10 values
# per line of the input file and there are 41*26 = 1066 values per variable.
# 1066/10 = 107 lines.
      for (kk in 1:4) {
        if (firstRead == TRUE){
          inDataWindE[,kk] <- scan(inFile, what=numeric(), nmax=lm)
          inDataWindN[,kk] <- scan(inFile, what=numeric(), nmax=lm)
          inDataSLP[,kk] <- scan(inFile, what=numeric(), nmax=lm)
          firstRead <- FALSE
        } else {
# So after the first read we need to skip 2 variables (temp and RH) of 107
# lines each before reading the E wind component again
          inDataWindE[,kk] <- scan(inFile, skip=107*2, what=numeric(), nmax=lm)
          inDataWindN[,kk] <- scan(inFile, what=numeric(), nmax=lm)
          inDataSLP[,kk] <- scan(inFile, what=numeric(), nmax=lm)
        }
      }
# Calculate the mean SLP 
      SLPDailyMean[,ll] <- (inDataSLP[,1] + inDataSLP[,2] + 
        inDataSLP[,3] + inDataSLP[,4]) / 4
      windEDailyMean[,ll] <- (inDataWindE[,1] + inDataWindE[,2] +
        inDataWindE[,3] + inDataWindE[,4]) / 4
      windNDailyMean[,ll] <- (inDataWindN[,1] + inDataWindN[,2] + 
        inDataWindN[,3] + inDataWindN[,4]) / 4
    }

    writeBin(as.vector(SLPDailyMean), outSLPDailyMean)
    writeBin(as.vector(windEDailyMean), outWindEDailyMean)
    writeBin(as.vector(windNDailyMean), outWindNDailyMean)

    SLPMonthlySum <- vector(mode="numeric", length=lm)
    windEMonthlySum <- vector(mode="numeric", length=lm)
    windNMonthlySum <- vector(mode="numeric", length=lm)
    for (mm in 1:days){
      SLPMonthlySum <- SLPMonthlySum + SLPDailyMean[,mm]
      windEMonthlySum <- windEMonthlySum + windEDailyMean[,mm]
      windNMonthlySum <- windNMonthlySum + windNDailyMean[,mm]
    }

    SLPMonthlyMean[,ii] <- SLPMonthlySum/days
    windEMonthlyMean[,ii] <- windEMonthlySum/days
    windNMonthlyMean[,ii] <- windNMonthlySum/days

    close(inFile)
    close(outSLPDailyMean)
    close(outWindEDailyMean)
    close(outWindNDailyMean)
  }

  writeBin(as.vector(SLPMonthlyMean), outSLPMonthlyMean)
  writeBin(as.vector(windEMonthlyMean), outWindEMonthlyMean)
  writeBin(as.vector(windNMonthlyMean), outWindNMonthlyMean)

  close(outSLPMonthlyMean)
  close(outWindEMonthlyMean)
  close(outWindNMonthlyMean)

# Calculate global means
  slpMMs <- vector(mode="numeric", length=12)
  windEMMs <- vector(mode="numeric", length=12)
  windNMMs <- vector(mode="numeric", length=12)
  
  for (mm in 1:12){
    slpMMs[mm] <- mean(SLPMonthlyMean[,mm])
    windEMMs[mm] <- mean(windEMonthlyMean[,mm])
    windNMMs[mm] <- mean(windNMonthlyMean[,mm])
  }

  globalSLPMonthlyMean[((jj*12)-11):(jj*12)] <- slpMMs
  globalWindEMonthlyMean[((jj*12)-11):(jj*12)] <- windEMMs
  globalWindNMonthlyMean[((jj*12)-11):(jj*12)] <- windNMMs

}

# Produce some nice monthly mean plots of the 1st month
SLP <- SLPMonthlyMean[,1]
dim(SLP) <- c(41,26)
image.plot(lon, lat, SLP, col=jet.colors(12), main="Mean SLP January 1999",
  xlab="Longitude", ylab="Latitude")

windE <- windEMonthlyMean[,1]
dim(windE) <- c(41,26)
x11()
image.plot(lon, lat, windE, col=jet.colors(12), 
  main="Mean E Wind January 1999", xlab="Longitude", ylab="Latitude")

windN <- windNMonthlyMean[,1]
dim(windN) <- c(41,26)
x11()
image.plot(lon, lat, windN, col=jet.colors(12), 
  main="Mean N Wind January 1999", xlab="Longitude", ylab="Latitude")
