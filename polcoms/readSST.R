###############################################################################
# Variation of readzet.r which calculates the 2D daily mean SSTs from the
# POLCOMMS model output. 
# 
# The SST arrays are output twice a day from the model at 2am and 2pm
#
# Simon Holgate, March 2007
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
# intersect(which(isea==86),which(jsea==90)) => ip=17701
#> (-5.55-xmin)/xres
#[1] 85.69998
#> (50.1-ymin)/yres
#[1] 89.90001

library(fields)
# Use jet colors for images
source('~/bin/RScripts/jet.colors.R')
lat<-seq(from=40.11111,by=1/9,length=224)
lon<-seq(from=-19.83333,by=1/6,length=198)

l<-198
m<-224
lm<-l*m
n<-66

# Number of months to iterate over this will be determined by the number of
# file names in inname
inname<-c( 
"SST.S12run405.jan",
"SST.S12run405.feb",
"SST.S12run405.mar",
"SST.S12run405.apr",
"SST.S12run405.may",
"SST.S12run405.jun",
"SST.S12run405.jul",
"SST.S12run405.aug",
"SST.S12run405.sep",
"SST.S12run405.oct",
"SST.S12run405.nov",
"SST.S12run405.dec" 
)

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

numYears<-length(yearsArray)
newlynMonthlyMean<-array(NA,dim=c(12,numYears))

globalSSTMonthlyMean <- vector(mode="numeric", length=(numYears*12))

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

  nmons<-length(inname)

  SSTMonthlyMean<-array(NA,dim=c(lm,12))
  SSTDailyMean<-array(NA,dim=c(lm,31))

  outMonthlyMean <- file(paste("monthlyMeans/sst.mm.",
    yearsArray[jj],".dat",sep=""), "wb")

  for (ii in 1:nmons) {

    fname <- paste(inname[ii], substr(as.character(yearsArray[jj]),3,4),sep="")
    inFile <- file(fname, "r")

    outDailyMean <- file(paste("dailyMeans/sst.dm.",monthArray[ii],".",   
      yearsArray[jj],".dat",sep=""), "wb")

    if ( leapFlag == 1) {
  # Leap year
      days<-leapYrDaysArray[ii]
    } else {
  # Regular year
      days<-daysArray[ii]
    }

# Read data into a table, ignoring the first line
    for (ll in 1:days) {
      inData <- array(NA, dim=c(lm,2))
      for (kk in 1:2) {
        inData[,kk] <- scan(inFile, skip=1, what=integer(), nmax=lm)
      }
# Calculate the mean SST and divide by 1000 to put into Celsius
      SSTDailyMean[,ll] <- (inData[,1]+inData[,2])/2/1000
    }

    writeBin(as.vector(SSTDailyMean), outDailyMean)

    monthlySum <- vector(mode="numeric", length=lm)
    for (mm in 1:days){
      monthlySum <- monthlySum + SSTDailyMean[,mm]
    }

    SSTMonthlyMean[,ii] <- monthlySum/days

    close(inFile)
    close(outDailyMean)
  }

  writeBin(as.vector(SSTMonthlyMean), outMonthlyMean)
  close(outMonthlyMean)

  sstMMs <- vector(mode="numeric", length=12)
  for (mm in 1:12){
    sstMMs[mm] <- mean(SSTMonthlyMean[,mm])
  }

  globalSSTMonthlyMean[((jj*12)-11):(jj*12)] <- sstMMs

}
